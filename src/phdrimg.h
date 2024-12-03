/*
 *  phdrimg.h
 *  panlib
 *
 *  Classes for building high dynamic-range images from multiple exposures.
 *
 *  Include after pancine.h
 *
 *  Created by gward on Fri Sep 07 2001.
 *  Copyright (c) 2006 Anyhere Software. All rights reserved.
 *
 */

#ifndef _PHDRIMG_H_
#define _PHDRIMG_H_

#define P_MAXHDR	50		// maximum number of exposures
#define P_MAXEXP	5		// max. poly exponent
#define P_ROTHRESH	0.07f		// rotation threshold (degrees)

struct PPatchSet;			// a set of patch exposure values

// A simple polynomial evaluation class
struct PPolynomial {
	int		degree;				// polynomial order
	double		coef[P_MAXEXP+1];		// coefficients
			// Evaluate polynomial at x
	double		Eval(double x) const {
				int	i = degree;
				double	res = coef[i];
				while (i-- > 0) res = x*res + coef[i];
				return res;
			}
			// Evaluate derivative at x
	double		Deriv(double x) const {
				int	i = degree;
				double	res = i*coef[i];
				while (--i > 0) res = x*res + i*coef[i];
				return res;
			}
			// Get float array field
	bool		Get(DBField *fp) const;
			// Assign from float array field
	bool		Set(const DBField &f);
};

// Assign polynomial from a float array field
inline bool
PPolynomial::Set(const DBField &f)
{
	float	coeff[P_MAXEXP+1];
	int	nc = f.Get(coeff, P_MAXEXP+1);
	degree = -1;
	if ((nc <= 0) | (f.GetNV() > P_MAXEXP+1))
		return false;
	degree = nc - 1;
	for (int i = nc; i--; )
		coef[i] = coeff[i];
	return true;
}

// Assign a float array field from a polynomial
inline bool
PPolynomial::Get(DBField *fp) const
{
	float	coeff[P_MAXEXP+1];
	if (degree < 0)
		return false;
	for (int i = degree; i >= 0; i--)
		coeff[i] = coef[i];
	fp->Set(coeff, degree + 1);
	return true;
}

// Read polynomial from input stream
inline istream &
operator>>(istream &ci, PPolynomial &pl)
{
	int	i;
	pl.degree = -1;
	ci >> pl.degree;
	if ((pl.degree < 0) | (pl.degree > P_MAXEXP)) {
		ci.clear(ios::failbit);
		return ci;
	}
	for (i = pl.degree; i >= 0; i--) {
		ci >> pl.coef[i];
		if (ci.fail())
			break;
	}
	return ci;
}

// Write polynomial to output stream
inline ostream &
operator<<(ostream &co, const PPolynomial &pl)
{
	co << pl.degree;
	for (int i = pl.degree; i >= 0; i--)
		co << ' ' << pl.coef[i];
	co << endl;
	return co;
}

// normalize uby8 value [0,256) to [0,1)
#define PbyteVal(v)	((1./256.)*(double(v)+.5))

// Structure for low dynamic-range image
struct PHDExposure {
	PCacheImage *	pci;			// cache image holder
	int		xoff, yoff;		// pixel offset from nominal image
	float		degCW;			// rotation from nominal image
	ABitMap2	usepix;			// map of usable pixels
	float		stonits;		// original sample-to-nits value
	ImgStruct	imr;			// resident image
	int		refCount;		// image reference count
	int		minv[3];		// mimimum usable values
	int		pct17[3];		// 17th percentile values
	int		median[3];		// 50th percentile values
	int		pct83[3];		// 83rd percentile values
	int		maxv[3];		// maximum usable values
			PHDExposure() {
				stonits = 0;
				xoff = yoff = 0; degCW = 0;
				refCount = 0;
				imr.img = NULL; imr.mbase = NULL;
				imr.csp = NULL;
				minv[0] = minv[1] = minv[2] =
				maxv[0] = maxv[1] = maxv[2] = -1;
				pci = new PCacheImage;
			}
			~PHDExposure() {
				PfreeImage(&imr);
				pci->HolderRelease();
			}
			// Check that image is resident
	bool		GetReady() {
				if (refCount++ > 0) return true;
				if (!pci->Ready()) return false;
				ImgReader *	ir = pci->GetReader();	
				imr.csp = &ir->cs;
				imr.xres = imr.yres = 0;
				return pci->GetCropped(&imr);
			}
			// Release resident image for now
	void		Release() {
				if (refCount <= 0) return;
				if (!--refCount) PfreeImage(&imr);
			}
			// Release resources and close image
	void		Done() {
				PfreeImage(&imr);
				pci->CloseImage();
				usepix.NewBitMap(0,0);
				xoff = yoff = 0; refCount = 0;
				minv[0] = minv[1] = minv[2] =
				minv[0] = minv[1] = minv[2] = -1;
			}
			// Check to see if pixel (x,y) is usable
	bool		Usable(int x, int y) const {
				return usepix.Check(x,y);
			}
			// Compute usable value range and map (call first)
	bool		ComputeUsability();
			// Compute position offset for alignment (call second)
	bool		ComputeOffset(PHDExposure *that);
			// Add to exposure patch list for this exposure
	int		Add2Patches(PPatchSet *ps);
			// Check if 24-bit RGB value is out of range
	bool		UnderRange(const uby8 *pp) const {
				return( (pp[0] < minv[0]) |
					(pp[1] < minv[1]) |
					(pp[2] < minv[2]) );
			}
	bool		OverRange(const uby8 *pp) const {
				return(	(pp[0] > maxv[0]) |
					(pp[1] > maxv[1]) |
					(pp[2] > maxv[2]) );
			}
private:
			// Average rectangular patch
	bool		AvgPatch(uby8 rgb[3], const ImgRect &r) const;
			// Find a new, light patch
	bool		NewLightPatch(uby8 rgb[3], ImgRect *r,
					const PPatchSet *ps, int lvl) const;
};

// Polynomial response functions and ordered exposure adjustments
struct PHDSoln {
	PPolynomial	ply[3];			// response functions
	float		expadj[P_MAXHDR];	// exposure adjustments
	void		Init();
};

// Flags for our solution
enum PHDFlags {
		PHDFresponse = 0x1,		// solve for responses
		PHDFexposure = 0x2,		// solve for exposures
		PHDFalignment = 0x4,		// solve for alignment
		PHDFall = 0x7,			// solve for everything
		PHDFghosting = 0x8,		// remove ghosts after
		PHDFskipexp = 0x10,		// skip unneeded exposures
};

// Info for list of camera response functions
enum CDBFieldID {
		CDBFmake,			// camera make
		CDBFmodel,			// camera model
		CDBFversion,			// camera version
		CDBFred,			// red response function
		CDBFgreen,			// green response function
		CDBFblue,			// blue response function
		CDBFend				// terminator
};

// Camera database field information object
class CDBFieldInfo : public DBFieldInfo {
public:
			CDBFieldInfo();
};
extern const CDBFieldInfo	CDBFInfo;		// only need one

// Class for computing high dynamic-range image
class PHDImageMaker {
	PHDExposure	ipi[P_MAXHDR];		// input images
	int		nipi;			// number of input images
	int		ord[P_MAXHDR];		// exposure order: long to short
	int		ustart, uend;		// useful exposure limit
	PHDSoln		soln;			// responses and exposures
	int		haveSoln;		// what have we solved for?
	PHDExposure &	Inp(int i) {
				return ipi[ord[i]];
			}
	const PHDExposure &
			Inp(int i) const {
				return ipi[ord[i]];
			}
	double		MaxExpAdj(int i1, int i0) const;
	double		CompErr(const PHDSoln *sl, const PPatchSet *ps) const;
	bool		FindPoly(PHDSoln *sl, const int n,
					const PPatchSet *ps) const;
	bool		FindExposures(PHDSoln *sl, const PPatchSet *ps) const;
	bool		FindSoln(PHDSoln *sl, const int n,
					const PPatchSet *ps) const;
	bool		SortExposures();
	bool		SkipExposure(int n);
	int		GetBestExposure(ImgStruct *ibp = NULL);
public:
	int		solveFlags;		// what we're solving for
						// progress report function
	int		(*reportProgress)(const char *, int);
			PHDImageMaker() {
				nipi = 0; ord[0] = -1;
				ustart = uend = -1;
				solveFlags = PHDFall; haveSoln = 0;
				reportProgress = NULL;
			}
			// Allocate next cache image (call first!)
	PCacheImage *	NewCacheImage() {
				if (nipi >= P_MAXHDR) return NULL;
				ClearSoln();
				return ipi[nipi++].pci;
			}
			// Get pointer to ith cache image (unordered)
	const PCacheImage *
			GetCacheImage(int i) const {
				if ((i < 0) | (i >= nipi)) return NULL;
				return ipi[i].pci;
			}
			// Get number of exposures
	int		GetNI() const {
				return nipi;
			}
			// What remains to be solved?
	int		Unsolved() const {
				return (solveFlags & ~haveSoln & PHDFall);
			}
			// Compute our solution
	bool		Compute();
			// Clear our solution (may clear assigned responses!)
	void		ClearSoln(int clrWhat = PHDFall);
			// Write out polynomial response functions
	bool		WriteResponse(ostream *ostr) {
				if (!Compute()) return false;
				*ostr << soln.ply[0];
				*ostr << soln.ply[1];
				*ostr << soln.ply[2];
				return !ostr->bad();
			}
	bool		WriteResponse(const char *fname) {
				ofstream	ofs(fname, ios::out|ios::trunc);
				if (!ofs.is_open()) return false;
				return WriteResponse(&ofs);
			}
			// Read and use polynomial response functions
	bool		ReadResponse(istream *istr) {
				ClearSoln(PHDFexposure|PHDFresponse);
				*istr >> soln.ply[0];
				*istr >> soln.ply[1];
				*istr >> soln.ply[2];
				if (istr->fail()) return false;
				haveSoln |= PHDFresponse;
				solveFlags &= ~PHDFresponse;
				return true;
			}
	bool		ReadResponse(const char *fname) {
				ifstream	ifs(fname);
				if (!ifs.is_open()) return false;
				return ReadResponse(&ifs);
			}
			// Save response function to camera response record
	bool		GetResponse(DBRecord *crp);
			// Set response function from camera response record
	bool		SetResponse(const DBRecord &cr);
			// Get polynomial for indicated primary
	const PPolynomial *
			GetPolynomial(int pri) {
				if ((pri < 0) | (pri >= 3)) return NULL;
				if (!Compute()) return NULL;
				return &soln.ply[pri];
			}
			// Evaluate response curve for indicated primary
	bool		GetResponse(float resp[256], int pri, float *der=NULL);
			// Get corrected stonits for ith ordered exposure
	float		GetStoNits(int i = -1) {
				if (!Compute()) return -1.f;
				if (i < 0) i = nipi/2;
				return Inp(i).stonits*soln.expadj[i];
			}
			// Compute and return high dynamic-range image
	bool		GetHDImage(ImgStruct *ims, ImgStruct *ivs = NULL);
			// Close images and free data
	void		Done() {
				while (nipi > 0) ipi[--nipi].Done();
			}
};

// Compute and remove lens flare from a high dynamic range image
extern bool	PHDremoveFlare(ImgStruct *ims, int (*pbar)(const char *, int) = NULL);

#endif	// ! _PHDRIMG_H_
