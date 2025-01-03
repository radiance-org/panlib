/*
 *  phdrimg.cpp
 *  panlib
 *
 *  Implementation of class to compute high dynamic-range image.
 *
 *  Based on Mitsunaga and Nayar paper "Radiometric Self Calibration."
 *  Limits are placed on exposure adjustments, since
 *  this part of the algorithm seems to add to convergence problems.
 *
 *  Created by Greg Ward on Fri Sep 07 2001.
 *  Copyright (c) 2006 Anyhere Software. All rights reserved.
 *
 */

#include <math.h>
#include "pancine.h"
#include "phdrimg.h"
#include "gaussjord.h"
#ifdef phd_debug
#include "imgwriter.h"
#endif

#define P_MAXNOISE	30			// maximum expected noise value
#define P_LINRESP	19			// response function lower limit
#define P_SAFEMAX	232			// safe maximum value

#define P_MAXERR	0.3			// acceptable RMS sol'n error
#define P_PTARGET	50			// target # patches
#define P_MAXPATCH	(8*P_PTARGET)		// max. # patches
#define P_PWIDTH	12			// patch width (pixels)
#define P_PVARIANCE	0.5			// max. patch variance
// Applying a saturation limit just seems to make matters worse...
// #define P_PSATURATION   0.4			// max. patch saturation

// Pixel variance measure equal to variance over mean (hey, it works!)
#define p_variance(sum,sum2,n) 	((double)(sum2)/(sum) - (double)(sum)/(n))

static const int		safe_minmax[2] = {P_MAXNOISE, P_SAFEMAX};
static const int		safe_margin[2] = {2, -16};
static const float		safe_tail[2] = {.002f, .0005f};

static const PPolynomial	linearPoly1 = {1, {0., 1.}};
						// RGB luminance coef.
static const double		ccoef[3] = {.256, .678, .066};

#define P_WTMIN			1e-4f		// minimum acceptable weight

// A set of patch exposure values
struct PPatchSet {
	struct PPatch {
		uby8		avg[3];				// averaged RGB
		uby8		valid;				// usable
	}		pval[P_MAXPATCH][P_MAXHDR];	// patch values
	ImgRect		ppos[P_MAXPATCH];		// patch positions
	int		pcount;				// patch count
	int		nexp;				// number of exposures
			PPatchSet() {
				pcount = 0; nexp = 0;
			}
};

// Compare two polynomials to see if they are within epsilon over [0-255]
static bool
samePoly(const PPolynomial &p1, const PPolynomial &p2, const double eps = .005)
{
	if (p1.degree != p2.degree)
		return false;
	for (int i = 256; i--; ) {
		double	v1 = p1.Eval(PbyteVal(i));
		if (absval(p2.Eval(PbyteVal(i)) - v1) > eps*(v1 > eps ? v1 : eps))
			return false;
	}
	return true;
}

// Crude brightness estimate for 24-bit RGB value (0-255 range)
static inline int
rgb_bright(const uby8 rgb[3])
{
	return ( 54*rgb[0] + 183*rgb[1] + 19*rgb[2] ) >> 8;
}

// Crude saturation estimate for 24-bit RGB value
static inline double
rgb_saturation(const uby8 rgb[3])
{
	int     vmin = (rgb[0] < rgb[1]) ? rgb[0] : rgb[1];
	if (rgb[2] < vmin)
		vmin = rgb[2];
	if (!vmin)
		return 1.;
	return 1. - (double)vmin/rgb_bright(rgb);
}

// Check if two patches overlap
static inline bool
patchOverlap(const ImgRect &r1, const ImgRect &r2)
{
	if (r1.xright <= r2.xleft) return false;
	if (r1.ybottom <= r2.ytop) return false;
	if (r2.xright <= r1.xleft) return false;
	if (r2.ybottom <= r1.ytop) return false;
	return true;
}

const CDBFieldInfo	CDBFInfo;		// Camera response info

// Initialize Camera field information (order must match enum CDBFieldID)
CDBFieldInfo::CDBFieldInfo()
{
	static const char	eemsg[] = "enumeration error in CDBFieldInfo";

	if (Field("Make", DBDTstring, "Camera manufacturer",
			0) != CDBFmake)
		DMESG(DMCassert, eemsg);
	if (Field("Model", DBDTstring, "Camera model",
			0) != CDBFmodel)
		DMESG(DMCassert, eemsg);
	if (Field("Version", DBDTstring, "Camera version",
			0) != CDBFversion)
		DMESG(DMCassert, eemsg);
	if (Field("Red", DBDTfloat, "Red response coefficients",
			DBFFarray) != CDBFred)
		DMESG(DMCassert, eemsg);
	if (Field("Green", DBDTfloat, "Green response coefficients",
			DBFFarray) != CDBFgreen)
		DMESG(DMCassert, eemsg);
	if (Field("Blue", DBDTfloat, "Blue response coefficients",
			DBFFarray) != CDBFblue)
		DMESG(DMCassert, eemsg);
}

// Compute value range for each primary and usability map
bool
PHDExposure::ComputeUsability()
{
	int			i, j;
						// check if we're done already
	if (minv[0] < maxv[0])
		return true;
	if (!PmatchColorSpace(&pci->GetReader()->cs, &ICS_sRGB, PICMptype)) {
		DMESG(DMCparameter, "HDR image exposures must be 24-bit RGB");
		return false;
	}
	if (!GetReady()) {
		DMESG(DMCresource, "Cannot load input exposure");
		return false;
	}
						// compute RGB histogram
	uby8			extrem8[3][2] = {{0,255}, {0,255}, {0,255}};
	unsigned long		hcount[3][256];
	memset(hcount, 0, sizeof(hcount));
	PcomputeHisto(extrem8, hcount[0], 256, &imr, PHCrgb3);
						// accumulate safe tails
	for (j = 3; j--; ) {
		long	cnt;
		cnt = long(safe_tail[0]*imr.xres*imr.yres + .5f);
		for (i = 0; i < safe_minmax[0]; ++i)
			if ((cnt -= hcount[j][i]) <= 0) break;
		minv[j] = i + safe_margin[0];
		cnt = long(0.01f*safe_tail[0]*imr.xres*imr.yres + 1.f);
		for (i = 0; i < 256; i++)	// don't pass 1% of safe tail
			if ((cnt -= hcount[j][i]) <= 0) break;
		if (i > minv[j]) minv[j] = i;
		cnt = long(safe_tail[1]*imr.xres*imr.yres + .5f);
		for (i = 255; i > safe_minmax[1]; --i)
			if ((cnt -= hcount[j][i]) <= 0) break;
		maxv[j] = i + safe_margin[1];
		cnt = long(0.01f*safe_tail[1]*imr.xres*imr.yres + 1.f);
		for (i = 256; i--; )		// don't pass 1% of safe tail
			if ((cnt -= hcount[j][i]) <= 0) break;
		if (i < maxv[j]) maxv[j] = i;
						// compute median value
		cnt = (unsigned long)imr.xres*imr.yres/2;
		for (i = 255; cnt > 0; --i)
			cnt -= hcount[j][i];
		median[j] = ++i;
						// compute 17th percentile
		cnt = (unsigned long)imr.xres*imr.yres/6;
		for (i = 0; cnt > 0; i++)
			cnt -= hcount[j][i];
		pct17[j] = i;
						// compute 83rd percentile
		cnt = (unsigned long)imr.xres*imr.yres/6;
		for (i = 255; cnt > 0; i--)
			cnt -= hcount[j][i];
		pct83[j] = ++i;
	}
#ifdef phd_debug
	fprintf(stderr, "Restricted range for RGB is [%d,%d] [%d,%d] [%d,%d]\n",
			minv[0], maxv[0], minv[1], maxv[1],
			minv[2], maxv[2]);
#endif
						// compute usability map
	usepix.NewBitMap(imr.xres, imr.yres, true);
	for (j = imr.yres; j--; ) {
		const uby8 *	pp = ProwPtr(&imr, j);
		for (i = 0; i < imr.xres; ++i, pp += 3)
			if (OverRange(pp) || UnderRange(pp))
				usepix.Reset(i,j);
	}
	Release();				// all done
	return true;
}

// Average a patch if we can
bool
PHDExposure::AvgPatch(uby8 rgb[3], const ImgRect &r) const
{
	ImgRect	myrect = r;
	if ((myrect.xleft += xoff) < 0)
		myrect.xleft = 0;
	if ((myrect.xright += xoff) > imr.xres)
		myrect.xright = imr.xres;
	if ((myrect.ytop += yoff) < 0)
		myrect.ytop = 0;
	if ((myrect.ybottom += yoff) > imr.yres)
		myrect.ybottom = imr.yres;
	const int	cnt = (myrect.xright - myrect.xleft) *
				(myrect.ybottom - myrect.ytop);
	if (cnt <= P_PWIDTH*P_PWIDTH/2)
		return false;
	int	rgbsum[3] = {0, 0, 0};
	for (int y = myrect.ytop; y < myrect.ybottom; y++) {
		int			x = myrect.xleft;
		const uby8 *	pp = PpixP(&imr, x, y, 3);
		for ( ; x < myrect.xright; x++, pp += 3) {
			if (!Usable(x,y))
				return false;
			rgbsum[0] += pp[0];
			rgbsum[1] += pp[1];
			rgbsum[2] += pp[2];
		}
	}
	rgb[0] = (rgbsum[0] + (cnt>>1)) / cnt;
	rgb[1] = (rgbsum[1] + (cnt>>1)) / cnt;
	rgb[2] = (rgbsum[2] + (cnt>>1)) / cnt;
	return true;
}

// Find a new, light patch
bool
PHDExposure::NewLightPatch(uby8 rgb[3], ImgRect *r,
				const PPatchSet *ps, int lvl) const
{
	if (lvl >= maxv[1])
		return false;
	int	n_iter = (imr.xres/P_PWIDTH - 1)*(imr.yres/P_PWIDTH - 1);
	while (n_iter--) {
		// should we restrict search to center to avoid vignetting?
		// irrelevant if all exposures at same aperture as recommmended
		const int	xleft = rand() % (imr.xres - P_PWIDTH);
		const int	ytop = rand() % (imr.yres - P_PWIDTH);
		long		sum[3] = {0, 0, 0};
		long		sum2[3] = {0, 0, 0};
		int		j;
		int		i;
		r->xright = (r->xleft = xleft - xoff) + P_PWIDTH;
		r->ybottom = (r->ytop = ytop - yoff) + P_PWIDTH;
		for (i = ps->pcount; i--; )
			if (patchOverlap(*r, ps->ppos[i]))
				goto nextpatch;	// overlaps existing patch
		for (j = ytop + P_PWIDTH; j-- > ytop; ) {
			const uby8 *	pp = PpixP(&imr, xleft, j, 3);
			for (i = xleft; i < xleft + P_PWIDTH; i++, pp += 3) {
				int	v;
				if (!Usable(i,j))
					goto nextpatch;
				v = pp[0]; sum[0] += v; sum2[0] += v*v;
				v = pp[1]; sum[1] += v; sum2[1] += v*v;
				v = pp[2]; sum[2] += v; sum2[2] += v*v;
			}
		}
		if ((sum[0] <= 0) | (sum[1] <= 0) | (sum[2] <= 0))
			continue;
		rgb[0] = (sum[0] + P_PWIDTH*P_PWIDTH/2) / (P_PWIDTH*P_PWIDTH);
		rgb[1] = (sum[1] + P_PWIDTH*P_PWIDTH/2) / (P_PWIDTH*P_PWIDTH);
		rgb[2] = (sum[2] + P_PWIDTH*P_PWIDTH/2) / (P_PWIDTH*P_PWIDTH);
					// bias towards brighter patches
		if (rgb_bright(rgb) - lvl < rand()%(maxv[1]-lvl))
			continue;
#ifdef P_PSATURATION
					// bias towards neutral patches
		if (rgb_saturation(rgb) > P_PSATURATION/RAND_MAX*rand())
			continue;
#endif
#ifdef P_PVARIANCE
					// must meet variance requirements
		if (p_variance(sum[0], sum2[0], P_PWIDTH*P_PWIDTH) > P_PVARIANCE)
			continue;
		if (p_variance(sum[1], sum2[1], P_PWIDTH*P_PWIDTH) > P_PVARIANCE)
			continue;
		if (p_variance(sum[2], sum2[2], P_PWIDTH*P_PWIDTH) > P_PVARIANCE)
			continue;
#endif
		return true;		// this patch is OK
	nextpatch:;
	}
	return false;			// didn't find one after max. iterations
}

// Add new exposure to patch list
int
PHDExposure::Add2Patches(PPatchSet *ps)
{
	int	nnew = 0;
	int	nappended = 0;
	int	i;
					// get image
	if (!GetReady())
		return 0;
					// append to existing patches
	int	lightlvl = 0;
	for (i = 0; i < ps->pcount; i++)
		if ((ps->pval[i][ps->nexp].valid =
				AvgPatch(ps->pval[i][ps->nexp].avg,
						ps->ppos[i]))) {
			int	v = rgb_bright(ps->pval[i][ps->nexp].avg);
			if (v > lightlvl)
				lightlvl = v;
			++nappended;
		}
					// add new patches
	int	ptarget = P_PTARGET - nappended;
	if (nappended > P_PTARGET*lightlvl/256)
		ptarget = P_PTARGET - P_PTARGET*lightlvl/256;
	while ((ps->pcount < P_MAXPATCH) & (nnew < ptarget)) {
		if (!NewLightPatch(ps->pval[ps->pcount][ps->nexp].avg,
				&ps->ppos[ps->pcount], ps, lightlvl)) {
			DMESG(DMCwarning, "Trouble finding HDR patches");
			break;
		}
		ps->pval[ps->pcount][i=ps->nexp].valid = true;
		while (i--)		// clear other exposures
			ps->pval[ps->pcount][i].valid = false;
		ps->pcount++;		// added new patch
		++nnew;
	}
	Release();			// done with image
	ps->nexp++;			// exposure's been added
	return nnew;
}

// Find least squares polynomial solution for patch set
bool
PHDImageMaker::FindPoly(PHDSoln *sl, const int n, const PPatchSet *ps) const
{
	int	q0, q1, p, i, j, k;
	DASSERT(uend-ustart == ps->nexp);
					// initialize
	sl->ply[0].degree = sl->ply[1].degree = sl->ply[2].degree = n;
					// find polynomial for each primary
	for (i = 3; i--; ) {
		double	mat[P_MAXEXP*P_MAXEXP];
		double	rhs[P_MAXEXP];
					// compute matrix coefficients
		for (j = n*n; j--; )
			mat[j] = 0.;
		for (j = n; j--; )
			rhs[j] = 0.;
		for (q1 = uend; --q1 > ustart; )
		    for (q0 = q1; q0-- > ustart; ) {
			const double	erat = sl->expadj[q1] *
						Inp(q1).stonits /
						(sl->expadj[q0] *
						Inp(q0).stonits);
			for (p = ps->pcount; p--; ) {
				if (!ps->pval[p][q1].valid |
						!ps->pval[p][q0].valid)
					continue;
				const double	m0 =
					PbyteVal(ps->pval[p][q0].avg[i]);
				const double	m1 =
					PbyteVal(ps->pval[p][q1].avg[i]);
				double	e0 = 1.;
				double	e1 = 1.;
				double	dvec[P_MAXEXP+1];
				for (k = 0; k < n; k++) {
					dvec[k] = e0 - erat*e1;
					e0 *= m0; e1 *= m1;
				}
				dvec[n] = e0 - erat*e1;
				for (j = 0; j < n; j++)
					for (k = 0; k < n; k++)
						mat[j*n+k] += dvec[j] *
							(dvec[k]-dvec[n]);
				for (j = 0; j < n; j++)
					rhs[j] -= dvec[j]*dvec[n];
			}
		    }
					// solve linear system
		if (!GJsolveLinearEq(n, mat, rhs, sl->ply[i].coef))
			return false;
		sl->ply[i].coef[n] = 1.;
		for (j = n; j--; )
			sl->ply[i].coef[n] -= sl->ply[i].coef[j];
	}
	return true;
}

// Recompute exposure adjustments based on patch set
bool
PHDImageMaker::FindExposures(PHDSoln *sl, const PPatchSet *ps) const
{
	float	erat[P_MAXHDR];
	int	q0, q1, p, i;
	DASSERT(uend-ustart == ps->nexp);
	for (q1 = uend; --q1 > ustart; ) {
		double	sum = 0.;
		int	cnt = 0;
		q0 = q1 - 1;		// compute ratio between q and q-1
		for (p = ps->pcount; p--; ) {
			if (!ps->pval[p][q1].valid |
					!ps->pval[p][q0].valid)
				continue;
			for (i = 3; i--; ) {
				sum += sl->ply[i].Eval(
					PbyteVal(ps->pval[p][q1].avg[i]))
					/ sl->ply[i].Eval(
					PbyteVal(ps->pval[p][q0].avg[i]));
				++cnt;
			}
		}
		if (cnt > 0)
			erat[q1] = sum / (double)cnt;
		else
			erat[q1] = sl->expadj[q0]*Inp(q0).stonits /
					(sl->expadj[q1]*Inp(q1).stonits);
	}
					// get final adjustments from ratios
	for (q1 = ustart+1; q1 < uend; q1++) {
		double	maxadj = MaxExpAdj(ord[q1], ord[q0]);
		q0 = q1 - 1;
		sl->expadj[q1] = sl->expadj[q0]*Inp(q0).stonits /
					(erat[q1]*Inp(q1).stonits);
		if (sl->expadj[q1] > maxadj)
			sl->expadj[q1] = maxadj;
		else if (sl->expadj[q1] < 1./maxadj)
			sl->expadj[q1] = 1./maxadj;
	}
	return true;
}

// Initialize solution
void
PHDSoln::Init()
{
	int	i;
	for (i = 3; i--; )
		ply[i] = linearPoly1;
	for (i = P_MAXHDR; i--; )
		expadj[i] = 1.f;
}

// Find nth degree polynomial solution and relative exposure adjustments
bool
PHDImageMaker::FindSoln(PHDSoln *sl, const int n, const PPatchSet *ps) const
{
	int	niter = 50;
					// initialize
	sl->Init();
					// find optimal solution set
	PPolynomial	lastPly[3];
	do {
		lastPly[0] = sl->ply[0];
		lastPly[1] = sl->ply[1];
		lastPly[2] = sl->ply[2];
					// find polynomial for each primary
		if (!FindPoly(sl, n, ps))
			return false;
					// recalibrate exposure adjustments
		if (!FindExposures(sl, ps))
			return false;
		if (!niter--) {
			DMESGF(DMCwarning, "Poor convergence for order %d fit", n);
			break;
		}
	} while ( !( samePoly(sl->ply[0],lastPly[0]) &&
			samePoly(sl->ply[1],lastPly[1]) &&
			samePoly(sl->ply[2],lastPly[2]) ) );
	return true;
}

// Compute maximum exposure adjustment between exposure i1 and i0
double
PHDImageMaker::MaxExpAdj(int i1, int i0) const
{
	if (i1 == i0)
		return 1.;
	double	ma = 1.12;
	float	aper1, aper0;
	float	time1, time0;
	if (!PDBgetField(ipi[i1].pci->ircd,PDBFexptime).Get(&time1))
		return ma;
	if (!PDBgetField(ipi[i0].pci->ircd,PDBFexptime).Get(&time0))
		return ma;
	if (!PDBgetField(ipi[i1].pci->ircd,PDBFaperture).Get(&aper1))
		return ma;
	if (!PDBgetField(ipi[i0].pci->ircd,PDBFaperture).Get(&aper0))
		return ma;
	ma = 1.04;
	if (time1 != time0)
		ma += .04;
	if (aper1 != aper0)
		ma += .25;
	return ma;
}

// Compute relative RMS error associated with the given solution
double
PHDImageMaker::CompErr(const PHDSoln *sl, const PPatchSet *ps) const
{
	double	err = 0.;
	int	cnt = 0;
	int	q, p;
	DASSERT(uend-ustart == ps->nexp);
	for (q = uend; --q > ustart; ) {
		const double	s2n1 = sl->expadj[q]*Inp(q).stonits;
		const double	s2n0 = sl->expadj[q-1]*Inp(q-1).stonits;
		for (p = ps->pcount; p--; ) {
			if (!ps->pval[p][q].valid | !ps->pval[p][q-1].valid)
				continue;
			double	rgb1[3], rgb0[3];
			for (int i = 3; i--; ) {
				rgb1[i] = s2n1 * sl->ply[i].Eval(
					PbyteVal(ps->pval[p][q].avg[i]));
				rgb0[i] = s2n0 * sl->ply[i].Eval(
					PbyteVal(ps->pval[p][q-1].avg[i]));
			}
			double	rf = 2./( ccoef[0]*(rgb1[0]+rgb0[0]) +
						ccoef[1]*(rgb1[1]+rgb0[1]) +
						ccoef[2]*(rgb1[2]+rgb0[2]) );
			double	r = rf*(rgb1[0] - rgb0[0]);
			double	g = rf*(rgb1[1] - rgb0[1]);
			double	b = rf*(rgb1[2] - rgb0[2]);
			err += ccoef[0]*r*r + ccoef[1]*g*g + ccoef[2]*b*b;
			++cnt;
#if 0
			// Swapping greens in plot shows convergence
			fprintf(stderr, "\t%f\t%f\t",
					PbyteVal(ps->pval[p][q-1].avg[1]),
					rgb1[1]/s2n0);
			fprintf(stderr, "\t%f\t%f\n",
					PbyteVal(ps->pval[p][q].avg[1]),
					rgb0[1]/s2n1);
#endif
		}
	}
	return sqrt(err / (double)cnt);
}

// Get exposure values and order light to dark
bool
PHDImageMaker::SortExposures()
{
	if (ord[0] >= 0)
		return true;			// already sorted
	int	n;
	for (n = 0; n < nipi; n++) {
		if (!PDBgetField(ipi[n].pci->ircd,
				PDBFstonits).Get(&ipi[n].stonits)) {
			DMESGF(DMCdata, "%s: needs exposure calibration",
					ipi[n].pci->GetReader()->file);
			return false;
		}
		int	i;
		for (i = n; i--; ) {		// insertion sort
			if (ipi[n].stonits > ipi[ord[i]].stonits)
				break;
			ord[i+1] = ord[i];
		}
		ord[i+1] = n;
	}
	int	n2close = 0;
	for (n = nipi; --n > 0; )
		n2close += (ipi[ord[n]].stonits < 1.14f*ipi[ord[n-1]].stonits);
	DTESTF(n2close, DMCwarning, "%d exposure(s) redundant", n2close);
	return true;
}

// Did exposure capture the highlights in all three channels?
static inline bool
caughtHighlight(const PHDExposure &ipe)
{
	static const int	clipped_max = safe_minmax[1] + safe_margin[1];

	return (ipe.maxv[0] <= clipped_max) &
			(ipe.maxv[1] <= clipped_max) &
			(ipe.maxv[2] <= clipped_max);
}

// Did exposure capture the shadows in all three channels?
static inline bool
caughtShadow(const PHDExposure &ipe)
{
	static const int	clipped_min = safe_minmax[0] + safe_margin[0];

	return (ipe.minv[0] >= clipped_min) &
			(ipe.minv[1] >= clipped_min) &
			(ipe.minv[2] >= clipped_min);
}

// Shall we skip this exposure because it's redundant?
bool
PHDImageMaker::SkipExposure(int n)
{
	if (ustart < 0) {			// initialize
		ustart = 0; uend = nipi;
	}
	if ((ustart > n) | (n >= uend))
		return true;
	if (!(solveFlags & PHDFskipexp))
		return false;			// user wants them all
	if (!Inp(n).ComputeUsability())
		return true;			// shouldn't happen!
#if 0			// can't skip middle exposure -- breaks code
	if (n > 0 && Inp(n).stonits/Inp(n-1).stonits < 1.05f) {
		DMESGF(DMCtrace, "Exposure %d redundant", n);
		return true;
	}
#endif
	if (n > ustart && caughtHighlight(Inp(n))) {
		if (!Inp(n-1).ComputeUsability())
			return true;		// shouldn't happen, either
		if (caughtHighlight(Inp(n-1))) {
			DMESGF(DMCtrace, "Exposure %d too dark", n);
			uend = n;
			return true;
		}
	}
	if (n < uend-1 && caughtShadow(Inp(n))) {
		if (!Inp(n+1).ComputeUsability())
			return true;		// shouldn't happen, threether
		if (caughtShadow(Inp(n+1))) {
			DMESGF(DMCtrace, "Exposure %d too light", n);
			ustart = n+1;
			return true;
		}
	}
	return false;				// provides useful range
}

static int
intcmp(const void *p1, const void *p2)
{
	return *(const int *)p1 - *(const int *)p2;
}

static int
floatcmp(const void *p1, const void *p2)
{
	float	diff = *(const float *)p1 - *(const float *)p2;
	return (diff > 0 ? 1 : diff < 0 ? -1 : 0);
}

// Compute camera response function and alignment for multiple exposures
bool
PHDImageMaker::Compute()
{
	static const char	cancelMsg[] = "User cancelled HDR solver";
	if (!SortExposures())
		return false;
	if (!Unsolved())
		return true;			// nothing left to do
	if (nipi < 2 || ipi[ord[nipi-1]].stonits < 1.19f*ipi[ord[0]].stonits) {
						// nothing we can do
		if (Unsolved() & (PHDFresponse |
					(nipi>1)*(PHDFalignment|PHDFghosting))) {
			DMESG(DMCparameter,
				"Insufficient exposures to compute HDR image");
			return false;
		}
		haveSoln |= solveFlags;
		return true;
	}
	DMESG(DMCtrace, "Begin HDR solver...");
	if ((solveFlags & (PHDFresponse|PHDFexposure|PHDFghosting)) ==
				(PHDFresponse|PHDFexposure|PHDFghosting))
		DMESG(DMCwarning, "Finding exposure with movement is problematic");
						// get exposure patch values
	int		n;
	PPatchSet	patchList;
	Inp(0).xoff = Inp(0).yoff = 0;
	for (n = 0; n < nipi; n++) {
		if (reportProgress && !(*reportProgress)("Analyzing Images",
								100*n/nipi)) {
			DMESG(DMCtrace, cancelMsg);
			return false;
		}
		if (SkipExposure(n))
			continue;
		if (!Inp(n).ComputeUsability())
			return false;
		if (n > ustart && Unsolved() & PHDFalignment) {
			DMESG(DMCtrace, "Computing offset for image...");
			if (!Inp(n).ComputeOffset(&Inp(n-1)))
				return false;
		}
		if (!(Unsolved() & (PHDFresponse|PHDFexposure)))
			continue;
		DMESG(DMCtrace, "Finding patches in image...");
		Inp(n).Add2Patches(&patchList);
		DMESGF(DMCtrace, "%d patches total", patchList.pcount);
#ifdef phd_debug
	{				// write out image with patch locations
		int	i, j;
		char	fname[64];
		Inp(n).GetReady();
		ImgStruct	myi = Inp(n).imr;
		myi.img = NULL;
		PrenderImageI(&myi, &Inp(n).imr);
		Inp(n).Release();
		if ((Inp(n).xoff != 0) | (Inp(n).yoff != 0))
			PshiftImageB(&myi, -Inp(n).xoff,
					-Inp(n).yoff);
		for (i = patchList.pcount; i--; ) {
			if (!patchList.pval[i][n].valid) continue;
			uby8	*pp;
#if 1
			pp = PpixP(&myi, patchList.ppos[i].xleft,
					patchList.ppos[i].ytop, 3);
			for (j = 0; j < P_PWIDTH; j++, pp += 3) {
				pp[0] = pp[1] = pp[2] = 255;
				pp[P_PWIDTH*myi.rowsize+0] =
				pp[P_PWIDTH*myi.rowsize+1] =
				pp[P_PWIDTH*myi.rowsize+2] = 255;
			}
			pp = PpixP(&myi, patchList.ppos[i].xleft,
					patchList.ppos[i].ytop, 3);
			for (j = 0; j < P_PWIDTH; j++, pp += myi.rowsize) {
				pp[0] = pp[1] = pp[2] = 0;
				pp[P_PWIDTH*3+0] =
				pp[P_PWIDTH*3+1] =
				pp[P_PWIDTH*3+2] = 0;
			}
#else
			for (j = 0; j < P_PWIDTH; j++) {
				pp = PpixP(&myi, patchList.ppos[i].xleft,
						patchList.ppos[i].ytop+j, 3);
				for (int k = 0; k < P_PWIDTH; k++, pp += 3) {
					pp[0] = 255; pp[1] = 0; pp[2] = 0;
				}
			}
#endif
		}
		sprintf(fname, "patch%d.jpg", n);
		DMESGF(DMCtrace, "Writing image with patch locations to '%s'...",
				fname);
		ImgWriteBuf	iwb;
		PsetWriteBuf(&iwb, &myi);
		iwb.info.flags = IIFquality;
		iwb.info.quality = 100;
		(*IWInterfaceJPEG.WriteImage)(fname, &iwb);
		PfreeImage(&myi);
	}
#endif
	}
	if (reportProgress && !(*reportProgress)(NULL, 100)) {
		DMESG(DMCtrace, cancelMsg);
		return false;
	}
						// complete operations
	if (Unsolved() & PHDFalignment) {
		int	xomed, yomed;		// find median offset
		float	rotmed;			// and rotation
		int	ilist[P_MAXHDR];
		float	flist[P_MAXHDR];
		int	i;
		for (n = 0, i = ustart; i < uend; i++)
			ilist[n++] = Inp(i).xoff;
		qsort(ilist, n, sizeof(ilist[0]), intcmp);
		xomed = ilist[n/2];
		for (n = 0, i = ustart; i < uend; i++)
			ilist[n++] = Inp(i).yoff;
		qsort(ilist, n, sizeof(ilist[0]), intcmp);
		yomed = ilist[n/2];
		for (n = 0, i = ustart; i < uend; i++)
			flist[n++] = Inp(i).degCW;
		qsort(flist, n, sizeof(flist[0]), floatcmp);
		rotmed = flist[n/2];
		for (i = 0; i < nipi; i++) {	// nullify med. offset & rot.
			ipi[i].xoff -= xomed;
			ipi[i].yoff -= yomed;
			ipi[i].degCW -= rotmed;
		}
		haveSoln |= PHDFalignment;
		if (!Unsolved())
			return true;		// nothing left to do
	}
	if (Unsolved() == PHDFexposure) {       // solve exposures only
		if (!FindExposures(&soln, &patchList))
			return false;
		haveSoln |= PHDFexposure;
		return true;
	}
						// find best poly solution
	double		thisErr, bestErr = 1.;
	PHDSoln		thisSoln;
	soln.ply[0].degree = soln.ply[1].degree = soln.ply[2].degree = 0;
	for (n = nipi; n--; )
		thisSoln.expadj[n] = soln.expadj[n];
						// try each polynomial degree
	for (n = 1; n <= P_MAXEXP; n++) {
		DMESGF(DMCtrace, "Finding polynomials of degree %d...", n);
		if ( (Unsolved() == PHDFresponse) ?
				!FindPoly(&thisSoln, n, &patchList) :
				!FindSoln(&thisSoln, n, &patchList) ) {
			DMESGF(DMCwarning, "No solution of degree %d", n);
			continue;
		}
		thisErr = CompErr(&thisSoln, &patchList);
#ifdef phd_debug
{						// print solution to stderr
	int	j, k;
	fprintf(stderr, "\nExposure adjustment factors:\n");
	for (j = ustart+1; j < uend; j++)
		fprintf(stderr, "\t%f", thisSoln.expadj[j]);
	fprintf(stderr, "\nPolynomial fits:\n");
	for (k = 0; k < 3; k++) {
		fprintf(stderr, "%d: ", k);
		for (j = thisSoln.ply[k].degree; j > 1; j--)
			fprintf(stderr, "%f*x^%d + ", thisSoln.ply[k].coef[j], j);
		fprintf(stderr, "%f*x + %f\n", thisSoln.ply[k].coef[1],
				thisSoln.ply[k].coef[0]);
	}
}
#endif
		DMESGF(DMCtrace, "RMS error is %f", thisErr);
		if (thisErr < bestErr) {
			bestErr = thisErr;
			soln = thisSoln;
		}
	}
						// check that we found solution
	if (bestErr <= P_MAXERR) {
		haveSoln |= Unsolved();
		return true;
	}
	DMESG(DMCdata, "Cannot solve for response function");
	return false;
}

// Clear all or part of our solution
void
PHDImageMaker::ClearSoln(int clrWhat)
{
	int	i;
	if (clrWhat & PHDFresponse)
		for (i = 3; i--; )
			soln.ply[i] = linearPoly1;
	if (clrWhat & PHDFexposure)
		for (i = P_MAXHDR; i--; )
			soln.expadj[i] = 1.f;
	if (clrWhat & PHDFalignment)
		for (i = P_MAXHDR; i--; ) {
			ipi[i].xoff = ipi[i].yoff = 0;
			ipi[i].degCW = 0;
		}
	ord[0] = -1;
	ustart = uend = -1;
	haveSoln &= ~clrWhat;
}

// Save response function to camera response record
bool
PHDImageMaker::GetResponse(DBRecord *crp)
{
	if (crp == NULL)
		return false;
	if (!Compute())
		return false;
	const DBRecord &	srcr = ipi[0].pci->ircd;
	bool			ok = true;
	DBField			dbf;
	crp->Init(&CDBFInfo);
	ok &= crp->SetField(CDBFmake, srcr);
	ok &= crp->SetField(CDBFmodel, srcr);
	ok &= crp->SetField(CDBFversion, srcr);
	ok &= soln.ply[0].Get(&dbf) && crp->SetField(CDBFred, dbf);
	ok &= soln.ply[1].Get(&dbf) && crp->SetField(CDBFgreen, dbf);
	ok &= soln.ply[2].Get(&dbf) && crp->SetField(CDBFblue, dbf);
	return ok;
}

// Set response function from camera response record
bool
PHDImageMaker::SetResponse(const DBRecord &cr)
{
	if (cr.GetFieldInfo() != &CDBFInfo)
		return false;
	haveSoln &= ~PHDFresponse;
	if (!soln.ply[0].Set(cr[CDBFred]))
		return false;
	soln.ply[1].Set(cr[CDBFgreen]);
	soln.ply[2].Set(cr[CDBFblue]);
	if ((soln.ply[1].degree != soln.ply[0].degree) |
			(soln.ply[2].degree != soln.ply[1].degree))
		return false;
	haveSoln |= PHDFresponse;
	solveFlags &= ~PHDFresponse;
	return true;
}

// Evaluate response curve (and derivative) for indicated primary
bool
PHDImageMaker::GetResponse(float resp[256], int pri, float *der)
{
	const PPolynomial *	pp = GetPolynomial(pri);
	if (pp == NULL)
		return false;
	float	ftemp;
	float *	dp = &ftemp;
	int	itoe = 256;
	while (itoe > P_LINRESP) {
		--itoe;
		if (der)
			dp = der + itoe;
		double	tv = PbyteVal(itoe);
		if ((resp[itoe] = (float)pp->Eval(tv)) <= 0 ||
				(*dp = (float)pp->Deriv(tv)) <= 0) {
			++itoe;
			break;
		}
	}
	int		i;			// force linear toe
	const float	tder = resp[itoe]/PbyteVal(itoe);
	for (i = itoe; i--; )
		resp[i] = PbyteVal(i) * tder;
	if (der)
		for (i = itoe; i--; )
			der[i] = tder;
	return true;
}

// Identify best image exposure, optionally linking to it
int
PHDImageMaker::GetBestExposure(ImgStruct *ibp)
{
	int	curBest = nipi/2;		// find starting point
	int	dir = 1;
	while (curBest < nipi-1 && SkipExposure(curBest))
		++curBest;
	if (curBest >= nipi) {
		dir = -1;
		curBest = nipi/2;
		do
			if (--curBest < 0)
				return -1;	// should never happen
		while (SkipExposure(curBest));
	}
	if (!Inp(curBest).ComputeUsability())
		return -1;			// should never happen
	uint32	curTotal = Inp(curBest).usepix.SumTotal();
	bool	improved = false;		// now, see if we can do better
	uint32	triTotal;
	while (!SkipExposure(curBest+dir) &&
			Inp(curBest+dir).ComputeUsability() &&
			(triTotal = Inp(curBest+dir).usepix.SumTotal())
					>= curTotal) {
		curBest += dir;
		curTotal = triTotal;
		improved = true;
	}
	if (!improved)				// look the other way?
		while (!SkipExposure(curBest-dir) &&
			Inp(curBest-dir).ComputeUsability() &&
			(triTotal = Inp(curBest-dir).usepix.SumTotal())
					>= curTotal) {
			curBest -= dir;
			curTotal = triTotal;
			improved = true;
		}
	if (!ibp)				// no image to link?
		return curBest;
	if (haveSoln & PHDFalignment) {		// realign to best?
		for (int i = nipi; i--; ) {
			if (i == ord[curBest]) continue;
			ipi[i].xoff -= Inp(curBest).xoff;
			ipi[i].yoff -= Inp(curBest).yoff;
			ipi[i].degCW -= Inp(curBest).degCW;
		}
		Inp(curBest).xoff = Inp(curBest).yoff = 0;
		Inp(curBest).degCW = 0;
	}
	if (!Inp(curBest).GetReady()) {
		DMESG(DMCresource, "Cannot load image");
		return -1;
	}
						// link loaded image
	PlinkImage(ibp, &Inp(curBest).imr);
	Inp(curBest).Release();
	return curBest;
}

// Perform linear interpolation along a row or column of pixels
static void
pLinterp(float *pp, const int wid, const int step, const int nc)
{
	const float *	p0 = pp;
	const float *	p1 = pp + (wid-1)*step;
	const float	mul = 1.f/float(wid);
	float		sf;
	int		i, j;
	
	for (i = 1; i < wid-1; i++) {
		pp += step;
		sf = mul*float(i);
		for (j = nc; j--; )
			pp[j] = (1.f-sf)*p0[j] + sf*p1[j];
	}
}

// Remove bad pixels detected in variance image
static int
cleanBadPixels(ImgStruct *ims, const ImgStruct *ivs)
{
	const int		ncomp = ImgPixelLen[ims->csp->format];
	const int		psiz = ncomp * sizeof(float);
	const int		rowlen = ims->rowsize / sizeof(float);
	const int		SKrad = 2;
	const int		SKwid = 2*SKrad + 1;
	static const float	SKImgData[SKwid][SKwid] = {
					-8.f, -4.f, .0f, -4.f, -8.f,
					-4.f,  .0f, .0f,  .0f, -4.f,
					 .0f,  .0f, 1.f,  .0f,  .0f,
					-4.f,  .0f, .0f,  .0f, -4.f,
					-8.f, -4.f, .0f, -4.f, -8.f,
				};
	const float		SKthresh = 0.1f;
	int			nBad = 0;
	ImgStruct		spikeKern, spikeImg;
	int			x, y, i;
						// detect spikes in variance
	spikeKern.csp = &ICS_Y;
	spikeKern.xres = spikeKern.yres = SKwid;
	spikeKern.rowsize = SKwid*sizeof(float);
	spikeKern.mbase = NULL;
	spikeKern.img = const_cast<uby8 *>((const uby8 *)SKImgData);
	spikeImg.csp = &ICS_Y;
	spikeImg.img = NULL;
	if (!PconvolveImage(&spikeImg, ivs, &spikeKern))
		return -1;
	for (y = SKrad; y < spikeImg.yres-SKrad; y++) {
		const float *	sp = (const float *)ProwPtr(&spikeImg,y) + SKrad;
		for (x = SKrad; x < spikeImg.xres-SKrad; x++, sp++)
			if (*sp >= SKthresh) {	// interpolate from corners
				pLinterp((float *)PpixP(ims,x-SKrad,y-SKrad,psiz),
						SKwid, rowlen, ncomp);
				pLinterp((float *)PpixP(ims,x+SKrad,y-SKrad,psiz),
						SKwid, rowlen, ncomp);
				for (i = -SKrad; i <= SKrad; i++)
					pLinterp((float *)PpixP(ims,x-SKrad,y+i,psiz),
						SKwid, ncomp, ncomp);
				++nBad;
			}
	}
	PfreeImage(&spikeImg);
	return nBad;
}

// Maximum of two vaues
static inline float
fmaxi(float a, float b)
{
	return (a > b) ? a : b;
}

// Minimum of two values
static inline float
fmini(float a, float b)
{
	return (a > b) ? b : a;
}

// Compute weight reduction based on reference image pixel
static float
weightAdjust(const uby8 rp[3], int rng, const float rc[3], const float pc[3])
{
	float	d;

	if (rng > 0) {				// overrange condition
		if ((pc[0] >= .9f*rc[0]) & (pc[1] >= .9f*rc[1]) &
				(pc[2] >= .9f*rc[2]))
			return 0.0256f/(0.0256f + 0.1f*0.1f);
		d = 1.f - fmini(pc[0]/rc[0], fmini(pc[1]/rc[1], pc[2]/rc[2]));
		return 0.0256f / (0.0256f + d*d);
	}
	if (rng < 0) {				// underrange condition
		if ((pc[0] <= 1.25f*rc[0]) & (pc[1] <= 1.25f*rc[1]) &
				(pc[2] <= 1.25f*rc[2]))
			return 0.0256f/(0.0256f + 0.25f*0.25f);
		if (fmini(rc[0], fmini(rc[1], rc[2])) <= 0)
			return 1.f;		// punt
		d = 1.f - fmaxi(pc[0]/rc[0], fmaxi(pc[1]/rc[1], pc[2]/rc[2]));
		return 0.0256f / (0.0256f + d*d);
	}
	float	e2 = PbyteVal(rgb_bright(rp));
	if (e2 > 0.85f)
		e2 = 0.058f + (e2 - 0.85f)*0.68f;
	else
		e2 = 0.04f + (1.f - e2)*0.12f;
	e2 *= e2;
	float	d2 = 0;				// find maximum difference
	for (int j = 3; j--; ) {
		d = 1.f - pc[j]/rc[j];
		d *= d;
		if (d2 < d) d2 = d;
	}
	return e2 / (e2 + d2);
}

// Compute and return high dynamic-range image
bool
PHDImageMaker::GetHDImage(ImgStruct *ims, ImgStruct *ivs)
{
	static const char	cancelMsg[] = "User cancelled HDR builder";
	const float		minderi = .2f;
	ImgColorSpace		tempCS;
	float			resp[3][256];
	float			deri[3][256];
	float			mult[P_MAXHDR];
	int			refExp = -1;
	int			j, n, x, y;
						// check color space
	if (ims == NULL)
		return false;
	if (ims->csp == NULL) {
		if (ims->img != NULL) {
			DMESG(DMCparameter, "Destination missing color space");
			return false;
		}
	} else if (ImgPixelLen[ims->csp->format] != 3 ||
			!PmatchColorSpace(ims->csp, &ICS_RGB709,
						PICMdtype|PICMgamma)) {
		DMESG(DMCparameter, "Only support HDR in linear float");
		return false;
	}
						// compute response curves
	if (!GetResponse(resp[0],0,deri[0]) ||
			!GetResponse(resp[1],1,deri[1]) ||
			!GetResponse(resp[2],2,deri[2]))
		return false;
						// limit derivatives
	for (n = 256; n--; )
		for (j = 3; j--; )
			if (deri[j][n] < minderi) deri[j][n] = minderi;
						// get calibration factors
	for (n = nipi; n--; )
		mult[n] = GetStoNits(n)/GetStoNits();
						// start progress
	if (reportProgress && !(*reportProgress)("Combining Images", 0)) {
		DMESG(DMCtrace, cancelMsg);
		return false;
	}
						// allocate destination
	for (n = nipi; n--; ) {
		if (!SkipExposure(n))
			break;
		if (reportProgress && !(*reportProgress)(NULL, 15*(nipi-n)/nipi)) {
			DMESG(DMCtrace, cancelMsg);
			return false;
		}
	}
						// usability always 1 ahead
	if (!Inp(n).ComputeUsability())
		return false;
	if (reportProgress && !(*reportProgress)(NULL, 15)) {
		DMESG(DMCtrace, cancelMsg);
		return false;
	}
	if (ims->csp == NULL) {			// set up destination CS
		PcopyCS(&tempCS, Inp(n).imr.csp);
		PrealCS(&tempCS);
		ims->csp = &tempCS;
	}
	if ((ims->xres <= 0) | (ims->yres <= 0)) {
		ims->xres = Inp(n).imr.xres;
		ims->yres = Inp(n).imr.yres;
	}
	if (!PnewImage(ims, (double)Inp(n).imr.yres/Inp(n).imr.xres))
		return false;
	if ((ims->xres != Inp(n).imr.xres) | (ims->yres != Inp(n).imr.yres)) {
		DMESG(DMCparameter, "Cannot resample HDR image");
		PfreeImage(ims);
		return false;
	}
						// allocate weight sum image
	ImgStruct	wsum, ivar, iref;
	wsum.xres = ims->xres; wsum.yres = ims->yres;
	wsum.csp = &ICS_RGB709;
	wsum.img = NULL;
	if (!PnewImage(&wsum, .0)) {
		PfreeImage(ims);
		return false;
	}
	memset(&ivar, 0, sizeof(ivar));
	memset(&iref, 0, sizeof(iref));
#define freeImages()		PfreeImage(&wsum); PfreeImage(&ivar); \
				PfreeImage(&iref)
						// allocate variance image
	if (nipi > 1 && (ivs != NULL) | !(solveFlags & PHDFalignment)) {
		ivar = *ims;
		ivar.img = NULL;
		if (!PsetImage(&ivar, Pblack)) {
			freeImages();
			PfreeImage(ims);
			return false;
		}
	}
						// zero image and weights
	PclearImage(ims, NULL);
	PclearImage(&wsum, NULL);
						// get best exposure if needed
	if (solveFlags & PHDFghosting && (refExp = GetBestExposure(&iref)) < 0)
		DMESG(DMCwarning, "Cannot identify reference exposure");
						// accumulate each exposure
	for (n = nipi; n-- > 0; ) {
		if (SkipExposure(n))
			continue;
		if (n > 0 && !SkipExposure(n-1))
			Inp(n-1).ComputeUsability();
		PHDExposure &	ipe = Inp(n);
		DMESGF(DMCtrace, "Adding in exposure %d...", n);
		if (!ipe.GetReady()) {
			DMESG(DMCresource, "Cannot load image");
			freeImages();
			PfreeImage(ims);
			return false;
		}
		if ((ims->xres != ipe.imr.xres) | (ims->yres != ipe.imr.yres)) {
			DMESG(DMCdata, "Exposure resolution mismatch");
			freeImages();
			PfreeImage(ims);
			return false;
		}
		float	swt[3][256];		// compute weighting functions
		for (j = 3; j--; ) {
			float		lastceil = 0;
			int		minv;
			if (n > ustart) {	// avoid noise but ensure overlap
				minv = Inp(n-1).maxv[j] - P_MAXNOISE/2;
				minv *= (minv > 0);
				lastceil = mult[n-1]*resp[j][minv];
			}
			for (minv = P_MAXNOISE; minv > ipe.minv[j]; minv--)
				if (mult[n]*resp[j][minv+P_MAXNOISE/2] < lastceil)
					break;
			if ((n == refExp) | (minv < ipe.minv[j]))
				minv = ipe.minv[j];
			const int	mid = (ipe.maxv[j] + minv)/2;
			for (int i = 256; i--; ) {
						// longer exposures contribute less noise
				swt[j][i] = resp[j][i]*mult[uend-1]/(deri[j][i]*mult[n]);
				if (i <= P_MAXNOISE && swt[j][i] > 0.95f*swt[j][i+1])
					swt[j][i] = 0.95f*swt[j][i+1];
				if (swt[j][i] < 1.05f*P_WTMIN)
					swt[j][i] = 1.05f*P_WTMIN;
				if ((n <= ustart) & (i <= mid))
					continue;
				if ((n >= uend-1) & (i >= mid))
					continue;
				float	f;	// mult. by hat function
				f = float(i-minv) / float(ipe.maxv[j] - minv);
				f = 2.f*f - 1.f;
				f *= f; f *= f; f *= f*f;
				swt[j][i] *= 1.f - f;
				if (swt[j][i] < 0)
					swt[j][i] = 0;
			}
		}
		ImgStruct	thisImg;	// rotate image if needed
		if (fabs(ipe.degCW) > P_ROTHRESH) {
			thisImg = ipe.imr;
			thisImg.img = NULL;
			if (!ProtateImage(&thisImg, &ipe.imr, -ipe.degCW,
						PSfast, Pblack)) {
				ipe.Release();
				freeImages();
				PfreeImage(ims);
				return false;
			}
		} else				// otherwise link & release
			PlinkImage(&thisImg, &ipe.imr);
		ipe.Release();
		for (y = ims->yres; y--; ) {	// compute weighted samples
			if (y + ipe.yoff < 0)
				break;
			if (y + ipe.yoff >= thisImg.yres)
				continue;
			const uby8 *	sp = PpixP(&thisImg, ipe.xoff,
							y+ipe.yoff, 3);
			float *		dp = (float *)ProwPtr(ims, y);
			float *		wp = (float *)ProwPtr(&wsum, y);
			float *		ip = NULL;
			if (ivar.img)
				ip = (float *)ProwPtr(&ivar, y);
			for (x = 0; x < ims->xres; x++, sp += 3, dp += 3, wp += 3) {
				if (x + ipe.xoff >= thisImg.xres)
					break;
				if (x + ipe.xoff < 0) {
					if (ip) ip += 3;
					continue;
				}
				float	wt[3], pc[3];
				float	wtlim = 0;
				for (j = 3; j--; ) {
					pc[j] = mult[n] * resp[j][ sp[j] ];
					if ((wt[j] = swt[j][ sp[j] ]) > wtlim)
						wtlim = wt[j];
				}
						// check against reference exp.
				if ((refExp >= 0) & (n != refExp)) {
					const uby8 *	rp = PpixP(&iref,x,y,3);
					int		rng = 0;
					if (!Inp(refExp).Usable(x,y))
						rng = Inp(refExp).OverRange(rp) ?
								1 : -1;
					float		rc[3];
					for (j = 3; j--; )
						rc[j] = mult[refExp] *
							    resp[j][ rp[j] ];
					const float	de = weightAdjust(rp,
								rng, rc, pc);
					wt[0] *= de; wt[1] *= de; wt[2] *= de;
					wtlim *= de;
				}
				wtlim *= 0.05f;	// ameliorates primary bias
				for (j = 3; j--; ) {
					if (wt[j] < wtlim) wt[j] = wtlim;
					wp[j] += wt[j];
					dp[j] += wt[j] * pc[j];
				}
				if (ip)		// accumulate sum of squares
					for (j = 0; j < 3; j++)
						*ip++ += wt[j]*pc[j]*pc[j];
			}
		}
		PfreeImage(&thisImg);		// done with this exposure
		if (reportProgress && !(*reportProgress)(NULL,
					15 + 70*(uend-n)/(uend-ustart))) {
			DMESG(DMCtrace, cancelMsg);
			freeImages();
			PfreeImage(ims);
			return false;
		}
	}
						// divide pixels by weights
	const float	mino = 0.5f * mult[ustart] *
				resp[1][ Inp(ustart).minv[1] ];
	const int	nextrow = ims->rowsize/sizeof(float);
	for (y = ims->yres; y--; ) {		// bottom to top
		float *		dp = (float *)ProwPtr(ims, y);
		const float *	wp = (float *)ProwPtr(&wsum, y);
		float *		ip = NULL;
		if (ivar.img)
			ip = (float *)ProwPtr(&ivar, y);
		for (x = 0; x < 3*ims->xres; x++, wp++, dp++) {
			if (*wp > P_WTMIN) {	// good primary
				*dp /= *wp;
			} else {		// else fill from neighbors
				if ((x >= 3) & (y < ims->yres-1))
					*dp = .495f*(dp[-3] + dp[nextrow]);
			}
						// check floor
			if (*dp < mino) {
				*dp = mino;
				if (ip) *ip++ = 0;
				continue;
			}
			if (!ip) continue;
						// compute variance
			if (*wp > P_WTMIN) {
				*ip /= *wp * *dp * *dp;
				*ip++ -= 1.f;
			} else
				*ip++ = 0;
		}
	}
	PfreeImage(&wsum);
						// filter bad pixels?
	if (!(solveFlags & PHDFalignment) & (nipi > 1)) {
		if (reportProgress && !(*reportProgress)(
					"Cleaning Bad Pixels", 93)) {
			DMESG(DMCtrace, cancelMsg);
			freeImages();
			PfreeImage(ims);
			return false;
		}
		cleanBadPixels(ims, &ivar);
	}
	if (reportProgress && !(*reportProgress)(NULL, 100)) {
		DMESG(DMCtrace, cancelMsg);
		freeImages();
		PfreeImage(ims);
		return false;
	}
						// save variance image?
	if ((ivs != NULL) & (ivar.img != NULL)) {
		if ((ivs->csp == NULL) | (ivs->xres <= 0) | (ivs->yres <= 0) ||
				(ivs->img == NULL) &
				(ivs->xres == ivar.xres) &
				(ivs->yres == ivar.yres) &
				PmatchColorSpace(ivs->csp, ivar.csp, PICMall))
			PlinkImage(ivs, &ivar);
		else
			PrenderImageI(ivs, &ivar);
	}
	freeImages();
						// need to convert CS?
	if (!PmatchColorSpace(ims->csp, Inp(ustart).imr.csp,
						PICMformat|PICMprims)) {
		const ImgColorSpace *	cstarget = ims->csp;
		ImgColorSpace		cscurrent;
		PcopyCS(&cscurrent, cstarget);
		cscurrent.format = Inp(ustart).imr.csp->format;
		memcpy(cscurrent.chroma, Inp(ustart).imr.csp->chroma,
					sizeof(cscurrent.chroma));
		ims->csp = &cscurrent;
		if (!PconvertColorSpace(ims, cstarget, 1.f)) {
			DMESG(DMCwarning, "Cannot convert HDR color space");
			ims->csp = cstarget;	// lie
		} else
			DMESG(DMCtrace, "Converted HDR color space");
		// XXX Should correct variance image (ivs) color space?
	}
	DMESG(DMCtrace, "Done building HDR image");
	return true;
#undef freeImages
}
