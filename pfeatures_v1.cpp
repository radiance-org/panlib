/*
 *  pfeatures.cpp
 *  pan
 *
 *  Modules for computing and matching image features using U-BRIEF descriptors:
 *	"BRIEF: Computing a local binary descriptor very fast"
 *	Calonder et al., Proc. ECCV 2010.
 *
 *  Created by Greg Ward on 2/24/15.
 *  Copyright 2015 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dmessage.h"
#include "panimage.h"
#include "pfeatures.h"
#include "gaussjord.h"
#include "interp2d.h"

#if defined(_WIN32) || defined(_WIN64)
#define drand48()	(((double)rand())*(1./RAND_MAX))
#endif

// BRIEF sampling pair
static struct BriefSamp {
	short		xo1, yo1;
	short		xo2, yo2;
}	ubrief[BRIEF_BITS];		// U-BRIEF sampling pattern

// convert 0-1 value into normal distribution with mean=0 & st.dev=1
static inline double
unif2norm(double x)
{
	double	sgn = 1;

	if (x < .5)
		sgn = -1;
	else
		x = 1. - x;

	if (x <= 0)
		x = 1e-6;

	x = sqrt(-2.*log(x));

	return sgn*(x - (2.515517+x*(.802853+x*.010328))) /
			(1. + x*(1.432788+x*(.189269+x*.001308)));
}

static inline int
iround(double x)
{
	return int(x + .5) - int(x < -.5);
}

// check if patch offset is out of range
static inline bool
badRange(int o)
{
	return (o < -BRIEF_DIA/2) | (o >= BRIEF_DIA/2);
}

// initialize U-BRIEF sampling pattern II (internal)
static void
initBRIEF()
{
	const double	stdev = BRIEF_DIA/5.;
	int		i = BRIEF_BITS;

	while (i--)
	    do {

		ubrief[i].xo1 = iround(stdev*unif2norm(drand48()));
		ubrief[i].yo1 = iround(stdev*unif2norm(drand48()));
		ubrief[i].xo2 = iround(stdev*unif2norm(drand48()));
		ubrief[i].yo2 = iround(stdev*unif2norm(drand48()));

	    } while (badRange(ubrief[i].xo1) | badRange(ubrief[i].yo1) |
			badRange(ubrief[i].xo2) | badRange(ubrief[i].yo2) ||
			(ubrief[i].xo1 == ubrief[i].xo2) &
				(ubrief[i].yo1 == ubrief[i].yo2));
}

// compute image derivative in the given direction (internal)
static bool
derivatImage(PanImage *ib, const PanImage &im, int dx, int dy)
{
	int	x, y;

	if (!ib->Init(im.Width(), im.Height(), ICS_Y, PICzero))
		return false;
	for (y = dy; y < im.Height()-dy; y++) {
		const float	*ip0 = im.GetRowFloat(y-dy) + dx;
		const float	*ip1 = im.GetRowFloat(y+dy) + dx;
		float		*dp = ib->RowFloat(y + dy/2) + dx/2;
		for (x = dx; x < im.Width()-dx; x++, ip0++, ip1++)
			*dp++ = ip1[dx] - ip0[-dx];
	}
	*ib *= 0.5f/(float)sqrt(dx*dx + dy*dy);
	return true;
}

// identify corners using Harris detector (internal)
static bool
cornerDetect(ABitMap2 *fmap, const PanImage &imxx, const PanImage &imyy,
		const PanImage &imxy, double drad)
{
	const float	thresh = 4e-4f;
	PanImage	crnIm(imxy.Width(), imxy.Height(), ICS_Y);
	int		x, y;

	if (!fmap->NewBitMap(imxy.Width(), imxy.Height()) || !crnIm.Init(PICzero))
		return false;
						// get candidates and values
	for (y = imxy.Height(); y--; ) {
		const float *	pxx = (const float *)imxx[y];
		const float *	pyy = (const float *)imyy[y];
		const float *	pxy = (const float *)imxy[y];
		float *		cscan = (float *)crnIm[y];
		for (x = 0; x < imxy.Width(); x++, pxx++, pyy++, pxy++) {
			float	v = (*pxx * *pyy - *pxy * *pxy)/(*pxx + *pyy + 1e-6f);
			if (v < thresh) continue;
			cscan[x] = v;
			fmap->Set(x, y);
		}
	}
	PanImage	dilIm(ICS_Y);		// dilate to identify local maxima
	if (!PdilateImage(dilIm.Img(), crnIm.GetImg(), drad))
		return false;
	for (x = y = 0; fmap->Find(&x, &y); x++)
		if (((const float *)dilIm[y])[x] > ((const float *)crnIm[y])[x])
			fmap->Reset(x, y);
						// eliminate clusters
	ABitMap2	clusterMap(*fmap);
	clusterMap.Expand(-1.2);
	clusterMap.Expand(1.6);
	return fmap->ClearBitsFrom(clusterMap);
}

// Normalize float Y image dividing by local blur
static bool
normGray(PanImage *imp, float blur)
{
	PanImage	bim(imp->Width(), imp->Height(), ICS_Y);

	if (!PblurImage(bim.Img(), imp->GetImg(), blur))
		return false;
					// leave black pixels alone
	for (int y = bim.Height(); y--; ) {
		const float *	ip = (const float *)(*imp)[y];
		float *		p = (float *)bim[y];
		for (int n = bim.Width(); n--; p++)
			if (*ip++ <= 5e-6f)
				*p = 1.f;
	}
	*imp /= bim;
	return true;
}

// Normalize image to emphasize edges, returning float gray map
bool
PnormEdges(PanImage *imp, const PanImage &im, float blur)
{
	static const float	gcoef[3] = {0.21264f, 0.71517f, 0.07219f};

	if (!imp || im.Ready() < PICread || blur <= 1.f)
		return false;
	if (imp->Ready() != PICready) {
		imp->Init(ICS_Y);
	} else if ((imp->Width() != im.Width()) | (imp->Height() != im.Height())) {
		DMESG(DMCparameter, "Resolution mismatch in PnormEdges");
		return false;
	} else if (!PmatchColorSpace(imp->GetCS(), &ICS_Y, PICMall)) {
		DMESG(DMCparameter, "Output color space must be float gray");
		return false;
	}
	if (im.NComp() == 1)		// processing single channel input?
		return imp->Load(im) && normGray(imp, blur);

	if (im.NComp() != 3) {		// else assume tristimulus
		DMESGF(DMCparameter, "Input cannot have %d components", im.NComp());
		return false;
	}
					// find best channel discriminators
	PanImage	dilIm(ICS_RGB709), eroIm(ICS_RGB709);
	if (!PdilateImage(dilIm.Img(), im.GetImg(), blur) ||
			!PerodeImage(eroIm.Img(), im.GetImg(), blur))
		return false;
	for (int y = im.Height(); y--; ) {
		const float *	ergb = (const float *)eroIm[y];
		float *		drgb = (float *)dilIm[y];
		for (int n = im.Width(); n--; ergb += 3, drgb += 3) {
			float	gf = gcoef[0]*ergb[0] + gcoef[1]*ergb[1]
						+ gcoef[2]*ergb[2];
			if (gf <= 5e-6f) {
				drgb[0] = drgb[1] = drgb[2] = 1.f;
				continue;
			}
			gf *= 0.1f;
			int	bestc = 0;
			float	bestv = drgb[0]/(ergb[0] + gf);
			int	i = 3;
			while (--i) {
				float	candv = drgb[i]/(ergb[i] + gf);
				if (candv <= bestv) continue;
				bestc = i;
				bestv = candv;
			}		// only one channel is non-zero
			for (i = 3; i--; )
				drgb[i] = (i==bestc)/(gcoef[i]*drgb[i] + 1e-7f);
		}
	}
	eroIm.Init();
	dilIm ^= blur;			// orig * blurred channel coefficients
	if (!eroIm.Load(im))
		return false;
	eroIm *= dilIm;
	dilIm.Init();
	return imp->Load(eroIm);	// combines to optimized edge channel
}

// Compute feature (corner) bitmap
bool
PdetectFeatures(ABitMap2 *fmap, const PanImage &im, float blur)
{
	const float	blur_kernel = blur + 1.f;
	PanImage	imx, imy, imxy;
	PanImage	cornIm;
					// normalize edge information
	if (!PnormEdges(&imxy, im, 3.f*blur_kernel))
		return false;
					// compute derivative images
	if (!derivatImage(&imx, imxy, iround(blur_kernel), 0))
		return false;
	if (!derivatImage(&imy, imxy, 0, iround(blur_kernel)))
		return false;
	imxy = imx;			// square derivatives
	imxy *= imy;
	imx *= imx;
	imy *= imy;
	imxy ^= blur_kernel;		// blur them! (can use box filter)
	imx ^= blur_kernel;
	imy ^= blur_kernel;
					// detect corners
	return cornerDetect(fmap, imx, imy, imxy, 7.f*blur_kernel);
}

// Identify features in the given image, return count
int
PanFeatureList::FindFeatures(const PanImage &im, float blur)
{
	PanImage	imY;
	ABitMap2	featMap;

	Clear();
					// detect features (in color if provided)
	if (!PdetectFeatures(&featMap, im, blur))
		return -1;
					// get grayscale version of image
	if (blur > 0.5f || !PmatchColorSpace(im.GetCS(), &ICS_Y, PICMall)
			|| !imY.LinkRO(im)) {
		imY.Init(ICS_Y);
		if (!imY.Load(im))
			return -1;
	}
	if (blur > 0.5f)		// apply smoothing if requested
		imY ^= blur;
					// avoid positions in margins
	xres = imY.Width();
	yres = imY.Height();
	featMap.ClearRect(0, 0, xres, BRIEF_DIA/2);
	featMap.ClearRect(0, yres-BRIEF_DIA/2, xres, BRIEF_DIA/2);
	featMap.ClearRect(0, BRIEF_DIA/2, BRIEF_DIA/2, yres-BRIEF_DIA);
	featMap.ClearRect(xres-BRIEF_DIA/2, BRIEF_DIA/2, BRIEF_DIA/2, yres-BRIEF_DIA);
					// count up what's left
	if (!(nfeats = featMap.SumTotal()))
		return 0;
					// allocate & compute features
	feat = new ImgFeature [nfeats];
	if (!ubrief[0].xo1 && !ubrief[0].yo1 && !ubrief[0].xo2 && !ubrief[0].yo2)
		initBRIEF();
	int	n, x, y;		// evaluate U-BRIEF descriptors
	for (n = x = y = 0; featMap.Find(&x, &y); x++, n++) {
		feat[n].x = x; feat[n].y = y;
		feat[n].descr.NewBitMap(BRIEF_BITS);
		for (int i = BRIEF_BITS; i--; ) {
			const BriefSamp	&bs = ubrief[i];
			if (((const float *)imY[y+bs.yo1])[x+bs.xo1] <
					((const float *)imY[y+bs.yo2])[x+bs.xo2])
				feat[n].descr.Set(i);
		}
	}
	return nfeats;
}

// Remove possibly confounding features from our list
int
PanFeatureList::CullSimilar(const double rad, const int thresh)
{
	ABitMap	cullMap(nfeats);
	uint32	i;
	int	j;			// identify similar descriptors
	for (i = 0; i < nfeats-1; i++)
	    for (j = i; ++j < nfeats; ) {
	    	if (rad >= 1 && feat[i].Dist2(feat[j]) > rad*rad)
			continue;	// outside specified radius
		if (feat[i].Fdist(feat[j]) <= thresh) {
			cullMap.Set(i);
			cullMap.Set(j);
		}
	    }
	const int	nc = cullMap.SumTotal();
	if (!nc)
		return 0;
	if (nc == nfeats) {
		Clear();
		return nc;
	}
	ImgFeature *	ofeat = feat;	// copy to new list
	nfeats -= nc;
	feat = new ImgFeature [nfeats];
	i = j = 0;
	while (cullMap.Find(&i, false))
		feat[j++].Take(&ofeat[i++]);
	delete [] ofeat;		// free old list
	return nc;
}

// structure for sorting out features
struct FeatEntry {
	int	fi;	// feature index
	int	di;	// feature distance
};

// comparison function puts closest feature at head
static int
feat_cmp(const void *p1, const void *p2)
{
	return (*(const FeatEntry *)p1).di - (*(const FeatEntry *)p2).di;
}

// Find closest n features
int
PanFeatureList::BestFeatures(int fndx[], int n, const ImgFeature &mf,
				const int thresh) const
{
	if (!fndx)
		return 0;
	memset(fndx, 0xff, n*sizeof(fndx[0]));
	if (n > nfeats)
		n = nfeats;
	if (n <= 0)
		return 0;
	int	i;
	if (n == 1) {			// special case
		int	best_di = thresh+1;
		for (i = nfeats; i--; ) {
			int	di = mf.Fdist(feat[i]);
			if (di >= best_di)
				continue;
			fndx[0] = i;
			best_di = di;
		}
		return (fndx[0] >= 0);
	}
					// else use qsort()
	FeatEntry *	sarr = new FeatEntry [nfeats];
	int		ngood = 0;
	for (i = nfeats; i--; )
		if ((sarr[ngood].di = mf.Fdist(feat[i])) <= thresh)
			sarr[ngood++].fi = i;
	qsort(sarr, ngood, sizeof(FeatEntry), feat_cmp);
	if (n > ngood)
		n = ngood;
	for (i = n; i--; )
		fndx[i] = sarr[i].fi;
	delete [] sarr;
	return n;
}

// Find nearest n features (within threshold)
int
PanFeatureList::NearFeatures(int fndx[], int n, const ImgFeature &mf,
				const int thresh) const
{
	if (!fndx)
		return 0;
	memset(fndx, 0xff, n*sizeof(fndx[0]));
	if (n > nfeats)
		n = nfeats;
	if (n <= 0)
		return 0;
	int	i;
	if (n == 1) {			// special case
		int	best_di = xres*xres + yres*yres;
		for (i = nfeats; i--; ) {
			if (thresh < BRIEF_BITS &&
					mf.Fdist(feat[i]) > thresh)
				continue;
			int	di = mf.Dist2(feat[i]);
			if (di >= best_di)
				continue;
			fndx[0] = i;
			best_di = di;
		}
		return (fndx[0] >= 0);
	}
					// else use qsort()
	FeatEntry *	sarr = new FeatEntry [nfeats];
	int		ngood = 0;
	for (i = nfeats; i--; ) {
		if (thresh < BRIEF_BITS && mf.Fdist(feat[i]) > thresh)
			continue;
		sarr[ngood].di = mf.Dist2(feat[i]);
		sarr[ngood++].fi = i;
	}
	qsort(sarr, ngood, sizeof(FeatEntry), feat_cmp);
	if (n > ngood)
		n = ngood;
	for (i = n; i--; )
		fndx[i] = sarr[i].fi;
	delete [] sarr;
	return n;
}

// Copy operator
PanFeatureList &
PanFeatureList::operator=(const PanFeatureList &orig)
{
	if (this == &orig)
		return *this;
	if (orig.nfeats > 0) {
		if (orig.nfeats != nfeats) {
			delete [] feat;
			feat = new ImgFeature [orig.nfeats];
		}
		for (nfeats = 0; nfeats < orig.nfeats; nfeats++ )
			feat[nfeats] = orig.feat[nfeats];
	} else
		Clear();
	xres = orig.xres;
	yres = orig.yres;
	return *this;
}

#define	NFEAT2CHECK		8	// maximum # features to match

// compute perspective transform from matched features (internal)
static bool
computePerspXform(float xfm[8], const PanFeatureList &fL,
		const ImgFeature *fmatch[][NFEAT2CHECK], const uby8 useMatch[])
{
	int	nOK = 0;
	int	i, j;
					// copy valid feature matches
	for (i = fL.Count(); i--; )
		nOK += (useMatch[i] < NFEAT2CHECK);
	if (nOK < 4)
		return false;
	float *	amat = new float [8*2*nOK];
	float *	bvec = new float [2*nOK];
	for (i = j = 0; j < fL.Count(); j++) {
		if (useMatch[j] >= NFEAT2CHECK)
			continue;
		const ImgFeature *	ffrom = fL.GetFeature(j);
		const ImgFeature *	fto = fmatch[j][useMatch[j]];
		float *	rowp = amat + 8*2*i;
		rowp[0] = ffrom->x;
		rowp[1] = ffrom->y;
		rowp[2] = 1;
		rowp[3] = rowp[4] = rowp[5] = 0;
		rowp[6] = -fto->x*ffrom->x;
		rowp[7] = -fto->x*ffrom->y;
		rowp += 8;
		rowp[0] = rowp[1] = rowp[2] = 0;
		rowp[3] = ffrom->x;
		rowp[4] = ffrom->y;
		rowp[5] = 1;
		rowp[6] = -fto->y*ffrom->x;
		rowp[7] = -fto->y*ffrom->y;
		bvec[2*i] = fto->x;
		bvec[2*i+1] = fto->y;
		++i;
	}
	bool	ok = GJsolveLeastSq(8, 2*nOK, amat, bvec, xfm);
	delete [] amat;
	delete [] bvec;
	return ok;
}

// map an (x,y) image position using a projective transform (internal)
static inline void
perspMap(float xyout[2], const float xfm[8], double xin, double yin)
{
	double	sca = xfm[6]*xin + xfm[7]*yin + 1.;

	if (sca == 0) {
		xyout[0] = xyout[1] = 0;
		return;
	}
	sca = 1./sca;
	xyout[0] = float(sca*(xfm[0]*xin + xfm[1]*yin + xfm[2]));
	xyout[1] = float(sca*(xfm[3]*xin + xfm[4]*yin + xfm[5]));
}

static inline double
ddist2(double x0, double y0, double x1, double y1)
{
	return (x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1);
}

// Compute warp grid (source positions) to register feature sets
bool
PregisterWarp(float warpGrid[][2], int whres, int wvres,
		const PanFeatureList &dstFL, const PanFeatureList &srcFL,
		const int thresh, PfeatureUsedMethod *cb, void *udp)
{
	if ((warpGrid == NULL) | (whres < 2) | (wvres < 2)) {
		DMESG(DMCparameter, "Warping grid must be at least 2x2");
		return false;
	}
	if ((dstFL.Count() < 4) | (srcFL.Count() < 4)) {
		DMESG(DMCdata, "Need at least 4 features from each image");
		return false;
	}
					// match up features
	const double		xsca = srcFL.Width()/double(dstFL.Width());
	const double		ysca = srcFL.Height()/double(dstFL.Height());
//	const ImgFeature *	fmatch[dstFL.Count()][NFEAT2CHECK];
//	uby8			useMatch[dstFL.Count()];
	const ImgFeature 	*(*fmatch)[NFEAT2CHECK] = new const ImgFeature * [dstFL.Count()][NFEAT2CHECK];
	uby8 *			useMatch = new uby8 [dstFL.Count()];
	double			avglength2 = 0;
	int			nOK = 0;
	const ImgFeature *	imfp;
	double			d2;
	int			i, j;
					// check scaling
	DTESTF((xsca*ysca > 1.44) | (xsca*ysca < 0.694), DMCwarning,
			"%.0f%% scale difference may affect feature matching",
			100.*(sqrt(xsca*ysca) - 1.));
					// start by taking best matches
	memset(useMatch, 0, sizeof(useMatch[0])*dstFL.Count());
	for (i = dstFL.Count(); i--; ) {
		imfp = dstFL.GetFeature(i);
		if (srcFL.BestFeatures(fmatch[i], NFEAT2CHECK, *imfp, thresh) > 0) {
			avglength2 += ddist2(xsca*imfp->x, ysca*imfp->y,
					fmatch[i][0]->x, fmatch[i][0]->y);
			++nOK;
		} else
			useMatch[i] = NFEAT2CHECK;
	}
					// initial preen based on average length
	for (i = dstFL.Count(); i--; ) {
		if (useMatch[i] >= NFEAT2CHECK)
			continue;
		imfp = dstFL.GetFeature(i);
		d2 = ddist2(xsca*imfp->x, ysca*imfp->y,
					fmatch[i][0]->x, fmatch[i][0]->y);
		while (d2*nOK > 4.*avglength2) {
			avglength2 -= d2;
			if (++useMatch[i] >= NFEAT2CHECK ||
					!fmatch[i][useMatch[i]]) {
				useMatch[i] = NFEAT2CHECK;
				--nOK;
				break;
			}
			avglength2 += d2 = ddist2(xsca*imfp->x, ysca*imfp->y,
				fmatch[i][useMatch[i]]->x,
				fmatch[i][useMatch[i]]->y);
		}
	}
	if (nOK < 4) {
		DMESG(DMCdata, "Need at least 4 matching features");
		delete [] fmatch;
		delete [] useMatch;
		return false;
	}
	avglength2 *= 1./nOK;
	float	perspXform[8];		// fit initial perspective transform
	if (!computePerspXform(perspXform, dstFL, fmatch, useMatch)) {
		DMESG(DMCdata, "Could not solve overall perspective transform");
		delete [] fmatch;
		delete [] useMatch;
		return false;
	}
	float	xyMapped[2];		// prefer matches that fit transform
	for (i = dstFL.Count(); i--; ) {
		if (useMatch[i] >= NFEAT2CHECK)
			continue;
		imfp = dstFL.GetFeature(i);
		perspMap(xyMapped, perspXform, imfp->x, imfp->y);
		double	bestd2 = 4.*avglength2;
		for (j = 0; j < NFEAT2CHECK && fmatch[i][j]; j++)
			if ((d2 = ddist2(xyMapped[0], xyMapped[1],
				  fmatch[i][j]->x, fmatch[i][j]->y)) < bestd2) {
				useMatch[i] = j;
				bestd2 = d2;
			}
					// check if feature is hopeless outlier
		if (bestd2 > 2.*avglength2) {
			useMatch[i] = NFEAT2CHECK;
			--nOK;
		}
	}
	if (nOK < 4) {
		DMESG(DMCdata, "Too many outliers");
		delete [] fmatch;
		delete [] useMatch;
		return false;
	}
	sprintf(dmessage_buf,"Computing %dx%d warp grid based on %d feature matches",
				whres, wvres, nOK);
	DMESG(DMCtrace, dmessage_buf);
	if (cb != NULL)			// report which features we're using
	    for (i = 0; i < dstFL.Count(); i++)
		for (j = 0; j < NFEAT2CHECK && fmatch[i][j]; j++) {
		    (*cb)(dstFL[i], *fmatch[i][j], (j==useMatch[i]), udp);
		    if (j == useMatch[i])
			break;		// skip worse matches we didn't use
		}
					// recompute mapping w/o outliers
	if (!computePerspXform(perspXform, dstFL, fmatch, useMatch)) {
		DMESG(DMCdata, "Could not resolve perspective transform");
		delete [] fmatch;
		delete [] useMatch;
		return false;
	}
	for (j = wvres; j--; )		// compute canonical warp grid values
		for (i = whres; i--; )
			perspMap(warpGrid[j*whres+i], perspXform,
						dstFL.Width()*i/(whres-1.),
						dstFL.Height()*j/(wvres-1.));
	if (nOK == 4 || (whres == 2) & (wvres == 2)) {
		delete [] fmatch;
		delete [] useMatch;
		return true;		// no interior to adjust
	}
					// else realign feature positions
	INTERP2 *	interp = interp2_alloc(nOK);
	if (interp == NULL) {
		DMESG(DMCmemory, "Not enough memory for 2-D interpolant");
		return false;
	}				// skipping outliers for interpolant
	float	(*xyAdj)[2] = new float [nOK][2];
	for (i = j = 0; j < dstFL.Count(); j++) {
		if (useMatch[j] >= NFEAT2CHECK)
			continue;
		imfp = dstFL.GetFeature(j);
		interp->spt[i][0] = imfp->x;
		interp->spt[i][1] = imfp->y;
		perspMap(xyMapped, perspXform, imfp->x, imfp->y);
		imfp = fmatch[j][useMatch[j]];
		xyAdj[i][0] = imfp->x - xyMapped[0];
		xyAdj[i][1] = imfp->y - xyMapped[1];
		++i;
	}
					// set interpolant resolution
	if (whres*dstFL.Height() > wvres*dstFL.Width())
		interp2_spacing(interp, .5*dstFL.Width()/whres);
	else
		interp2_spacing(interp, .5*dstFL.Height()/wvres);
	float *	wt = new float [nOK];
	for (j = wvres; j--; ) {	// use interpolant to adjust positions
		for (i = whres; i--; ) {
			int	n = interp2_weights(wt, interp,
						dstFL.Width()*i/(whres-1.),
						dstFL.Height()*j/(wvres-1.));
			while (n-- > 0) {
				warpGrid[j*whres+i][0] += wt[n]*xyAdj[n][0];
				warpGrid[j*whres+i][1] += wt[n]*xyAdj[n][1];
			}
			if (j < wvres-1 && warpGrid[(j+1)*whres+i][1] <=
						warpGrid[j*whres+i][1])
				break;
			if (i < whres-1 && warpGrid[j*whres+i+1][0] <=
						warpGrid[j*whres+i][0])
				break;
		}
		if (i >= 0) break;
	}
	interp2_free(interp);		// clean up & return
	delete [] wt;
	delete [] xyAdj;
	delete [] fmatch;
	delete [] useMatch;
	if (j >= 0) {
		DMESG(DMCdata, "Fold-over in warp grid");
		return false;
	}
	return true;
}
