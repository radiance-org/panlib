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

// BRIEF sampling pair
static struct BriefSamp {
	short		xo1, yo1;
	short		xo2, yo2;
}	ubrief[BRIEF_BITS];		// U-BRIEF sampling pattern


static void
markPix(const PanImage &im, const ABitMap2 &map, const char *fname)
{
	const PixelVal	redP = PsRGB(255,0,0);
	PanImage	rgbIm(ICS_sRGB);
	int		x, y;
	if (!rgbIm.Load(im))
		return;
	for (x = y = 0; map.Find(&x, &y); x++)
		rgbIm.SetPixel(x, y, redP);
	rgbIm.Write(fname);
}
	
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

// check that patch offset is in out of range
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
		const float	*ip0 = im.GetRowFloat(y-dy);
		const float	*ip1 = im.GetRowFloat(y+dy);
		float		*dp = ib->RowFloat(y);
		for (x = dx; x < im.Width()-dx; x++, ip0++, ip1++)
			*dp++ = .5f*(ip1[dx] - ip0[-dx]);
	}
	return true;
}

// identify corners using Harris detector (internal)
static bool
cornerDetect(ABitMap2 *fmap, const PanImage &imxx, const PanImage &imyy,
		const PanImage &imxy, const float kappa, const float thresh)
{
	if (!fmap->NewBitMap(imxy.Width(), imxy.Height()))
		return false;

	for (int y = imxy.Height(); y--; ) {
		const float *	pxx = imxx.GetRowFloat(y);
		const float *	pyy = imyy.GetRowFloat(y);
		const float *	pxy = imxy.GetRowFloat(y);
		for (int x = 0; x < imxy.Width(); x++, pxx++, pyy++, pxy++)
			if (*pxx * *pyy - *pxy * *pxy -
					kappa*(*pxx + *pyy)*(*pxx + *pyy) > thresh)
				fmap->Set(x, y);
	}
	return true;
}

// count bits in rectangle (internal)
static inline int
countRect(const ABitMap2 &bm, const int x0, const int y0, const int w, const int h)
{
	int	cnt = 0;

	for (int y = y0+h; y-- > y0; )
		for (int x = x0+w; x-- > x0; )
			cnt += bm.Check(x, y);
	return cnt;
}

// reduce clusters to just their centers (internal)
static void
clusterCenters(ABitMap2 *fmap)
{
	ABitMap2	omap = *fmap;
	int		x, y;
					// remark cluster centers
	fmap->ClearBitMap();
	for (x = y = 0; omap.Find(&x, &y); x++) {
		int	cnt = 1, siz = 1;
		do {
			const int	nsiz = siz + (siz>>3) + 1;
			const int	linc = countRect(omap, x-(nsiz-siz), y, nsiz-siz, nsiz);
			const int	rinc = countRect(omap, x+siz, y, nsiz-siz, nsiz);
			int		ncnt = cnt + countRect(omap, x, y+siz, siz, nsiz-siz);
			ncnt += (linc > rinc) ? linc : rinc;
			if (ncnt == cnt || ncnt < nsiz*nsiz>>2)
				break;	// getting too sparse
			if (linc > rinc) {
				x -= nsiz - siz;
				x *= (x > 0);
			}
			siz = nsiz;
			cnt = ncnt;
		} while (siz < BRIEF_DIA);
		omap.ClearRect(x, y, siz, siz);
		fmap->Set(x+(siz>>1), y+(siz>>1));
	}
}

// Compute feature (corner) bitmap
bool
PdetectFeatures(ABitMap2 *fmap, const PanImage &im)
{
	const float	blur_kernel = 1.f;
	const float	kappa = 0.05f;
	const float	thresh = 0.0003f;
	PanImage	imx, imy, imxy(ICS_Y);
					// argument check
	if (!(im.Ready() & PICread)) {
		DMESG(DMCparameter, "Uninitialized image in PdetectFeatures");
		return false;
	}
					// convert to Y & normalize
	if (!imxy.Load(im))
		return false;
	imxy /= imxy.GetAverage().v.f[0];
					// compute derivative images
	if (!derivatImage(&imx, imxy, 1, 0))
		return false;
	if (!derivatImage(&imy, imxy, 0, 1))
		return false;
	imxy = imx;			// square derivatives
	imxy *= imy;
	imx *= imx;
	imy *= imy;
	imxy ^= blur_kernel;		// blur them (can use box filter)
	imx ^= blur_kernel;
	imy ^= blur_kernel;
					// detect corners
	if (!cornerDetect(fmap, imx, imy, imxy, kappa, thresh))
		return false;
markPix(im, *fmap, "/tmp/cornerMap.tif");
					// eliminate loners
	fmap->Expand(-1.2);
					// reduce map (center each cluster)
	clusterCenters(fmap);
markPix(im, *fmap, "/tmp/featureMap.tif");
	return true;
}

// Identify features in the given image, return count
int
PanFeatureList::FindFeatures(const PanImage &im, float blur)
{
	PanImage	imY;
	ABitMap2	featMap;
					// convert image to Y color space
	Clear();
	if (blur > 0.5f || !PmatchColorSpace(im.GetCS(), &ICS_Y, PICMall)
			|| !imY.LinkRO(im)) {
		imY.Init(ICS_Y);
		if (!imY.Load(im))
			return -1;
	}
					// detect features
	if (!PdetectFeatures(&featMap, imY))
		return -1;
					// elide positions too near borders
	xres = imY.Width();
	yres = imY.Height();
	featMap.ClearRect(0, 0, xres, BRIEF_DIA/2);
	featMap.ClearRect(0, yres-BRIEF_DIA/2, xres, BRIEF_DIA/2);
	featMap.ClearRect(0, BRIEF_DIA/2, BRIEF_DIA/2, yres-BRIEF_DIA);
	featMap.ClearRect(xres-BRIEF_DIA/2, BRIEF_DIA/2, BRIEF_DIA/2, yres-BRIEF_DIA);
					// count up what's left
	if (!(nfeats = featMap.SumTotal()))
		return 0;
	if (blur > 0.5f)		// apply smoothing if requested
		imY ^= blur;
					// allocate & compute features
	feat = new ImgFeature [nfeats];
	if (!ubrief[0].xo1 && !ubrief[0].yo1 && !ubrief[0].xo2 && !ubrief[0].yo2)
		initBRIEF();
	int	n, x, y;		// evaluate U-BRIEF descriptors
	for (n = x = y = 0; featMap.Find(&x, &y); x++, n++) {
		feat[n].x = x; feat[n].y = y;
		for (int i = BRIEF_BITS; i--; ) {
			const BriefSamp	&bs = ubrief[i];
			if (imY.GetRowFloat(y+bs.yo1)[x+bs.xo1] <
					imY.GetRowFloat(y+bs.yo2)[x+bs.xo2])
				feat[n].descr.Set(i);
		}
	}
	return nfeats;
}

// Remove possibly confounding features from our list
int
PanFeatureList::CullSimilar(const int thresh)
{
	if (thresh >= BRIEF_BITS)
		return 0;
	ABitMap	cullMap(nfeats);
	uint32	i;
	int	j;			// identify similar descriptors
	for (i = 0; i < nfeats-1; i++)
	    for (j = i; ++j < nfeats; )
		if (PfeatureDistance(feat[i], feat[j]) <= thresh) {
			cullMap.Set(i);
			cullMap.Set(j);
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
		feat[j++] = ofeat[i++];
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
PanFeatureList::BestFeatures(const ImgFeature *farr[], int n,
				const ImgFeature &mf, const int thresh) const
{
	if (farr == NULL)
		return 0;
	memset(farr, 0, sizeof(farr[0])*n);
	if (n > nfeats)
		n = nfeats;
	if (n <= 0)
		return 0;
	int	i;
	if (n == 1) {			// special case
		int	best_di = thresh;
		for (i = nfeats; i--; ) {
			int	di = PfeatureDistance(mf, feat[i]);
			if (di >= best_di)
				continue;
			farr[0] = &feat[i];
			best_di = di;
		}
		return (best_di < thresh);
	}
					// else use qsort()
	FeatEntry	sarr[nfeats];
	int		ngood = 0;
	for (i = 0; i < nfeats; i++)
		if ((sarr[ngood].di = PfeatureDistance(mf, feat[i])) < thresh)
			sarr[ngood++].fi = i;
	if (ngood > 1)
		qsort(sarr, ngood, sizeof(FeatEntry), feat_cmp);
	if (n > ngood)
		n = ngood;
	for (i = 0; i < n; i++)
		farr[i] = &feat[sarr[i].fi];
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

#define	NFEAT2CHECK		4	// maximum # features to match

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
	float	amat[2*nOK][8];
	float	bvec[2*nOK];
	for (i = j = 0; j < fL.Count(); j++) {
		if (useMatch[j] >= NFEAT2CHECK)
			continue;
		const ImgFeature *	ffrom = fL.GetFeature(j);
		const ImgFeature *	fto = fmatch[j][useMatch[j]];
		float *	rowp = amat[2*i];
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
	return GJsolveLeastSq(8, 2*nOK, (float *)amat, bvec, xfm);
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
dist2(double x0, double y0, double x1, double y1)
{
	return (x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1);
}

// Compute warp grid (source positions) to register features
bool
PregisterWarp(float warpGrid[][2], int whres, int wvres,
		const PanFeatureList &dstFL, const PanFeatureList &srcFL,
		int thresh, PfeatureUsedMethod *cb, void *udp)
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
	const ImgFeature *	fmatch[dstFL.Count()][NFEAT2CHECK];
	uby8			useMatch[dstFL.Count()];
	double			avglength2 = 0;
	int			nOK = 0;
	const ImgFeature *	imfp;
	double			d2;
	int			i, j;
					// start by taking best matches
	memset(useMatch, 0, sizeof(useMatch[0])*dstFL.Count());
	for (i = dstFL.Count(); i--; ) {
		imfp = dstFL.GetFeature(i);
		if (srcFL.BestFeatures(fmatch[i], NFEAT2CHECK, *imfp, thresh) > 0) {
			avglength2 += dist2(xsca*imfp->x, ysca*imfp->y,
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
		d2 = dist2(xsca*imfp->x, ysca*imfp->y,
					fmatch[i][0]->x, fmatch[i][0]->y);
		while (d2*nOK > 4.*avglength2) {
			avglength2 -= d2;
			if (++useMatch[i] >= NFEAT2CHECK ||
					!fmatch[i][useMatch[i]]) {
				useMatch[i] = NFEAT2CHECK;
				--nOK;
				break;
			}
			avglength2 += d2 = dist2(xsca*imfp->x, ysca*imfp->y,
				fmatch[i][useMatch[i]]->x,
				fmatch[i][useMatch[i]]->y);
		}
	}
	if (nOK < 4) {
		DMESG(DMCdata, "Need at least 4 matching features");
		return false;
	}
	avglength2 *= 1./nOK;
	float	perspXform[8];		// fit initial perspective transform
	if (!computePerspXform(perspXform, dstFL, fmatch, useMatch)) {
		DMESG(DMCdata, "Could not solve overall perspective transform");
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
			if ((d2 = dist2(xyMapped[0], xyMapped[1],
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
		return false;
	}
	sprintf(dmessage_buf,"Computing %dx%d warp grid based on %d feature matches",
				whres, wvres, nOK);
	DMESG(DMCtrace, dmessage_buf);
	if (cb != NULL)			// report which features we're using
	    for (i = 0; i < dstFL.Count(); i++)
		for (j = 0; j < NFEAT2CHECK && fmatch[i][j]; j++) {
		    (*cb)(*dstFL.GetFeature(i), *fmatch[i][j],
					 (j==useMatch[i]), udp);
		    if (j == useMatch[i])
			break;		// skip worse matches we didn't use
		}
					// recompute mapping w/o outliers
	if (!computePerspXform(perspXform, dstFL, fmatch, useMatch)) {
		DMESG(DMCdata, "Could not resolve perspective transform");
		return false;
	}
	for (j = wvres; j--; )		// compute canonical warp grid values
		for (i = whres; i--; )
			perspMap(warpGrid[j*whres+i], perspXform,
						dstFL.Width()*i/(whres-1.),
						dstFL.Height()*j/(wvres-1.));
	if (nOK == 4)
		return true;		// not overdetermined
	if ((whres == 2) & (wvres == 2))
		return true;		// no interior to adjust
					// else realign feature positions
	INTERP2 *	interp = interp2_alloc(nOK);
	if (interp == NULL) {
		DMESG(DMCmemory, "Not enough memory for 2-D interpolant");
		return false;
	}
	float	xyAdj[nOK][2];		// skipping outliers for interpolant
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
	for (j = wvres; j--; ) {	// use interpolant to adjust positions
		for (i = whres; i--; ) {
			float	wt[nOK];
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
	if (j >= 0) {
		DMESG(DMCdata, "Fold-over in warp grid");
		return false;
	}
	return true;
}
