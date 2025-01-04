/*
 *  pbilat.c
 *  pan
 *
 *  Fast bilateral filter modified from Chen et al., SIGGRAPH 2007
 *
 *  For floating point images, we strongly recommend creating a log
 *  image and using that as the control input.
 *
 *  Created by Greg Ward on 2/4/09.
 *  Copyright 2023 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include "dmessage.h"
#include "pimage.h"

#if defined(_WIN32) || defined(_WIN64)
#undef expf
#define expf(x)		(float)exp((double)(x))
#undef sqrtf
#define sqrtf(x)	(float)sqrt((double)(x))
#undef fabsf
#define fabsf(x)	(float)fabs((double)(x))
#define drand48()	(((double)rand())*(1./RAND_MAX))
#endif

#ifndef true
#define true	1
#define false	0
#endif

#define	SRCH_MULT	2.5f		/* multiplier for search radius */
#define INTENS_RAD	3		/* intensity sampling radius */
#define INTENS_SIGMA	(.5f*INTENS_RAD)	/* target intensity sigma */
#define MAX_LEVELS	2500		/* maximum levels for optimized ops */
#define	BOXFILT_SIGMA	1.5		/* blur tuning for bilateralFlatter() */
#define MAXBOXSAMP	5		/* maximum box filter samples */

extern double	tcos(double x);

#define  tsin(x)	tcos((x)-(M_PI/2.))

#define ctlVal(x,y)	((const float *)PpixP(ic,x,y,psiz))[prim]

/* Lookup table for negative exponent values */

#define	EXPTABLEN	1024		/* entry count for exponent lookup */
static float	expTab[EXPTABLEN];	/* exponent lookup table */
#define EXPTABMUL	-100.f		/* multiplier (roughly 1% quantization) */

static void
mkExpTab()
{
	int	i;
	
	for (i = EXPTABLEN*(expTab[0] < .5f); i--; )
		expTab[i] = expf((1.f/EXPTABMUL)*(i+.5f));
}

static float
texpf(float x)
{
	if ((x *= EXPTABMUL) >= EXPTABLEN)
		return .0f;

	return expTab[(x > 0)*(int)x];
}

/* Bilateral filter fall-back method (single channel, brute force) */
static int
bilateralFallback(ImgStruct *ib, const ImgStruct *ia, float sigma_s, float sigma_r,
				const ImgStruct *ic, int prim)
{
	const int	psiz = ImgPixelLen[ia->csp->format]*sizeof(float);
	const int	ctlSource = (ia->xres >= ib->xres);
	const float	xdens = (float)ib->xres/ia->xres;
	const float	ydens = (float)ib->yres/ia->yres;
	const int	srad = (int)(SRCH_MULT*sigma_s + .5f);
	int		xd, yd;

	for (yd = 0; yd < ib->yres; yd++) {
	     const float	ys_c = (yd+.5f)/ydens;
	     const int		ys = (int)ys_c;
	     for (xd = 0; xd < ib->xres; xd++) {
		const float	xs_c = (xd+.5f)/xdens;
		const int	xs = (int)xs_c;
		const float	cr = ctlSource ? ctlVal(xs,ys) : ctlVal(xd,yd);
		double		vsum = 0, wsum = 0;
		float		r, w;
		int		xx, yy, xrad;
		for (yy = ys-srad; yy <= ys+srad; yy++) {
		    if (yy < 0) continue;
		    if (yy >= ia->yres) break;
		    xrad = (int)(sqrtf(srad*srad - (yy-ys)*(yy-ys))+.5f);
		    for (xx = xs-xrad; xx <= xs+xrad; xx++) {
			if (xx < 0) continue;
			if (xx >= ia->xres) break;
			r = ctlSource ? ctlVal(xx,yy) :
				ctlVal((int)((xx+.5)*xdens),(int)((yy+.5)*ydens));
			w = texpf(-.5f*( ((xx+.5f-xs_c)*(xx+.5f-xs_c) +
					(yy+.5f-ys_c)*(yy+.5f-ys_c)) /
						(sigma_s*sigma_s) +
					(r-cr)*(r-cr) / (sigma_r*sigma_r) ));
			vsum += w * ((const float *)PpixP(ia,xx,yy,psiz))[prim];
			wsum += w;
		    }
		}
		if (wsum > 0) {
		    ((float *)PpixP(ib,xd,yd,psiz))[prim] = (float)(vsum/wsum);
		    continue;
		}				/* should never happen */
		((float *)PpixP(ib,xd,yd,psiz))[prim] =
				((const float *)PpixP(ia,xs,ys,psiz))[prim];
	    }
	}
	return true;
}

/* Interpolated LUT version of BOXFILT_SIGMA*sqrt(-2.*log(drand48())) */
static float
randomRadius()
{
#define	RTABSIZ	64
	static float	radTab[RTABSIZ];
	double		x = drand48();
	int		i;

	if (radTab[0] == .0f) {		/* build look-up table on first call */
		for (i = RTABSIZ; --i; )
			radTab[i] = BOXFILT_SIGMA*sqrt(-2.*log(i*(1./RTABSIZ)));
		radTab[0] = radTab[1]*radTab[1]/radTab[2];
	}
	i = x *= RTABSIZ;
	if (i >= RTABSIZ-1) i = RTABSIZ-2;
	x -= i;
	return (1.-x)*radTab[i] + x*radTab[i+1];
#undef RTABSIZ
}

/* Gaussian disk random sampling function */
static int
gaussianSamp(int xy[2], int width, int height, double sx, double sy)
{
	double	r = randomRadius();
	double	t = (2.*M_PI)*(drand48() - .5);
	double	v;

	v = width*sx + r*tcos(t);
	if ((v < 0) | (v >= width))
		return false;
	xy[0] = v;

	v = height*sy + r*tsin(t);
	if ((v < 0) | (v >= height))
		return false;
	xy[1] = v;
	
	return true;
}

/* Bilateral filter optimized for larger spatial kernels (sigma_s >= 20)
 *  when we make a stack of cheaply down-sampled images to interpolate
 *  using a fixed number of intensity samples.  Run-time does not depend
 *  on sigma_s and is also constant with sigma_r.
 */
static int
bilateralFlatter(ImgStruct *ib, const ImgStruct *ia, float sigma_s, float sigma_r,
				const ImgStruct *ic, float minmax[][2])
{
	const int	plen = ImgPixelLen[ia->csp->format];
	const int	psiz = plen*sizeof(float);
	const int	ctlSource = (ia->xres >= ib->xres);
	const float	xdens = (float)ib->xres/ia->xres;
	const float	ydens = (float)ib->yres/ia->yres;
	const int	rdu_width = (int)(1.1*BOXFILT_SIGMA*ic->xres/sigma_s + .5f);
	const int	rdu_height = (int)(1.1*BOXFILT_SIGMA*ic->yres/sigma_s + .5f);
	int		prim = plen;

    	sprintf(dmessage_buf,
			"Wide bilateral filter of %dx%d to %dx%d %d-channel image",
			ia->xres, ia->yres, ib->xres, ib->yres, plen);
	DMESG(DMCtrace, dmessage_buf);
   while (prim--) {			/* loop over primaries */
	int		nlevels = INTENS_SIGMA*(minmax[prim][1] -
						minmax[prim][0])/sigma_r;
	float		stepScale;
	ImgStruct	*sumStk, *cntStk;
	double		nsamp;
	int		xs, ys, xd, yd;
	int		i;
					/* compute needed # intensity levels */
	nlevels += !nlevels;
	stepScale = (nlevels-.03f)/(minmax[prim][1] - minmax[prim][0]);
					/* allocate box-filter image stack */
	sumStk = (ImgStruct *)calloc(2*nlevels, sizeof(ImgStruct));
	if (sumStk == NULL) {
		DMESG(DMCmemory, "calloc() failed in bilateralFlatter()");
		return false;
	}
	cntStk = sumStk + nlevels;	/* (piggy-back allocation) */
	for (i = 0; i < nlevels; i++) {
		sumStk[i].xres = rdu_width; sumStk[i].yres = rdu_height;
		sumStk[i].csp = &ICS_Y;
		cntStk[i].xres = rdu_width; cntStk[i].yres = rdu_height;
		cntStk[i].csp = &ICS_Y16;
		if (!PsetImage(sumStk+i, Pblack) ||
				!PsetImage(cntStk+i, Pblack)) {
			PfreeImage(sumStk+i);
			while (i--) {
				PfreeImage(sumStk+i);
				PfreeImage(cntStk+i);
			}
			free(sumStk);
			return false;
		}
	}
					/* Gaussian/box-filter image stack */
	nsamp = 32000.*(rdu_width*rdu_height)/(ia->xres*ia->yres);
	if (nsamp > MAXBOXSAMP)
		nsamp = MAXBOXSAMP;
	for (ys = 0; ys < ia->yres; ys++) {
		const float	*fp = (const float *)ProwPtr(ia, ys) + prim;
		for (xs = 0; xs < ia->xres; xs++, fp += plen) {
			double	ns;
			i = (int)(((ctlSource ? ctlVal(xs,ys) :
				ctlVal((int)((xs+.5f)*xdens),(int)((ys+.5f)*ydens)))
				- minmax[prim][0]) * stepScale);
			for (ns = nsamp; ns > .001; ns -= 1.) {
				int	xy[2];
				if (ns < .999 && ns <= drand48())
					break;
				if (!gaussianSamp(xy, rdu_width, rdu_height,
						(xs+.5)/ia->xres, (ys+.5)/ia->yres))
					continue;
				((float *)ProwPtr(sumStk+i, xy[1]))[xy[0]] += *fp;
				((unsigned short *)ProwPtr(cntStk+i, xy[1]))[xy[0]]++;
			}
		}
	}
					/* interpolate level planes */
	for (yd = 0; yd < ib->yres; yd++) {
		float	vpos = (yd + .5f)*rdu_height/ib->yres - .5f;
		int	vndx = vpos;
		float	*dstp = (float *)ProwPtr(ib, yd) + prim;
		if (vndx >= rdu_height-1) vndx = rdu_height-2;
		vpos -= vndx;
		ys = (int)((yd+.5f)/ydens);
		for (xd = 0; xd < ib->xres; xd++, dstp += plen) {
			float	hpos = (xd + .5f)*rdu_width/ib->xres - .5f;
			int	hndx = hpos;
			float	vsum = 0;
			float	wtsum = 0;
			float	wt;
			float	lvlf;
			int	lvl0;
			if (hndx >= rdu_width-1) hndx = rdu_width-2;
			hpos -= hndx;
			xs = (int)((xd+.5f)/xdens);
					/* sum over levels */
			lvlf = ctlSource ? ctlVal(xs,ys) : ctlVal(xd,yd);
			lvl0 = lvlf = stepScale*(lvlf - minmax[prim][0]);
			lvlf -= .5f;
			for (i = lvl0-INTENS_RAD; i <= lvl0+INTENS_RAD; i++) {
				const float		*fp;
				const unsigned short	*sp;
				float			lwt;
				if (i < 0) continue;
				if (i >= nlevels) break;
				lwt = texpf(-(i-lvlf)*(i-lvlf) *
					(.5f/INTENS_SIGMA/INTENS_SIGMA));
				fp = (const float *)ProwPtr(sumStk+i,
								vndx) + hndx;
				sp = (const unsigned short *)ProwPtr(cntStk+i,
								vndx) + hndx;
				wt = lwt*(1.f-hpos)*(1.f-vpos);
				vsum += fp[0] * wt;
				wtsum += sp[0] * wt;
				wt = lwt*hpos*(1.f-vpos);
				vsum += fp[1] * wt;
				wtsum += sp[1] * wt;
				fp = (const float *)ProwPtr(sumStk+i,
								vndx+1) + hndx;
				sp = (const unsigned short *)ProwPtr(cntStk+i,
								vndx+1) + hndx;
				wt = lwt*(1.f-hpos)*vpos;
				vsum += fp[0] * wt;
				wtsum += sp[0] * wt;
				wt = lwt*hpos*vpos;
				vsum += fp[1] * wt;
				wtsum += sp[1] * wt;
			}
			if (wtsum > 0) {
				*dstp = vsum / wtsum;
				continue;
			}		/* should never happen */
			*dstp = ((const float *)PpixP(ia,xs,ys,psiz))[prim];
		}
	}
	for (i = nlevels; i--; ) {	/* free box-filter stack */
		PfreeImage(sumStk+i);
		PfreeImage(cntStk+i);
	}
	free(sumStk);
	/* free(cntStk);	(cntStk piggy-backs allocation) */
    }
	return true;
}

/* Calculate bilateral filter on prepared floating point images:
 *  Colorspace of ib and ia match, and image ic has same pixel type.
 *  Our approach is to compute separable Gaussian filters in the two
 *  spatial dimensions and segment the signal dimension so that only
 *  a fixed number of signal samples is needed.  This results in an
 *  O(N*c) time complexity (N pixel samples, c related to sigma_s) and
 *  O(c*k) storage (k related to signal range over sigma_r).  An even
 *  faster version is called for large values of sigma_s, which is O(N).
 */
static int
bilateralFilter(ImgStruct *ib, const ImgStruct *ia, float sigma_s, float sigma_r,
				const ImgStruct *ic)
{
	const int	plen = ImgPixelLen[ia->csp->format];
	const int	psiz = plen*sizeof(float);
	const int	ctlSource = (ia->xres >= ib->xres);
	const float	xdens = (float)ib->xres/ia->xres;
	const float	ydens = (float)ib->yres/ia->yres;
	const int	srad = (int)(SRCH_MULT*sigma_s + .5f);
	const int	swid = 2*srad + 1;
	int		prim = plen;
	float		minmax[3][2];
	float **	lvlSum=NULL;
	float		*slwt=NULL, *colwta=NULL, *cwp;
	int		xd, yd;
	int		i;
					/* get extrema */
	minmax[0][0] = minmax[1][0] = minmax[2][0] = 1.f;
	minmax[0][1] = minmax[1][1] = minmax[2][1] = -1.f;
	if (!PcomputeHisto(minmax, NULL, 0, ic, (plen==3)?PHCrgb3:PHCrgb))
		return false;
	mkExpTab();			/* make exponent table */
	for (i = 0; i < plen; i++)	/* can any primary be optimized? */
		if (minmax[i][1] - minmax[i][0] <
				(MAX_LEVELS+1.)/INTENS_SIGMA*sigma_r)
			break;
					/* is spatial blur really wide? */
	if (!i & (sigma_s >= 20.f) & (ic->xres + ic->yres > 500) &&
				INTENS_SIGMA*(minmax[0][1] - minmax[0][0]) <
					sigma_s*sigma_s*sigma_r)
		return bilateralFlatter(ib, ia, sigma_s, sigma_r, ic, minmax);
	if (i < plen) {			/* create weight arrays */
		colwta = (float *)malloc(sizeof(float)*swid*ib->xres);
		slwt = (float *)malloc(sizeof(float)*swid);
		lvlSum = (float **)malloc(swid*sizeof(float *));
		if ((slwt == NULL) | (colwta == NULL) | (lvlSum == NULL))
			goto memerr;
		cwp = colwta;		/* assign column weight array */
		for (xd = 0; xd < ib->xres; xd++) {
			const float	xsc = (xd+.5f)/xdens;
			const int	xs0 = (int)xsc;
			for (i = -srad; i <= srad; i++) {
				const float	rx = xs0+i+.5f - xsc;
				*cwp++ = texpf(-rx*rx/(2.f*sigma_s*sigma_s));
			}
		}
	}
    	sprintf(dmessage_buf,
			"Bilateral filter of %dx%d to %dx%d %d-channel image",
			ia->xres, ia->yres, ib->xres, ib->yres, plen);
	DMESG(DMCtrace, dmessage_buf);
   while (prim--) {			/* loop over primaries */
	int	nlevels = INTENS_SIGMA*(minmax[prim][1] - minmax[prim][0])/sigma_r;
	float	stepScale;
	float *	lvlSumBase;

	if (nlevels > MAX_LEVELS) {	/* check number of levels */
		sprintf(dmessage_buf,
				"Intensity sigma in channel %d would require %d levels",
				prim, nlevels);
		DMESG(DMCtrace, dmessage_buf);
		if (bilateralFallback(ib, ia, sigma_s, sigma_r, ic, prim))
			continue;
		free(slwt);
		free(colwta);
		free(lvlSum);
		return false;		/* fall-back failed */
	}
	nlevels += !nlevels;
	stepScale = (nlevels-.03f)/(minmax[prim][1] - minmax[prim][0]);
	nlevels += 2*INTENS_RAD;	/* we sample past ends */
	minmax[prim][0] -= (INTENS_RAD+.01f)/stepScale;
	minmax[prim][1] += (INTENS_RAD+.01f)/stepScale;
					/* allocate level sum array */
	lvlSumBase = (float *)malloc(swid*nlevels*2 * sizeof(float));
	if (lvlSumBase == NULL)
		goto memerr;
	for (i = swid; i--; )
		lvlSum[i] = lvlSumBase + i*nlevels*2;
					/* compute output scanlines */
	for (yd = 0; yd < ib->yres; yd++) {
		const float	ysc = (yd+.5f)/ydens;
		const int	ys0 = (int)ysc;
		float *		dstp = (float *)ProwPtr(ib, yd) + prim;
		float *		lvlp;
		float		wt;
		int		right_xs;
		int		xs, ys;
					/* compute row weights at this yd */
		for (i = -srad; i <= srad; i++) {
			const float	ry = ys0+i+.5f - ysc;
			slwt[i+srad] = texpf(-ry*ry/(2.f*sigma_s*sigma_s));
		}

#define lvlPtr(i,s)	(lvlSum[i] + (int)(stepScale*((s)-minmax[prim][0]))*2)
#define add2level(x,y,xx,yy)	lvlp = lvlPtr( (xx)-(x)+srad, \
					ctlSource ? ctlVal(xx,yy) : \
			ctlVal((int)((xx+.5)*xdens),(int)((yy+.5)*ydens)) ); \
				wt = slwt[(yy)-(y)+srad]; \
				lvlp[0] += wt*((const float *)PpixP(ia,xx,yy,psiz))[prim]; \
				lvlp[1] += wt

					/* clear column level sums */
		for (i = srad; i < swid; i++)
			memset(lvlSum[i], 0, 2*sizeof(float)*nlevels);
		right_xs = 0;
					/* proceed across output scanline */
		for (xd = 0; xd < ib->xres; xd++, dstp += plen) {
			const int	xs0 = (int)((xd+.5)/xdens);
			float		wtsum = 0;
			float		lvlf;
			int		lvl0;
					/* synchronize column level sums */
			while (right_xs < xs0+srad) {
				lvlp = lvlSum[0];
				memmove(lvlSum, lvlSum+1, (swid-1)*sizeof(float *));
				memset(lvlp, 0, 2*sizeof(float)*nlevels);
				lvlSum[swid-1] = lvlp;
				if (++right_xs < ia->xres)
					for (ys = ys0-srad; ys <= ys0+srad; ys++) {
						if (ys < 0) continue;
						if (ys >= ia->yres) break;
						add2level(right_xs-srad, ys0,
							  right_xs, ys);
					}
			}
			*dstp = 0;
			cwp = colwta + xd*swid;
			lvlf = ctlSource ? ctlVal(xs0,ys0) : ctlVal(xd,yd);
			lvl0 = lvlf = stepScale*(lvlf - minmax[prim][0]);
			lvlf -= lvl0 + .5f;	/* sum over columns & levels */
			for (i = -INTENS_RAD; i <= INTENS_RAD; i++) {
				float	lwt = texpf(-(i-lvlf)*(i-lvlf) *
					(.5f/INTENS_SIGMA/INTENS_SIGMA));
				for (xs = xs0-srad; xs <= xs0+srad; xs++) {
					lvlp = lvlSum[xs-xs0+srad] + (lvl0+i)*2;
					wt = lwt * cwp[xs-xs0+srad];
					*dstp += wt * lvlp[0];
					wtsum += wt * lvlp[1];
				}
			}
			if (wtsum > 0) {
				*dstp /= wtsum;
				continue;
			}		/* should never happen */
			*dstp = ((const float *)PpixP(ia,xs0,ys0,psiz))[prim];
		}
#undef add2level
#undef lvlPtr
	}
	free(lvlSumBase);		/* end of primary loop */
    }
					/* final clean up */
	free(lvlSum);
	free(slwt);
	free(colwta);
	return true;
memerr:
	DMESG(DMCmemory, "malloc() failed in bilateralFilter()");
	return false;
}
#undef ctlVal

/* Compute bilateral filter:
 *  Source image ia is modulated by control image ic, which
 *  may be the same as ia, and is assumed the same if set to NULL.
 *  The dimensions of ic must match the larger of ia or ib.
 *  Sigma_s controls the sampling radius measured in input pixels;
 *  sigma_r controls the sampling radius in the signal domain
 *  using normalized units (e.g., an 8-bit quantum is 1./256.).
 */
int
PbilatFilter(ImgStruct *ib, const ImgStruct *ia, float sigma_s, float sigma_r,
				const ImgStruct *ic)
{
	int			retOK = false;
	ImgColorSpace		myControlCS, mySourceCS;
	ImgStruct		myControl, mySource;
	ImgStruct		dstImg;
						/* check arguments */
	if ((sigma_s <= 0.4f) | (sigma_r <= 0))	/* no point in BLF */
		return PsizeImage(ib, ia, PSbest);
	if (ia == NULL || (ia->csp == NULL) | (ia->img == NULL))
		return false;
	if (ib == NULL || ib->csp == NULL || (ib->xres <= 0) | (ib->yres <= 0))
		return false;
	if (ic == NULL || ic->img == ia->img)
		ic = ia;
	else if ((ic->img == NULL) | (ic->csp == NULL))
		return false;
	if (ib->xres > ia->xres ^ ib->yres > ia->yres) {
		DMESG(DMCparameter, "Ambiguous image dimensions");
		return false;
	}
	if ((ib->xres > ia->xres) ?
			(ib->xres != ic->xres) | (ib->yres != ic->yres) :
			(ia->xres != ic->xres) | (ia->yres != ic->yres)) {
		DMESG(DMCparameter, "Mismatched image dimensions");
		return false;
	}
						/* prepare source */
	if (ia->csp->dtype != IDTfloat) {
		if (ImgPixelLen[ia->csp->format] < ImgPixelLen[ib->csp->format])
			PcopyCS(&mySourceCS, ia->csp);
		else
			PcopyCS(&mySourceCS, ib->csp);
		PrealCS(&mySourceCS);
		mySource.csp = &mySourceCS;
		mySource.img = NULL;
		if (!PmapImage(&mySource, ia, 1.f))
			goto cleanup;
		if (ic == ia && PmatchColorSpace(ic->csp, &mySourceCS,
							PICMformat|PICMgamma))
			ic = &mySource;
		ia = &mySource;
	} else
		PcopyCS(&mySourceCS, ia->csp);
						/* prepare control channel */
	if (!PmatchColorSpace(ic->csp, &mySourceCS, PICMptype)) {
		PcopyCS(&myControlCS, &mySourceCS);
		myControlCS.gamma = ic->csp->gamma;
		myControlCS.logorig = ic->csp->logorig;
		myControl.csp = &myControlCS;
		myControl.img = NULL;
		if (!PmapImage(&myControl, ic, 1.f))
			goto cleanup;
		ic = &myControl;
	}
						/* prepare destination */
	if (ImgPixelSize(ib->csp) != ImgPixelSize(&mySourceCS)) {
		dstImg.csp = &mySourceCS;
		dstImg.xres = ib->xres; dstImg.yres = ib->yres;
		dstImg.img = NULL;
		if (!PnewImage(&dstImg, .0))
			goto cleanup;
	} else {
		if (PimagesOverlap(ia, ib)) {
			DMESG(DMCparameter, "Cannot filter image onto itself");
			return false;
		}
		if (ic != ia && PimagesOverlap(ic, ib)) {
			DMESG(DMCparameter, "Cannot use destination as control");
			return false;
		}
		if (!PnewImage(ib, .0) || !PlinkImage(&dstImg, ib))
			goto cleanup;
		dstImg.csp = &mySourceCS;
	}
						/* perform actual operation */
	retOK = bilateralFilter(&dstImg, ia, sigma_s, sigma_r, ic);
cleanup:
	if (ic == &myControl)
		PfreeImage(&myControl);
	if (ia == &mySource)
		PfreeImage(&mySource);
	if (retOK)				/* final conversion */
		retOK = PmapImage(ib, &dstImg, 1.f);	/* may be no-op */
	else
		PfreeImage(ib);
	PfreeImage(&dstImg);
	return retOK;
}

/* Noise removal routine -- modifies floating point input & output images */
static int
denoiseImage(ImgStruct *ib, ImgStruct *ia, float sigma_s, const float *noise_c)
{
	const int	nc = ImgPixelLen[ia->csp->format];
	float		srb[3], two_m[3];
	int		x, y, c;
	float		*fp;
						/* loop constants we'll need */
	for (c = nc; c--; ) {
		srb[c] = sqrtf(noise_c[2*c]);
		two_m[c] = 2.f/noise_c[2*c+1];
	}
						/* input -> conformal space */
	for (y = 0; y < ia->yres; y++) {
		fp = (float *)ProwPtr(ia, y);
		for (x = 0; x < ia->xres; x++)
			for (c = 0; c < nc; c++, fp++)
				if (*fp >= 0)
					*fp = two_m[c]*(sqrtf(noise_c[2*c] +
						*fp*noise_c[2*c+1]) - srb[c]);
				else
					*fp = -two_m[c]*(sqrtf(noise_c[2*c] -
						*fp*noise_c[2*c+1]) - srb[c]);
	}
						/* perform bilateral filter */
	if (!bilateralFilter(ib, ia, sigma_s, 1.f, ia))
		return false;
						/* output -> original space */
	for (y = 0; y < ib->yres; y++) {
		fp = (float *)ProwPtr(ib, y);
		for (x = 0; x < ib->xres; x++)
			for (c = 0; c < nc; c++, fp++)
				*fp *= fabsf(*fp)*.25f*noise_c[2*c+1] + srb[c];
	}
	return true;
}

/* Remove noise from image:
 *  Apply bilateral filter to remove noise from the given image (ia)
 *  which may be the same as destination image (ib).
 *  Input and output image dimensions must match in any case.
 *  The filter sampling radius in pixels is specified as sigma_s.
 *  A linear noise model is specified in noise_c, which is dimensioned
 *  noise_c[nc*2], where nc is the number of components in ia.
 *  The variance due to noise for component c at
 *  level v equals (noise_c[2*c] + v*noise_c[2*c+1]).
 *  A non-float image is considered to have a 0-1 range.
 */
int
PdenoiseImage(ImgStruct *ib, const ImgStruct *ia, float sigma_s,
			const float *noise_c)
{
	int		retOK = false;
	ImgColorSpace	workCS;
	ImgStruct	srcImg, dstImg;
	int		c;
						/* check arguments */
	if (ia == NULL || (ia->csp == NULL) | (ia->img == NULL))
		return false;
	if (ib == NULL || ib->csp == NULL)
		return false;
	if (ib->img == NULL) {
		ib->xres = ia->xres;
		ib->yres = ia->yres;
	} else if ((ib->xres != ia->xres) | (ib->yres != ia->yres)) {
		DMESG(DMCparameter, "Mismatched image dimensions");
		return false;
	}
	if ((c = ImgPixelLen[ia->csp->format]) > 3) {
		DMESG(DMCparameter, "Unsupported color space for denoising");
		return false;
	}
	c *= 2;
	while (c--)
		if (c&1 ? (noise_c[c] <= .0f) : (noise_c[c] < .0f)) {
			DMESG(DMCparameter, "Illegal noise model parameters");
			return false;
		}
	if (sigma_s < .7f) {
		DMESG(DMCparameter, "Noise smoothing radius too small");
		return false;
	}
	PcopyCS(&workCS, ia->csp);		/* create working image */
	PrealCS(&workCS);
	srcImg.csp = &workCS;
	srcImg.img = NULL;
	if (!PmapImage(&srcImg, ia, 1.f))	/* may just be image copy */
		goto cleanup;
						/* prepare destination */
	if (ImgPixelSize(ib->csp) != ImgPixelSize(&workCS)) {
		dstImg.csp = &workCS;
		dstImg.xres = ib->xres; dstImg.yres = ib->yres;
		dstImg.img = NULL;
		if (!PnewImage(&dstImg, .0))
			goto cleanup;
	} else {
		if (!PnewImage(ib, .0) || !PlinkImage(&dstImg, ib))
			goto cleanup;
		dstImg.csp = &workCS;
	}
						/* remove noise */
	retOK = denoiseImage(&dstImg, &srcImg, sigma_s, noise_c);
cleanup:
	PfreeImage(&srcImg);
	if (retOK)
		retOK = PmapImage(ib, &dstImg, 1.f);
	else
		PfreeImage(ib);
	PfreeImage(&dstImg);
	return retOK;
}
