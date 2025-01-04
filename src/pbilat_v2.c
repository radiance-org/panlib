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
 *  Copyright 2009 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "dmessage.h"
#include "pimage.h"

#if defined(_WIN32) || defined(_WIN64)
#undef expf
#define expf(x)		(float)exp((double)(x))
#undef sinhf
#define sinhf(x)	(float)sinh((double)(x))
#undef asinhf
#define asinhf(x)	(float)asinh((double)(x))
#endif

#ifndef true
#define true	1
#define false	0
#endif

#define	SIGMA_MULT	4.52f		/* multiplier for search radius */
#define INTENS_RAD	2		/* intensity sampling radius */
#define INTENS_SIGMA	(.5f*INTENS_RAD)	/* target intensity sigma */
#define MAX_LEVELS	1024		/* maximum levels supported */

/* Calculate bilateral filter on prepared floating point images:
 *  Colorspace of ib and ia match, and image ic has same pixel type.
 *  Our approach is to compute separable Gaussian filters in the two
 *  spatial dimensions and segment the signal dimension so that only
 *  a fixed number of signal samples is needed.  This results in an
 *  O(N*M) time complexity (N samples, M related to sigma_s) and
 *  O(M) storage, not counting the buffers to hold our float images.
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
	const int	srad = (int)(SIGMA_MULT*sigma_s + .5f);
	const int	swid = 2*srad + 1;
	int		prim = plen;
	float		minmax[3][2];
	unsigned int	nlevels;
	double		stepScale;
	float **	lvlSum;
	float		*colwta, *cwp, *slwt;
	int		xd, yd;
	int		i;
					/* get extrema */
	minmax[0][0] = minmax[1][0] = minmax[2][0] = 1.f;
	minmax[0][1] = minmax[1][1] = minmax[2][1] = -1.f;
	if (!PcomputeHisto(minmax, NULL, 0, ic, (plen==3)?PHCrgb3:PHCrgb))
		return false;
					/* allocate weight arrays */
	colwta = (float *)malloc(sizeof(float)*swid*ib->xres);
	slwt = (float *)malloc(sizeof(float)*swid);
	lvlSum = (float **)malloc(swid*sizeof(float *));
	if ((slwt == NULL) | (colwta == NULL) | (lvlSum == NULL))
		goto memerr;
	cwp = colwta;			/* assign column weight array */
	for (xd = 0; xd < ib->xres; xd++) {
		const float	xsc = (xd+.5f)/xdens;
		const int	xs0 = (int)xsc;
		for (i = -srad; i <= srad; i++) {
			const float	rx = xs0+i+.5f - xsc;
			*cwp++ = expf(-rx*rx/(2.f*sigma_s*sigma_s));
		}
	}
    	sprintf(dmessage_buf,
			"Bilateral filter of %dx%d to %dx%d %d-channel image",
			ia->xres, ia->yres, ib->xres, ib->yres, plen);
	DMESG(DMCtrace, dmessage_buf);
					/* loop over primaries */
    while (prim--) {
					/* compute number of levels */
	nlevels = INTENS_SIGMA*(minmax[prim][1] - minmax[prim][0])/sigma_r + .5f;
	if ((nlevels < 2) | (nlevels > MAX_LEVELS)) {
		sprintf(dmessage_buf,
				"Intensity sigma in channel %d requires %u levels",
				prim, nlevels);
		DMESG(DMCparameter, dmessage_buf);
		free(slwt);
		free(colwta);
		free(lvlSum);
		return false;
	}
	stepScale = (nlevels-.03)/(minmax[prim][1] - minmax[prim][0]);
	nlevels += 2*INTENS_RAD;	/* we sample past ends */
	minmax[prim][0] -= (INTENS_RAD+.01)/stepScale;
	minmax[prim][1] += (INTENS_RAD+.01)/stepScale;
					/* allocate level arrays */
	for (i = 0; i < swid; i++)
		if ((lvlSum[i] = (float *)malloc(nlevels*2 *
						sizeof(float))) == NULL) {
			while (i--) free(lvlSum[i]);
			goto memerr;
		}
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
			slwt[i+srad] = expf(-ry*ry/(2.f*sigma_s*sigma_s));
		}

#define ctlVal(x,y)	((float *)PpixP(ic,x,y,psiz))[prim]
#define lvlPtr(i,s)	(lvlSum[i] + (int)(stepScale*((s)-minmax[prim][0]))*2)
#define add2level(x,y,xx,yy)	lvlp = lvlPtr( (xx)-(x)+srad, \
					ctlSource ? ctlVal(xx,yy) : \
			ctlVal((int)((xx+.5)*xdens),(int)((yy+.5)*ydens)) ); \
				wt = slwt[(yy)-(y)+srad]; \
				lvlp[0] += wt*((float *)PpixP(ia,xx,yy,psiz))[prim]; \
				lvlp[1] += wt

					/* clear column level sums */
		for (i = srad; i < swid; i++)
			memset(lvlSum[i], 0, 2*sizeof(float)*nlevels);
		right_xs = 0;
					/* proceed across output scanline */
		for (xd = 0; xd < ib->xres; xd++, dstp += plen) {
			const int	xs0 = (int)((xd+.5)/xdens);
			float		wtsum = 0;
			double		lvlf;
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
			lvlf -= lvl0 + .5;	/* sum over columns & levels */
			for (i = -INTENS_RAD; i <= INTENS_RAD; i++) {
				float	lwt = expf(-(i-lvlf)*(i-lvlf) *
					(1.f/(2.f*INTENS_SIGMA*INTENS_SIGMA)));
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
			}
					/* fill from above and left */
			if (yd > 0) {
				if (xd > 0) {
					*dstp += dstp[-ib->rowsize/sizeof(float)
							- plen] * .2f;
					*dstp += dstp[-plen] * .4f;
					wt = .4f;
				} else
					wt = 1.f;
				*dstp += wt * dstp[-ib->rowsize/sizeof(float)];
			} else if (xd > 0)
				*dstp = dstp[-plen];
		}
#undef add2level
#undef lvlPtr
#undef ctlVal
	}
	for (i = swid; i--; )		/* end of primary loop */
		free(lvlSum[i]);
    }
					/* final clean up */
	free(lvlSum);
	free(slwt);
	free(colwta);
	return true;
memerr:
	DMESG(DMCmemory, "Out of memory in bilateralFilter()");
	return false;
}

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
	if (sigma_r <= 0) {
		DMESG(DMCparameter, "Illegal intensity sigma");
		return false;
	}
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
	float		m[3], _m[3], m_b[3], b_m[3];
	int		x, y, c;
	float		*fp;
						/* loop constants we'll need */
	for (c = nc; c--; ) {
		m[c] = sqrtf(noise_c[2*c+1]);
		_m[c] = 1.f/m[c];
		m_b[c] = sqrtf(noise_c[2*c+1]/noise_c[2*c]);
		b_m[c] = 1.f/m_b[c];
	}
						/* input -> conformal space */
	for (y = 0; y < ia->yres; y++) {
		fp = (float *)ProwPtr(ia, y);
		for (x = 0; x < ia->xres; x++)
			for (c = 0; c < nc; c++, fp++)
				*fp = asinhf(*fp*m_b[c])*_m[c];
	}
						/* perform bilateral filter */
	if (!bilateralFilter(ib, ia, sigma_s, 1.f, ia))
		return false;
						/* output -> original space */
	for (y = 0; y < ib->yres; y++) {
		fp = (float *)ProwPtr(ib, y);
		for (x = 0; x < ib->xres; x++)
			for (c = 0; c < nc; c++, fp++)
				*fp = sinhf(*fp*m[c])*b_m[c];
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
		if (noise_c[c] <= .0f) {
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
