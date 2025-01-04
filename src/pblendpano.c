/*
 *  pblendpano.c
 *
 *  HF/LF image blender.
 *
 *  Created by Greg Ward on 2/16/06.
 *  Copyright 2006 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "system.h"
#include "dmessage.h"
#include "pimage.h"

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif
#if 0
#undef sqrtf
#define sqrtf(x)	(float)sqrt((double)(x))
#endif
#ifndef true
#define true	1
#define false	0
#endif

/*
 * The basic idea here is to blend the low frequencies and splice
 * the high frequencies at a prominent edge.  Finding the edge and
 * making sure the two boundaries align is most of the battle.  We
 * use a multiplicative decomposition of the blend region into low
 * and high frequencies using a Gaussian filter.  An edge convolution
 * filter is then applied to the high frequencies along the direction
 * of the blend to identify and align edges.  Starting from a single
 * matched feature point, we work our way up and down the boundary,
 * identifying an edge where we splice the high frequencies, while
 * the low frequencies are smoothly blended.  We make many calls to
 * the exorbitant pow() function, but our assumption is that the
 * blend region is small relative to the image, so we hope this is
 * OK.  In most cases, the convolution filter dominates, anyway.
 */
 
#define EDGETHRESH	0.75f			/* empirical edge threshold */

#ifndef DOWNSAMPRATE
#define DOWNSAMPRATE	16			/* downsampling ratio */
#endif
#ifndef EDGERAD
#define EDGERAD		(DOWNSAMPRATE/4)	/* edge averaging radius */
#endif
#ifndef EDGESRCH
#define EDGESRCH	(14.f*EDGERAD)		/* edge search radius */
#endif
#ifndef MAXSTRETCH
#define MAXSTRETCH	12			/* maximum stretch (%) */
#endif

#define EDGEFILTLEN	(int)(3*EDGERAD+1)

static float		edgeFilt[EDGEFILTLEN];

					/* blending flags */
#define BLPROD		0x1			/* multiply into output */
#define BLREV		0x2			/* reverse blend function */
#define BLSQRT		0x4			/* square root blending */

/* Separate image into low and high frequency versions */
static int
sepFreq(ImgStruct *hfi, ImgStruct *lfi, const ImgStruct *src)
{
	ImgStruct	tmpImg;
	int		x, y;
						/* copy/map source */
	hfi->xres = lfi->xres = src->xres;
	hfi->yres = lfi->yres = src->yres;
	if (!PmapImage(hfi, src, 1.f))
		return false;
						/* resample to get LF */
	tmpImg.csp = hfi->csp;
	tmpImg.xres = hfi->xres / DOWNSAMPRATE;
	tmpImg.yres = hfi->yres / DOWNSAMPRATE;
	tmpImg.img = NULL;
	if (!PsizeImage(&tmpImg, hfi, PSgaussian) ||
			!PsizeImage(lfi, &tmpImg, PSlinear)) {
		PfreeImage(&tmpImg);
		return false;
	}
	PfreeImage(&tmpImg);
						/* divide to get HF */
	DASSERT(hfi->csp->dtype == IDTfloat);
	DASSERT(lfi->csp->dtype == IDTfloat);
	for (y = hfi->yres; y-- > 0; ) {
		const float *	sp = (const float *)ProwPtr(lfi,y);
		float *		dp = (float *)ProwPtr(hfi,y);
		for (x = ImgPixelLen[hfi->csp->format]*hfi->xres; x-- > 0; )
			*dp++ /= *sp++;
	}
	return true;
}

/* Compute edge function and check if we have one above threshold */
static int
compEdge(float *edp, const float *inp, const int step, int len)
{
	int	hasEdge = 0;
	float	v, sum;
	int	i, j;

	if (edgeFilt[(int)EDGERAD] < 0.01f) {
		i = EDGEFILTLEN;		/* initialize edge filter */
		while (i--)
			edgeFilt[i] = i * exp(i*i*(-.5/EDGERAD/EDGERAD));
	}
						/* compute convolution */
	memset(edp, 0, sizeof(float)*len);
	for (i = EDGEFILTLEN-1; i < len-(EDGEFILTLEN-1); i++) {
		sum = 0;
		for (j = 1; j < EDGEFILTLEN; j++) {
			sum += v = edgeFilt[j]*inp[(i-j)*step];
			edp[i] += v;
			sum += v = edgeFilt[j]*inp[(i+j)*step];
			edp[i] -= v;
		}
		if (sum <= 0) continue;
		edp[i] *= (EDGEFILTLEN-1.f)/sum;
		hasEdge += (edp[i]*edp[i] > EDGETHRESH*EDGETHRESH);
	}
	return hasEdge;
}

/* Linearly interpolate float vector (color value) */
static void
linInterp(float *rv, float x, const float *v0, const float *v1, int nc)
{
	while (nc--)
		*rv++ = (1.f-x) * *v0++ + x * *v1++;
}

/* Blend and warp along a scanline (performs work of blendLine3) */
static void
blendInto(float *dstp, int flags, const float *inp, const int ncomp,
		const int step, int len, int stretch, int splice)
{
	double	bc;
	float	val[3], x;
	int	i, j;
	int	beg, end;
						/* check for splice */
	if (splice < 0) {
		beg = 0; end = len;
	} else if (flags & BLREV) {
		beg = splice; end = len;
	} else {
		beg = 0; end = splice;
	}
	if (flags & BLPROD)			/* prologue */
		dstp += beg*step;
	else
		for (i = beg; i-- > 0; dstp += step-ncomp)
			for (j = ncomp; j--; )
				*dstp++ = 1.f;
						/* compute each output pixel */
	for (i = beg; i < end; i++, dstp += step) {
		if (flags & BLREV) {
			x = (len-1)*(1.f - (len-1-(float)i)/(len-1+stretch));
			bc = M_PI*(double)i/(len-1.) - M_PI/2.;
		} else {
			x = (len-1)*(float)i/(len-1+stretch);
			bc = M_PI/2. - M_PI*(double)i/(len-1.);
		}
		if (x < 0) x = j = 0;
		else if (x > len-1) { x = len-1; j = len-2; }
		else if ((j = (int)x) >= len-1) --j;
		x -= (float)j;			/* stretch line */
		linInterp(val, x, inp+j*step, inp+(j+1)*step, ncomp);
		if (splice >= 0)		/* no weighting of splice */
			switch (flags & (BLPROD|BLSQRT)) {
			case 0:
				for (j = 0; j < ncomp; j++)
					dstp[j] = val[j];
				continue;
			case BLPROD:
				for (j = 0; j < ncomp; j++)
					dstp[j] *= val[j];
				continue;
			case BLSQRT:
				for (j = 0; j < ncomp; j++)
					dstp[j] = sqrtf(val[j]);
				continue;
			case BLPROD|BLSQRT:
				for (j = 0; j < ncomp; j++)
					dstp[j] *= sqrtf(val[j]);
				continue;
			}
		bc = 0.5 + 0.5*sin(bc);		/* else compute weight */
		if (flags & BLSQRT)
			bc *= 0.5;
		if (flags & BLPROD)
			for (j = 0; j < ncomp; j++)
				dstp[j] *= val[j] <= 0 ? 0. : pow(val[j], bc);
		else
			for (j = 0; j < ncomp; j++)
				dstp[j] = val[j] <= 0 ? 0. : pow(val[j], bc);
	}
	if (flags & BLPROD)
		return;
	for (i = len-end; i-- > 0; dstp += step-ncomp)	/* epilogue */
		for (j = ncomp; j--; )
			*dstp++ = 1.f;
}

/* Process next scanline in blend region (only works for 3 component scans) */
static int
blendLine3(float *dstp, int *edge, int *stretch,
		const float *begHF, const float *begLF,
		const float *endHF, const float *endLF,
		const int step, int len)
{	
	float	*begEF, *endEF, *comEF, *weightF;
	float	bestSum;
	int	i, j;
	int	newEdge, newStretch;
						/* compute edge functions */
	begEF = (float *)malloc(sizeof(float)*4*len);
	if (begEF == NULL) {
		DMESG(DMCmemory, "malloc() failed");
		return false;
	}
	endEF = begEF + len;
	comEF = endEF + len;
	weightF = comEF + len;
	if (!compEdge(begEF, begHF+1, step, len) |
			!compEdge(endEF, endHF+1, step, len)) {
		free(begEF);			/* smooth blend, same stretch */
		blendInto(dstp, 0, begLF, 3, step, len, *stretch, -1);
		blendInto(dstp, BLPROD|BLREV, endLF, 3, step, len, *stretch, -1);
		if (*edge < 0) {
			blendInto(dstp, BLPROD, begHF, 3, step, len, *stretch, -1);
			blendInto(dstp, BLPROD|BLREV, endHF, 3, step, len, *stretch, -1);
		} else {			/* soften change */
			blendInto(dstp, BLPROD|BLSQRT, begHF, 3, step, len, 
						*stretch, *edge);
			blendInto(dstp, BLPROD|BLSQRT|BLREV, endHF, 3, step, len,
						*stretch, *edge);
			blendInto(dstp, BLPROD|BLSQRT, begHF, 3, step, len,
						*stretch, -1);
			blendInto(dstp, BLPROD|BLSQRT|BLREV, endHF, 3, step, len,
						*stretch, -1);
		}
		*edge = -1;			/* we lost our edge */
		return true;
	}
						/* composite weight function */
	for (j = len; j-- > 0; )
		weightF[j] = .5 + .5*sin(2.*M_PI*(j+.5)/len - M_PI/2.);
	for (j = *edge > 0 ? len : 0; j-- > 0; )
		weightF[j] /= 1.f + (j - *edge)*(j - *edge)*(.5f/EDGESRCH/EDGESRCH);
						/* adjust stretch & edge */
	newStretch = *stretch;
	newEdge = *edge;
	bestSum = 0;
	for (i = -1; i <= 1; i++) {
		float	sum = 0;
		if (i < 0 && *stretch + i < -MAXSTRETCH*len/100)
			continue;
		if (i > 0 && *stretch + i > MAXSTRETCH*len/100)
			break;
		blendInto(comEF, 0, begEF, 1, 1, len, *stretch + i, len);
		blendInto(comEF, BLPROD|BLREV, endEF, 1, 1, len, *stretch + i, 0);
		for (j = len; j-- > 0; )
			sum += (comEF[j] *= weightF[j]);
		if (sum > bestSum) {
			newStretch = *stretch + i;
			bestSum = sum;
			sum = 0;
			for (j = len; j-- > 0; )
				if (comEF[j] > sum) {
					newEdge = j;
					sum = comEF[j];
				}
		}
	}
	free(begEF);
						/* blend low frequencies */
	blendInto(dstp, 0, begLF, 3, step, len, newStretch, -1);
	blendInto(dstp, BLPROD|BLREV, endLF, 3, step, len, newStretch, -1);
						/* splice high frequencies */
	if (abs(*edge - newEdge) <= (int)EDGESRCH) {
		blendInto(dstp, BLPROD, begHF, 3, step, len, newStretch, newEdge);
		blendInto(dstp, BLPROD|BLREV, endHF, 3, step, len, newStretch, newEdge);
	} else {				/* soften change */
		blendInto(dstp, BLPROD|BLSQRT, begHF, 3, step, len,
					newStretch, *edge);
		blendInto(dstp, BLPROD|BLSQRT|BLREV, endHF, 3, step, len,
					newStretch, *edge);
		blendInto(dstp, BLPROD|BLSQRT, begHF, 3, step, len,
					newStretch, newEdge);
		blendInto(dstp, BLPROD|BLSQRT|BLREV, endHF, 3, step, len,
					newStretch, newEdge);
	}
	*stretch = newStretch;			/* return new values */
	*edge = newEdge;
#if 0
dstp[step*newEdge] = 1.f;
dstp[step*newEdge+1] = 1000.f;
dstp[step*newEdge+2] = 1.f;
#endif
	return true;
}

/* Blend sections of a panorama:
 *	Input is two overlapping images in matching color spaces.
 *	The rectangle "anch" has its left-top coordinate
 *	set to a feature in "inpa" that matches "inpb" when the
 *	right-bottom corner of "inpb" is at the right-bottom corner of
 *	"anch".  (This is not necessarily a legal rectangle.)
 *	On completion, PblendPano() creates the image "blnd" to cover
 *	the overlapping region between "inpa" and "inpb".
 */
int
PblendPano(ImgStruct *blnd, const ImgRect *anch,
			const ImgStruct *inpa, const ImgStruct *inpb)
{
	ImgStruct	lfa, lfb, hfa, hfb;
	ImgStruct	subImg, iblnd;
	ImgRect		bra, brb;
						/* check arguments */
	if (blnd == NULL || blnd->csp == NULL)
		return false;
	if (anch == NULL)
		return false;
	if (inpa == NULL || (inpa->img == NULL) | (inpa->csp == NULL))
		return false;
	if (inpb == NULL || (inpb->img == NULL) | (inpb->csp == NULL))
		return false;
						/* compute blend regions */
	bra.xleft = anch->xright - inpb->xres;
	if (bra.xleft < 0) bra.xleft = 0;
	bra.xright = anch->xright;
	if (bra.xright > inpa->xres) bra.xright = inpa->xres;
	bra.ytop = anch->ybottom - inpb->yres;
	if (bra.ytop < 0) bra.ytop = 0;
	bra.ybottom = anch->ybottom;
	if (bra.ybottom > inpa->yres) bra.ybottom = inpa->yres;
	brb = bra;
	brb.xleft += inpb->xres - anch->xright;
	brb.xright += inpb->xres - anch->xright;
	brb.ytop += inpb->yres - anch->ybottom;
	brb.ybottom += inpb->yres - anch->ybottom;
						/* create HF & LF sources */
	iblnd.csp = hfa.csp = hfb.csp = lfa.csp = lfb.csp = &ICS_XYZ;
	iblnd.img = hfa.img = hfb.img = lfa.img = lfb.img = subImg.img = NULL;
	if (!PlinkSubimage(&subImg, inpa, &bra) ||
			!sepFreq(&hfa, &lfa, &subImg))
		goto fail;
	PfreeImage(&subImg);
	if (!PlinkSubimage(&subImg, inpb, &brb) ||
			!sepFreq(&hfb, &lfb, &subImg))
		goto fail;
	PfreeImage(&subImg);
						/* allocate destination */
	iblnd.xres = hfa.xres; iblnd.yres = hfa.yres;
	if (!PnewImage(&iblnd, .0))
		goto fail;
						/* blend scanlines */
	if (abs(inpb->xres - anch->xright)*(inpa->yres + inpb->yres) >
			abs(inpb->yres - anch->ybottom)*(inpa->xres + inpb->xres)) {
		const int	xstep = 3;
		const ImgStruct	*leftHF, *leftLF, *rightHF, *rightLF;
		int		xedge, xstretch, y;
		DASSERT(xstep == ImgPixelLen[hfa.csp->format]);
		if (inpb->xres < anch->xright) {
			rightHF = &hfb; rightLF = &lfb;
			leftHF = &hfa; leftLF = &lfa;
		} else {
			rightHF = &hfa; rightLF = &lfa;
			leftHF = &hfb; leftLF = &lfb;
		}
		xedge = anch->xleft - bra.xleft; xstretch = 0;
		for (y = anch->ytop - bra.ytop; y < hfa.yres; y++)
			if (!blendLine3((float *)ProwPtr(&iblnd,y),
					&xedge, &xstretch,
					(const float *)ProwPtr(leftHF,y),
					(const float *)ProwPtr(leftLF,y),
					(const float *)ProwPtr(rightHF,y),
					(const float *)ProwPtr(rightLF,y),
					xstep, iblnd.xres))
				goto fail;
		xedge = anch->xleft - bra.xleft; xstretch = 0;
		for (y = anch->ytop - bra.ytop - 1; y >= 0; y--)
			if (!blendLine3((float *)ProwPtr(&iblnd,y),
					&xedge, &xstretch,
					(const float *)ProwPtr(leftHF,y),
					(const float *)ProwPtr(leftLF,y),
					(const float *)ProwPtr(rightHF,y),
					(const float *)ProwPtr(rightLF,y),
					xstep, iblnd.xres))
				goto fail;
	} else {				/* vertical panorama */
		const int	ystep = hfa.rowsize/sizeof(float);
		const ImgStruct	*topHF, *topLF, *bottomHF, *bottomLF;
		int		yedge, ystretch, x;
		DASSERT(ystep*sizeof(float) == hfa.rowsize);
		DASSERT((hfa.rowsize == lfa.rowsize) &
			(lfa.rowsize == hfb.rowsize) &
			(hfb.rowsize == lfb.rowsize) &
			(lfb.rowsize == iblnd.rowsize));
		if (inpb->yres < anch->ybottom) {
			topHF = &hfa; topLF = &lfa;
			bottomHF = &hfb; bottomLF = &lfb;
		} else {
			topHF = &hfb; topLF = &lfb;
			bottomHF = &hfa; bottomLF = &lfa;
		}
		yedge = anch->ytop - bra.ytop; ystretch = 0;
		for (x = anch->xleft - bra.xleft; x < hfa.xres; x++)
			if (!blendLine3((float *)PpixPtr(&iblnd,x,0),
					&yedge, &ystretch,
					(const float *)PpixPtr(topHF,x,0),
					(const float *)PpixPtr(topLF,x,0),
					(const float *)PpixPtr(bottomHF,x,0),
					(const float *)PpixPtr(bottomLF,x,0),
					ystep, iblnd.yres))
				goto fail;
		yedge = anch->ytop - bra.ytop; ystretch = 0;
		for (x = anch->xleft - bra.xleft - 1; x >= 0; x--)
			if (!blendLine3((float *)PpixPtr(&iblnd,x,0),
					&yedge, &ystretch,
					(const float *)PpixPtr(topHF,x,0),
					(const float *)PpixPtr(topLF,x,0),
					(const float *)PpixPtr(bottomHF,x,0),
					(const float *)PpixPtr(bottomLF,x,0),
					ystep, iblnd.yres))
				goto fail;
	}
						/* copy/convert result */
	if (!PmapImage(blnd, &iblnd, 1.f))
		goto fail;
						/* clean up & return */
	PfreeImage(&iblnd);
	PfreeImage(&hfa); PfreeImage(&hfb);
	PfreeImage(&lfa); PfreeImage(&lfb);
	return true;
fail:
	PfreeImage(&subImg); PfreeImage(&iblnd);
	PfreeImage(&hfa); PfreeImage(&hfb);
	PfreeImage(&lfa); PfreeImage(&lfb);
	PfreeImage(blnd);
	return false;
}
