/*
 *  pidequant.c
 *  pan
 *
 *  Created by Greg Ward on 2/2/2009.
 *  Algorithm #1 by Peter Longhurst and Richard Webb, Dolby.
 *  Algorithm #2 by Greg Ward, Dolby.
 *  Copyright 2015 Anyhere Software and Dolby Laboratories.
 *  All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "system.h"
#include "dmessage.h"
#include "pimage.h"
#include "rtmath.h"			/* needed for uint16, uint32 */

#ifndef true
#define true	1
#define false	0
#endif

#define DEQ_OVER	(.75/256.)	/* allowed quanta delta */
#define DEQ_RAD		4.9f		/* image blur radius */

/* Dequantize scanline with selective filter */
static int
dequantScan(const uby8 *scan, int len, int y, void *udp)
{
	ImgStruct *	im = (ImgStruct *)udp;
	const float *	iscan = (const float *)scan;
	float *		oscan = (float *)ProwPtr(im, y);
	
	for (len *= ImgPixelLen[im->csp->format]; len--; oscan++, iscan++) {
		double	pdif = *oscan - *iscan;
		if ((pdif > DEQ_OVER) | (pdif < -DEQ_OVER))
			*oscan = *iscan;
	}
	return 0;
}

/* Byte to float image conversion with debanding filter */
int
PdequantizeImage(ImgStruct *ib, const ImgStruct *ia, float sf)
{
	const ImgColorSpace *	cstarget;
	ImgColorSpace		mySpace;
						/* check format */
	if (ia == NULL || (ia->img == NULL) | (ia->csp == NULL))
		return false;
	if (ia->csp->dtype != IDTubyte) {
		DMESG(DMCparameter, "Input image not 8-bits/channel");
		return false;
	}
	if (ib == ia)
		return false;
	if (ib == NULL || ib->csp == NULL)
		return false;
	if (ib->csp->dtype != IDTfloat) {
		DMESG(DMCparameter, "Output image not floating point");
		return false;
	}
	if (ib->img == NULL) {
		ib->xres = ia->xres; ib->yres = ia->yres;
	} else if ((ib->xres != ia->xres) | (ib->yres != ia->yres)) {
		DMESG(DMCparameter, "Mismatched image dimensions");
		return false;
	}
	cstarget = ib->csp;
	PcopyCS(&mySpace, ia->csp);		/* orig. primaries & gamma */
	mySpace.dtype = IDTfloat;
	ib->csp = &mySpace;
						/* convert to float & blur */
	if (!PblurImage(ib, ia, DEQ_RAD))
		return false;
						/* dequantization filter */
	if (PconvertImageCB(ia, ib->csp, 1.f, dequantScan, ib) < 0) {
		PfreeImage(ib);
		return false;
	}
						/* final color conversion */
	return PconvertColorSpace(ib, cstarget, sf);
}

/* Compute area average using summed area table */
static float
areaAvg(const uint32 *sat, const int xres, const int yres,
			const int xc, const int yc, int rad)
{
	while ((xc < rad) | (yc < rad) | (xc > xres-rad) | (yc > yres-rad))
		if (!--rad)
			return -1.f;

	return (sat[(yc-rad)*xres+xc-rad] + sat[(yc+rad-1)*xres+xc+rad-1] -
		sat[(yc-rad)*xres+xc+rad-1] - sat[(yc+rad-1)*xres+xc-rad]) /
		((2*rad-1)*(2*rad-1)*256.f) + (.5f/256.f);
}

/* Smooth large gradients in the specified color channel */
static int
dequantChannel(ImgStruct *ib, const ImgStruct *ia, const int c,
				const int quant, const int minrad)
{
	const int	ncomp = ImgPixelLen[ib->csp->format];
	const int	mxres = (ib->xres + minrad-1)/minrad;
	const int	myres = (ib->yres + minrad-1)/minrad;
	uby8 *		mintab = (uby8 *)malloc(myres*mxres);
	uby8 *		maxtab = (uby8 *)calloc(myres*mxres, 1);
	ImgStruct	radtab;
	uint32 *	sumtab;
	int		x, y, x2, y2;
	const uby8 *	pp;
						/* just to be sure... */
	DASSERT(ib->xres == ia->xres & ib->yres == ia->yres);
	DASSERT(ia->csp->dtype == IDTubyte & ib->csp->dtype == IDTfloat);
						/* get superpixel extrema */
	if ((mintab == NULL) | (maxtab == NULL)) {
		DMESG(DMCmemory, "Cannot allocate extrema tables");
		return false;
	}
	memset(mintab, 0xff, myres*mxres);	/* maxtab already zeroed */
	for (y = 0; y < ia->yres; y++) {
		uby8 *		mintp = mintab + (y*myres/ia->yres)*mxres;
		uby8 *		maxtp = maxtab + (y*myres/ia->yres)*mxres;
		int		mrem = minrad;
		pp = ProwPtr(ia,y) + c;
		for (x = 0; x < ia->xres; x++, pp += ncomp) {
			if (*pp < *mintp) *mintp = *pp;
			if (*pp > *maxtp) *maxtp = *pp;
			if (!--mrem) {
				++mintp; ++maxtp;
				mrem = minrad;
			}
		}
	}
	for (x = myres*mxres; x--; )		/* check if anything to do */
		if (maxtab[x] - mintab[x] <= 2*quant)
			break;
	if (x < 0) {				/* no gradient areas? */
		free(mintab); free(maxtab);
		return true;
	}
						/* compute radius table */
	radtab.xres = mxres; radtab.yres = myres;
	radtab.csp = &ICS_Y16;
	radtab.img = NULL;
	if (!PsetImage(&radtab, Pblack))
		return false;
	for (y = 0; y < myres; y++)
	    for (x = 0; x < mxres; x++) {
		const int	rmin = maxtab[y*mxres+x] - quant;
		const int	rmax = mintab[y*mxres+x] + quant;
		int		rad = 1;
		if (rmin > rmax)
			continue;
		while ((y-rad >= 0) & (y+rad < myres) &
				(x-rad >= 0) & (x+rad < mxres)) {
		    for (x2 = x-rad; x2 <= x+rad; x2++) {
			if ((mintab[(y-rad)*mxres+x2] < rmin) |
					(maxtab[(y-rad)*mxres+x2] > rmax))
				goto hitlimit;
			if ((mintab[(y+rad)*mxres+x2] < rmin) |
					(maxtab[(y+rad)*mxres+x2] > rmax))
				goto hitlimit;
		    }
		    for (y2 = y-rad+1; y2 < y+rad; y2++) {
			if ((mintab[y2*mxres+x-rad] < rmin) |
					(maxtab[y2*mxres+x-rad] > rmax))
				goto hitlimit;
			if ((mintab[y2*mxres+x+rad] < rmin) |
					(maxtab[y2*mxres+x+rad] > rmax))
				goto hitlimit;
		    }
		    ++rad;
		}
	    hitlimit:
		if (!--rad)
			continue;
		*(uint16 *)PpixPtr(&radtab,x,y) = rad*minrad;
	    }
	free(mintab); free(maxtab);		/* done with these */
	{					/* resample radius image */
		ImgStruct	itmp;
		itmp.xres = ia->xres; itmp.yres = ia->yres;
		itmp.csp = radtab.csp;
		itmp.img = NULL;
		if (!PsizeImage(&itmp, &radtab, PSlinear))
			return false;
		PfreeImage(&radtab);
		radtab = itmp;
	}
						/* compute summed area table */
	if ((sumtab = (uint32 *)malloc(sizeof(uint32)*ia->xres*ia->yres)) == NULL) {
		DMESG(DMCmemory, "Cannot allocate summed area table");
		PfreeImage(&radtab);
		return false;
	}
	for (y = 0; y < ia->yres; y++) {
		uint32 *	satp = sumtab + y*ia->xres;
		pp = ProwPtr(ia,y) + c;
		for (x = 0; x < ia->xres; x++, pp += ncomp, satp++) {
			*satp = *pp;
			if (x) {
				*satp += satp[-1];
				if (y)
					*satp += satp[-ia->xres] -
							satp[-ia->xres-1];
			} else if (y)
				*satp += satp[-ia->xres];
		}
	}
						/* apply contour filter */
	for (y = 0; y < ib->yres; y++) {
	    const uint16 *	rp = (const uint16 *)ProwPtr(&radtab,y);
	    for (x = 0; x < ib->xres; x++, rp++)
		if (*rp > 1) {
		    float	v = areaAvg(sumtab,ia->xres,ia->yres,x,y,*rp);
		    if (v >= .0f)
			((float *)PpixPtr(ib,x,y))[c] = v;
		}
	}			
	PfreeImage(&radtab);
	free(sumtab);
	return true;
}

/* Byte to float image conversion with variable-width debanding filter */
int
PdequantizeImage2(ImgStruct *ib, const ImgStruct *ia, float sf, int quant, int minrad)
{
	const ImgColorSpace *	cstarget;
	ImgColorSpace		mySpace;
	int			ncomp, c;
						/* check format */
	if (ia == NULL || (ia->img == NULL) | (ia->csp == NULL))
		return false;
	if (ia->xres*ia->yres >= 1L<<24) {
		DMESG(DMCwarning, "Image too large for PdequantizeImage2");
		return PdequantizeImage(ib, ia, sf);
	}
	if (ia->csp->dtype != IDTubyte) {
		DMESG(DMCparameter, "Input image not 8-bits/channel");
		return false;
	}
	if (ib == ia)
		return false;
	if (ib == NULL || ib->csp == NULL)
		return false;
	if (ib->csp->dtype != IDTfloat) {
		DMESG(DMCparameter, "Output image not floating point");
		return false;
	}
	ncomp = ImgPixelLen[ia->csp->format];
	if (ib->img == NULL) {
		ib->xres = ia->xres; ib->yres = ia->yres;
	} else if ((ib->xres != ia->xres) | (ib->yres != ia->yres)) {
		DMESG(DMCparameter, "Mismatched image dimensions");
		return false;
	}
						/* check other arguments */
	if (quant <= 0)
		quant = 1;
	if ((minrad <= 1) | (minrad >= ia->xres/2) | (minrad >= ia->yres/2))
		minrad = ia->yres/64 + 2*(ia->yres < 128);
						/* default conversion */
	cstarget = ib->csp;
	PcopyCS(&mySpace, ia->csp);		/* orig. primaries & gamma */
	mySpace.dtype = IDTfloat;
	ib->csp = &mySpace;
	if (!PmapImage(ib, ia, 1.f))
		return false;
	for (c = ncomp; c--; )			/* filter each channel */
		if (!dequantChannel(ib, ia, c, quant, minrad)) {
			PfreeImage(ib);
			return false;
		}
						/* final color conversion */
	return PconvertColorSpace(ib, cstarget, sf);
}
