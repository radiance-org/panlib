/*
 *  pconvolve.c
 *  panlib
 *
 *  Apply a convolution filter.
 *
 *  Second form (PsampleImage) allows image resizing/rotation and programmable
 *  sampling pattern.
 *
 *  Created by Greg Ward on 1/8/07.
 *  Copyright 2007 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include "system.h"
#include "dmessage.h"
#include "pimage.h"

#ifndef true
#define true	1
#define false	0
#endif

/* Scale and sum image into destination with offset */
static void
sumWithOffset(ImgStruct *ib, const ImgStruct *ia,
		const int xoff, const int yoff, const float *sf)
{
	const int	nc = ImgPixelLen[ib->csp->format];
	int		x, y, i;
					/* check for zero scaling */
	for (i = nc; i--; )
		if (sf[i] != 0) break;
	if (i < 0)
		return;
	for (y = 0; y < ib->yres; y++) {
		float *		dp = (float *)ProwPtr(ib, y);
		const float *	sp = (const float *)(y-yoff < 0 ? ia->img :
				     y-yoff < ia->yres ? ProwPtr(ia,y-yoff) :
						ProwPtr(ia,ia->yres-1));
		int		xend = ib->xres;
		if (xoff < 0) {
			sp += -xoff*nc;
			xend += xoff;
		}
		for (x = 0; x < xoff; x++)
			for (i = 0; i < nc; i++)
				*dp++ += sp[i] * sf[i];
		switch (nc) {		/* optimize overlapping section */
		case 1:
			for ( ; x < xend; x++)
				*dp++ += *sp++ * *sf;
			break;
		case 3:
			for ( ; x < xend; x++) {
				*dp++ += *sp++ * sf[0];
				*dp++ += *sp++ * sf[1];
				*dp++ += *sp++ * sf[2];
			}
			break;
		default:
			for ( ; x < xend; x++)
				for (i = 0; i < nc; i++)
					*dp++ += *sp++ * sf[i];
			break;
		}
		for (sp -= nc; x < ib->xres; x++)
			for (i = 0; i < nc; i++)
				*dp++ += sp[i] * sf[i];
	}
}

/* Convolve source image with filter kernel */
int
PconvolveImage(ImgStruct *ib, const ImgStruct *ia, const ImgStruct *ikern)
{
	ImgStruct	srcImg, dstImg, krnImg;
	ImgColorSpace	mySpace;
	int		kxrad, kyrad;
	int		x, y;
						/* check parameters */
	if (ib == NULL || ib->csp == NULL)
		return false;
	if (ia == NULL || (ia->csp == NULL) | (ia->img == NULL))
		return false;
	if (ikern == NULL || (ikern->csp == NULL) | (ikern->img == NULL) |
			((ikern->xres < 3) & (ikern->yres < 3)) |
			(ikern->xres >= ia->xres) | (ikern->yres >= ia->yres))
		return false;
	if (ib->img == NULL) {
		ib->xres = ia->xres; ib->yres = ia->yres;
	} else if ((ib->xres != ia->xres) | (ib->yres != ia->yres)) {
		DMESG(DMCparameter, "Mismatched image dimensions");
		return false;
	}
						/* need linear float images */
	PcopyCS(&mySpace, ib->csp);
	PrealCS(&mySpace);
	if (!PmatchColorSpace(ia->csp, &mySpace, PICMall)) {
		srcImg.csp = &mySpace;
		srcImg.img = NULL;
		if (!PmapImage(&srcImg, ia, 1.f)) {
			PfreeImage(ib);
			return false;
		}
	} else if (!PlinkImage(&srcImg, ia)) {
		PfreeImage(ib);
		return false;
	}
	if (!PmatchColorSpace(ib->csp, &mySpace, PICMall)) {
		dstImg.csp = &mySpace;
		dstImg.xres = ia->xres; dstImg.yres = ia->yres;
		dstImg.img = NULL;
		if (!PnewImage(&dstImg, .0)) {
			PfreeImage(ib);
			PfreeImage(&srcImg);
			return false;
		}
	} else {
		if (PimagesOverlap(ib, &srcImg)) {
			DMESG(DMCparameter, "Cannot convolve image onto itself");
			PfreeImage(ib);
			PfreeImage(&srcImg);
			return false;
		}
		if (!PnewImage(ib, .0) || !PlinkImage(&dstImg, ib)) {
			PfreeImage(&srcImg);
			return false;
		}
	}
	if (!PmatchColorSpace(ikern->csp, &mySpace, PICMall)) {
		krnImg.csp = &mySpace;
		krnImg.img = NULL;
		if (!PmapImage(&krnImg, ikern, 1.f)) {
			PfreeImage(ib);
			PfreeImage(&dstImg);
			PfreeImage(&srcImg);
			return false;
		}
	} else if (!PlinkImage(&krnImg, ikern)) {
		PfreeImage(ib);
		PfreeImage(&dstImg);
		PfreeImage(&srcImg);
		return false;
	}
						/* perform convolution */
	sprintf(dmessage_buf, "Convolving %dx%d %s image with %dx%d kernel",
			ia->xres, ia->yres, PdescribeCS(ia->csp, NULL),
			ikern->xres, ikern->yres);
	DMESG(DMCtrace, dmessage_buf);
	PclearImage(&dstImg, NULL);
	kxrad = (krnImg.xres-1)/2; kyrad = (krnImg.yres-1)/2;
	for (y = 0; y < krnImg.yres; y++)
		for (x = 0; x < krnImg.xres; x++)
			sumWithOffset(&dstImg, &srcImg, kxrad-x, kyrad-y,
					(const float *)PpixPtr(&krnImg,x,y));
	PfreeImage(&srcImg);			/* clean-up */
	PfreeImage(&krnImg);
	x = PmapImage(ib, &dstImg, 1.f);	/* no-op for linked image */
	PfreeImage(&dstImg);
	return x;
}

/*
 *  Apply sampling weights & measures from given upper-left start point.
 *  The wma[] array of lists gets tiled over destination image.  The tile's
 *  array dimensions are whres by wvres, and the uppper-left of the first
 *  tile is placed at the axleft and aytop coordinates in input ia.  Steps in
 *  the input for increments in output pixel are given by ahstep and avstep.
 *  Use axleft=aytop=0 & ahstep=avstep=NULL for "onto" image mapping.
 *  Can be used to sample one channel at a time using PHCred, etc.
 *  Fill border as needed with given pixel pfv.
 */
int
PsampleImage(ImgStruct *ib, const ImgStruct *ia,
		PweightMeasure *wma[], const int whres, const int wvres,
		const double axleft, const double aytop,
		const double ahstep[2], const double avstep[2],
		PHistoChan chan, PixelVal pfv)
{
	ImgColorSpace	mySpace;
	double		myHstep[2], myVstep[2];
	ImgStruct	srcImg, dstImg;
	int		nc, x, y;
						/* check arguments */
	if (ia == NULL || (ia->img == NULL) | (ia->csp == NULL))
		return false;
	if (!wma | (whres <= 0) | (wvres <= 0))
		return false;
	if (ib == NULL || (ib->xres < whres) | (ib->yres < wvres) | (ib->csp == NULL))
		return false;
						/* PHCrgb3 is insistent */
	if ((chan != PHCluminance) & (chan != PHCrgb) &&
		    (ia->csp->format != IPFrgb) & (ib->csp->format != IPFrgb)) {
		DMESG(DMCparameter, "PsampleImage() expects RGB image");
		return false;
	}
	if (!ahstep) {
		myHstep[0] = ia->xres/(double)ib->xres;
		myHstep[1] = 0;
		ahstep = myHstep;
	}
	if (!avstep) {
		myVstep[0] = 0;
		myVstep[1] = ia->yres/(double)ib->yres;
		avstep = myVstep;
	}
	PcopyCS(&mySpace, ia->csp);		/* get real working space */
	if (chan == PHCluminance)
		mySpace.format = IPFy;
	PrealCS(&mySpace);
	pfv = PconvertPixel(pfv, &mySpace);
						/* prepare real images */
	srcImg.csp = &mySpace; srcImg.img = NULL;
	if (PmatchColorSpace(&mySpace, ia->csp, PICMall) ?
			!PlinkImage(&srcImg, ia) : !PmapImage(&srcImg, ia, 1.f))
		return false;
	dstImg.csp = &mySpace, dstImg.img = NULL;
	dstImg.xres = ib->xres; dstImg.yres = ib->yres;
	if ((ib->img != NULL && chan > PHCblue &&
				PmatchColorSpace(&mySpace, ib->csp, PICMall) &&
				!PimagesOverlap(ib, &srcImg) &&
				!PlinkImage(&dstImg, ib)) ||
				!PsetImage(&dstImg, Pblack)) {
		PfreeImage(&srcImg);
		PfreeImage(ib);
		return false;
	}
	if (ib->img == NULL && PmatchColorSpace(&mySpace, ib->csp, PICMall) &&
			!PlinkImage(ib, &dstImg)) {
		PfreeImage(&dstImg);
		PfreeImage(&srcImg);
		return false;
	}
	nc = ImgPixelLen[mySpace.format];	/* step through destination */
	for (y = 0; y < dstImg.yres; y++) {
		const int	wmyo = (y % wvres) * whres;
		float		*dp = (float *)ProwPtr(&dstImg, y);
		for (x = 0; x < dstImg.xres; x++, dp += nc) {
			const PweightMeasure	*wmp = wma[wmyo + (x % whres)];
			const int		axc = (int)(axleft + (x+.5)*ahstep[0]
							   + (y+.5)*avstep[0]);
			const int		ayc = (int)(aytop + (x+.5)*ahstep[1]
							   + (y+.5)*avstep[1]);
			const float		*sp;
			if (!wmp) continue;
			while (wmp->wt != 0) {	/* add in weighted measures */
				int	ix = axc + wmp->mx;
				int	iy = ayc + wmp->my;
				int	c = nc;
				if ((0 <= ix) & (ix < srcImg.xres) &
						(0 <= iy) & (iy < srcImg.yres))
					sp = (const float *)PpixP(&srcImg,ix,iy,
							      sizeof(float)*nc);
				else if (pfv.csp)
					sp = pfv.v.f;
				else
					c = 0;
				while (c--)
					dp[c] += wmp->wt * sp[c];
				++wmp;		/* next in this position list */
			}
		}
	}
	PfreeImage(&srcImg);			/* done with source copy */
	if (chan <= PHCblue) {			/* single channel transfer? */
		PcopyCS(&mySpace, ib->csp);
		mySpace.format = dstImg.csp->format;
		if (!PconvertColorSpace(&dstImg, &mySpace, 1.f)) {
			PfreeImage(ib);
			return false;
		}
		x = PcopyComponent(ib, (ib->csp->format!=IPFy)*chan,
				&dstImg, (dstImg.csp->format!=IPFy)*chan);
	} else					/* else convert result */
		x = PmapImage(ib, &dstImg, 1.f);

	PfreeImage(&dstImg);			/* done with copy */
	return x;
}
