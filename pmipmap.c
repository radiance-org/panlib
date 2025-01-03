/*
 *  pmipmap.c
 *
 *  MIP map generation and sampling routines.
 *
 *  Created by Greg Ward on 6/1/08.
 *  Copyright 2008 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dmessage.h"
#include "pimage.h"

#ifdef DBG
#include "imgwriter.h"
#endif

#define MIN_RES		2			/* pyramid apex resolution */

#define GOLDEN_RATIO	1.6180339887498948482	/* (1 + sqrt(5))/2 */

						/* golden spiral chirality */
#define		mmChirMax	32000
#define		mmChirRed1	(int)(mmChirMax/GOLDEN_RATIO)
#define		mmChirRed2	(int)(mmChirRed1/GOLDEN_RATIO + 1)

static const ImgRect	mmChirality[4] = {
				{mmChirMax, mmChirRed1, 0, -mmChirRed2},
				{mmChirRed2, 0, mmChirMax, mmChirRed1},
				{-mmChirRed1, -mmChirMax, mmChirRed2, 0},
				{0, -mmChirRed2, -mmChirRed1, -mmChirMax}
			};

#ifndef true
#define true	1
#define false	0
#endif

/* Spiral rectangle inwards to next lvl */
static void
mmSpiral(ImgRect *r, int lvl)
{
	const ImgRect *	cp = &mmChirality[lvl&03];
	const int	width = r->xright - r->xleft;
	const int	height = r->ybottom - r->ytop;

	r->xleft += width * cp->xleft / mmChirMax;
	r->xright += width * cp->xright / mmChirMax;
	r->ytop += height * cp->ytop / mmChirMax;
	r->ybottom += height * cp->ybottom / mmChirMax;
}

/* Compute number of levels we want for this image */
static int
mmGetNLevels(const ImgStruct *im)
{
	int	dim, lvl;

	if (im == NULL)
		return 0;
	dim = im->xres < im->yres ? im->xres : im->yres;
	for (lvl = 0; dim >= MIN_RES; lvl++)
		dim = dim * mmChirRed1 / mmChirMax;
	return lvl;
}

/* Create a MIP map for fast and general resampling */
ImgStruct *
PcreateMIPmap(const ImgStruct *im, int linkOK)
{
	int		nlevels = mmGetNLevels(im);
	ImgRect		brect, rct;
	ImgStruct *	mip;
	ImgStruct	ibase;
	int		lvl;

	if (nlevels < 2)
		return NULL;
	if (im == NULL || (im->img == NULL) | (im->csp == NULL))
		return NULL;
	mip = (ImgStruct *)calloc(nlevels+1, sizeof(ImgStruct));
	if (mip == NULL) {
		DMESG(DMCmemory, "calloc() failed in PcreateMIPmap");
		return NULL;
	}
	ibase.img = NULL;			/* create base level */
	ibase.csp = im->csp;
	ibase.yres = im->yres;			/* first step to right */
	rct.xleft = 0; rct.xright = im->xres;
	rct.ytop = 0; rct.ybottom = im->yres;
	if (linkOK) {				/* link for first part */
		PlinkImage(mip, im);
		ibase.xres = im->xres*mmChirRed1/mmChirMax;
		if (!PsetImage(&ibase, Pblack)) {
			free(mip);
			return NULL;
		}
		rct.xright = ibase.xres;
		rct.ybottom = ibase.yres*mmChirRed1/mmChirMax;
	} else {				/* else allocate our own */
		ibase.xres = im->xres + im->xres*mmChirRed1/mmChirMax;
		if (!PsetImage(&ibase, Pblack)) {
			free(mip);
			return NULL;
		}
		if (!PlinkSubimage(mip, &ibase, &rct) ||
				!PmapImage(mip, im, 1.f)) {
			PfreeImage(&ibase);
			free(mip);
			return NULL;
		}
		mmSpiral(&rct, 0);
	} 
	brect.xleft = 0; brect.xright = ibase.xres;
	brect.ytop = 0; brect.ybottom = ibase.yres;
						/* downsample pyramid */
	for (lvl = 1; lvl < nlevels; lvl++) {
		PclipRect(&rct, &brect);	/* avoids round-off problems */
		if (!PlinkSubimage(&mip[lvl], &ibase, &rct) ||
				!PsizeImage(&mip[lvl], &mip[lvl-1], PSgaussian)) {
			while (lvl--)
				PfreeImage(&mip[lvl]);
			PfreeImage(&ibase);
			free(mip);
			return NULL;
		}
		mmSpiral(&rct, lvl);		/* next pyramid level */
	}
#ifdef DBG
{
	ImgWriteBuf	iwb;
	PsetWriteBuf(&iwb, &ibase);
	(*IWInterfaceTIFF.WriteImage)("/tmp/mipmap.tif", &iwb);
}
#endif
	PfreeImage(&ibase);			/* linked everywhere else */
	return mip;
}

/* Create a floating point MIP map from any source */
ImgStruct *
PcreateFMIPmap(const ImgStruct *im)
{
	ImgColorSpace	fltCS;
	ImgStruct	fltImg;
	ImgStruct *	pim;

	if (im == NULL || (im->img == NULL) | (im->csp == NULL))
		return NULL;
	if (im->csp->dtype == IDTfloat)		/* already floating point? */
		return PcreateMIPmap(im, false);
	PcopyCS(&fltCS, im->csp);
	PrealCS(&fltCS);
	fltImg.img = NULL;
	fltImg.csp = &fltCS;
	if (!PmapImage(&fltImg, im, 1.f))	/* change color space */
		return NULL;
	pim = PcreateMIPmap(&fltImg, true);	/* make image pyramid */
	PfreeImage(&fltImg);
	return pim;
}

/* Free allocated MIP map */
void
PfreeMIPmap(ImgStruct *mip)
{
	ImgStruct *	mp2;

	if (mip == NULL)
		return;
	for (mp2 = mip; mp2->img != NULL; mp2++)
		PfreeImage(mp2);
	free(mip);
}

/* Determine which pyramid level(s) to interpolate */
static int
mmGetCoef(float *coef, const ImgStruct *mip, float rad)
{
	int	mi;
						/* find level */
	for (mi = 0; (mip[mi+2].img != NULL) &
			(rad >= (float)(0.5*GOLDEN_RATIO)); mi++)
		rad *= (float)(1./GOLDEN_RATIO);
						/* "below" coefficient */
	*coef = (0.5*GOLDEN_RATIO -  rad)*(1./(0.5*GOLDEN_RATIO - 0.5));
						/* don't extrapolate */
	if (*coef < 0)
		*coef = 0;
	else if (*coef > 1.f)
		*coef = 1.f;
						/* return "below" level */
	return mi;
}

/* Render image from MIP map (region) */
int
PrenderMIPmap(ImgStruct *ib, const ImgStruct *mip, const ImgRect *r)
{
	ImgRect			myRect;
	ImgStruct		myImg;
	int			mi;
	float			coef;

	if (ib == NULL || (ib->csp == NULL) | (ib->xres <= 0) | (ib->yres <= 0))
		return false;
	if (mip == NULL || (mip->img == NULL) | (mip->csp == NULL))
		return false;
	if (r == NULL) {
		myRect.xleft = 0; myRect.xright = mip->xres;
		myRect.ytop = 0; myRect.ybottom = mip->yres;
	} else
		myRect = *r;
						/* make sure we're allocated */
	if (!PnewImage(ib, (double)(myRect.ybottom - myRect.ytop) /
					(myRect.xright - myRect.xleft)))
		return false;
						/* choose mipmap slab */
	mi = mmGetCoef(&coef, mip, (myRect.xright - myRect.xleft)*0.5/ib->xres);
	mi += (coef < 0.25f);			/* sharpness over speed */
	if ((mi > 0) & (r != NULL)) {		/* adjust rectangle */
		myRect.xleft = r->xleft*mip[mi].xres/mip->xres;
		myRect.xright = r->xright*mip[mi].xres/mip->xres;
		myRect.ytop = r->ytop*mip[mi].yres/mip->yres;
		myRect.ybottom = r->ybottom*mip[mi].yres/mip->yres;
		r = &myRect;
	}
	if (!PlinkSubimage(&myImg, &mip[mi], r)) {
		DMESG(DMCparameter, "Bad image area in PrenderMIPmap()");
		return false;
	}
	mi = PsizeImage(ib, &myImg, PScubic);	/* interpolate chosen slab */
	PfreeImage(&myImg);			/* free linked image */
	return mi;
}

/* Sample floating point MIP map at a pixel using the given filter radius */
int
PsampleFMIPmap(float *res, const ImgStruct *mip, float x, float y, float rad)
{
	static int	nreports_left = 20;
	static int	nwarnings_left = 10;
	int		ncomp;
	int		mi, i, c;
	float		coef[2];

	if (res == NULL)
		return false;
	if (mip == NULL || (mip->img == NULL) | (mip->csp == NULL))
		return false;
	if (mip->csp->dtype != IDTfloat) {
		if (nreports_left > 0) {
			DMESG(DMCparameter, "PsampleFMIPmap() requires float");
			--nreports_left;
		}
		return false;
	}
	ncomp = ImgPixelLen[mip->csp->format];
	for (c = ncomp; c--; )			/* clear result */
		res[c] = 0;
	if ((x < 0) | (x >= mip->xres) | (y < 0) | (y >= mip->yres)) {
		if (nwarnings_left > 0) {
			DMESG(DMCwarning, "Bad sample in PsampleFMIPmap()");
			--nwarnings_left;
		}
		return false;
	}
	mi = mmGetCoef(coef, mip, rad);
	coef[1] = 1.f - coef[0];
	for (i = 2; i--; ) {			/* trilinear interpolation */
		const ImgStruct *	sip = &mip[mi+i];
		const int		rlen = sip->rowsize / sizeof(float);
		float			xv, yv;
		int			x0, y0;
		const float *		pp0;

		if (coef[i] == 0)
			continue;
		xv = sip->xres*(x + .5)/mip->xres - .5;
		x0 = xv <= 0 ? 0 : (int)xv >= sip->xres-2 ?
						sip->xres-2 : (int)xv;
		xv = xv - x0;
		yv = sip->yres*(y + .5)/mip->yres - .5;
		y0 = yv <= 0 ? 0 : (int)yv >= sip->yres-2 ?
						sip->yres-2 : (int)yv;
		yv = yv - y0;
		pp0 = (const float *)sip->img + y0*rlen + x0*ncomp;
		for (c = ncomp; c--; )
		    res[c] += coef[i]*((1.f - xv)*((1.f - yv) * pp0[c] +
						yv * pp0[rlen+c]) +
					xv*((1.f - yv) * pp0[ncomp+c] +
						yv * pp0[rlen+ncomp+c])) ;
	}
	return true;
}
