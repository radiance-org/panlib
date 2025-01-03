/*
 *  phistomatch.c
 *  pan
 *
 *  Image histogram matching routines.
 *
 *  Created by Greg Ward on 6/9/09.
 *  Copyright 2009 Anyhere Software. All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "system.h"
#include "dmessage.h"
#include "pimage.h"
#include "tmprivat.h"

#ifndef true
#define true	1
#define false	0
#endif

/* Map histogram of one byte-based image to another */
static int
byteMapHisto(ImgStruct *ib, const ImgStruct *ia)
{
	static uby8	minmax[2] = {0, 0xff};
	const int	plen = ImgPixelLen[ia->csp->format];
	unsigned long	hista[256], tota, histb[256], totb;
	uby8		mapb2a[256];
	int		i, x, y;
	double		hcvt;
	unsigned long	n;

	if ((i = plen) > 3) i = 3;		/* map each primary */
	while (i--) {
		memset(hista, 0, sizeof(hista));
		tota = PcomputeHisto(minmax, hista, 256, ia, i);
		if (!tota)
			return false;
		memset(histb, 0, sizeof(histb));
		totb = PcomputeHisto(minmax, histb, 256, ib, i);
		if (!totb)
			return false;
		hcvt = (double)tota / totb;
		n = 0;
		mapb2a[0] = x = 0;
		for (y = 0; y < 255; mapb2a[++y] = x) {
			n += (unsigned long)(histb[y]*hcvt + .5);
			while (x < 255 && n >= hista[x])
				n -= hista[x++];
		}
		for (y = 0; y < ib->yres; y++) {
			uby8 *	p = ProwPtr(ib,y) + i;
			for (x = ib->xres; x--; p += plen)
				*p = mapb2a[*p];
		}
	}
	return true;
}

/* Map histogram of one short-based image to another */
static int
shortMapHisto(ImgStruct *ib, const ImgStruct *ia)
{
	static uint16	minmax[2] = {0, 0xffff};
	const int	plen = ImgPixelLen[ia->csp->format];
	unsigned long	hista[1<<16], tota, histb[1<<16], totb;
	uint16		mapb2a[1<<16];
	int		i, x, y;
	double		hcvt;
	unsigned long	n;

	if ((i = plen) > 3) i = 3;		/* map each primary */
	while (i--) {
		memset(hista, 0, sizeof(hista));
		tota = PcomputeHisto(minmax, hista, 1<<16, ia, i);
		if (!tota)
			return false;
		memset(histb, 0, sizeof(histb));
		totb = PcomputeHisto(minmax, histb, 1<<16, ib, i);
		if (!totb)
			return false;
		hcvt = (double)tota / totb;
		n = 0;
		mapb2a[0] = x = 0;
		for (y = 0; y < 0xffff; mapb2a[++y] = x) {
			n += (unsigned long)(histb[y]*hcvt + .5);
			while (x < 0xffff && n >= hista[x])
				n -= hista[x++];
		}
		for (y = 0; y < ib->yres; y++) {
			uint16 *	p = (uint16 *)ProwPtr(ib,y) + i;
			for (x = ib->xres; x--; p += plen)
				*p = mapb2a[*p];
		}
	}
	return true;
}

#undef	UVSCALE
#define UVSCALE		409.

/* Private struct used to accumulate transformed float image statistics */
typedef struct {
	TMstruct *	tms;			/* luminance histogram */
	unsigned long	uvhist[2][256];		/* CIE (u',v') histograms */
	TMbright *	ls;			/* encoded luminance buffer */
} LuvStatBuf;

/* Gather XYZ scanline statistics for floatMapHisto (return -1 to abort) */
static int
xyzStats(const uby8 *scan, int len, int y, void *udp)
{
	LuvStatBuf *	sbp = (LuvStatBuf *)udp;
	const float *	fp = (const float *)scan;
	int		x;
						/* add to luminance histogram */
	if (tmCvColors(sbp->tms, sbp->ls, TM_NOCHROM, (COLOR *)fp, len) != TM_E_OK)
		return -1;
	if (tmAddHisto(sbp->tms, sbp->ls, len, 1) != TM_E_OK)
		return -1;
						/* add to (u',v') histogram */
	for (x = 0; x < len; x++, fp += 3) {
		const double	den = fp[0] + 15.*fp[1] + 3.*fp[2];
		int		u, v;
		if (sbp->ls[x] < MINBRT)
			continue;
		u = (int)(UVSCALE*4.*fp[0]/den);
		v = (int)(UVSCALE*9.*fp[1]/den);
		if ((u < 0) | (u > 255) | (v < 0) | (v > 255))
			continue;
		sbp->uvhist[0][u]++;
		sbp->uvhist[1][v]++;
	}
	return len;
}

/* Gather luminance statistics for scanline (return -1 to abort) */
static int
lumStats(const uby8 *scan, int len, int y, void *udp)
{
	LuvStatBuf *	sbp = (LuvStatBuf *)udp;
						/* add to luminance histogram */
	if (tmCvLums(sbp->ls, (float *)scan, len) != TM_E_OK)
		return -1;
	if (tmAddHisto(sbp->tms, sbp->ls, len, 1) != TM_E_OK)
		return -1;
	return len;
}

/* Map histogram of any image type to float-based target */
static int
floatMapHisto(ImgStruct *ib, const ImgStruct *ia)
{
	int			ok = false;
	const ImgColorSpace *	orig_csp = ib->csp;
	LuvStatBuf		luva, luvb;
	uby8			mapuv[2][256];
	int			i, x, y, histamax, histbmax;
	double			hcvt;
	unsigned long		tota, totb, n;
	float *			fp;
						/* initialize stat struct's */
	luva.tms = tmInit(TM_F_BW|TM_F_NOSTDERR, NULL, .0);
	tmSetSpace(luva.tms, TM_XYZPRIM, 1.);
	luvb.tms = tmDup(luva.tms);
	luva.ls = (TMbright *)malloc(sizeof(TMbright)*ia->xres);
	luvb.ls = (TMbright *)malloc(sizeof(TMbright)*ib->xres);
	if ((luva.ls == NULL) | (luvb.ls == NULL)) {
		DMESG(DMCmemory, "malloc() failed in floatMapHisto");
		return false;
	}
						/* gather statistics */
	if (ib->csp->format == IPFy) {		/* luminance only */
		if (PconvertImageCB(ia, &ICS_Y, 1.f, lumStats, &luva) < 0)
			goto cleanup;
		if (PconvertImageCB(ib, &ICS_Y, 1.f, lumStats, &luvb) < 0)
			goto cleanup;
	} else {				/* separate chrominance */
		memset(luva.uvhist, 0, sizeof(luva.uvhist));
		if (PconvertImageCB(ia, &ICS_XYZ, 1.f, xyzStats, &luva) < 0)
			goto cleanup;
		if (PconvertColorSpace(ib, &ICS_XYZ, 1.f) < 0)
			goto cleanup;
		memset(luvb.uvhist, 0, sizeof(luvb.uvhist));
		if (PconvertImageCB(ib, &ICS_XYZ, 1.f, xyzStats, &luvb) < 0)
			goto cleanup;
		for (i = 2; i--; ) {		/* compute u & v mappings */
			tota = totb = 0;
			for (y = 256; y--; ) {
				tota += luva.uvhist[i][y];
				totb += luvb.uvhist[i][y];
			}
			hcvt = (double)tota / totb;
			n = 0;
			mapuv[i][0] = x = 0;
			for (y = 0; y < 255; mapuv[i][++y] = x) {
				n += (unsigned long)(luvb.uvhist[i][y]*hcvt + .5);
				while (x < 255 && n >= luva.uvhist[i][x])
					n -= luva.uvhist[i][x++];
			}
		}
	}
						/* compute luminance mapping */
	tota = 0;
	y = HISTI(luva.tms->hbrmax) + 1 - HISTI(luva.tms->hbrmin);
	while (y--) tota += luva.tms->histo[y];
	totb = 0;
	y = HISTI(luvb.tms->hbrmax) + 1 - HISTI(luvb.tms->hbrmin);
	while (y--) totb += luvb.tms->histo[y];
	hcvt = (double)tota / totb;
	luvb.tms->mbrmin = luvb.tms->hbrmin;
	luvb.tms->mbrmax = luvb.tms->hbrmax;
	luvb.tms->lumap = (TMAP_TYP *)malloc(sizeof(TMAP_TYP) *
				(luvb.tms->mbrmax - luvb.tms->mbrmin + HISTEP));
	if (luvb.tms->lumap == NULL) {
		DMESG(DMCmemory, "malloc() failed in floatMapHisto");
		goto cleanup;
	}
	histamax = HISTI(luva.tms->hbrmax) - HISTI(luva.tms->hbrmin);
	histbmax = HISTI(luvb.tms->hbrmax) - HISTI(luvb.tms->hbrmin);
	n = 0;
	luvb.tms->lumap[0] = x = 0;
	for (y = 0; y < histbmax; luvb.tms->lumap[++y*HISTEP] = x*HISTEP) {
		n += (unsigned long)(luvb.tms->histo[y]*hcvt + .5);
		while (x < histamax && n >= luva.tms->histo[x])
			n -= luva.tms->histo[x++];
	}
						/* fill in gaps */
	for (y = luvb.tms->mbrmax - luvb.tms->mbrmin; y > histbmax*HISTEP; )
		luvb.tms->lumap[y--] = luvb.tms->lumap[histbmax*HISTEP];
	for (y = histbmax*HISTEP; y >= HISTEP; y -= HISTEP)
		for (x = HISTEP; --x > 0; )
			luvb.tms->lumap[y-HISTEP+x] = ( x*luvb.tms->lumap[y] +
				(HISTEP-x)*luvb.tms->lumap[y-HISTEP] ) / HISTEP;
						/* apply mapping(s) */
	if (ib->csp->format == IPFy)		/* luminance only */
		for (y = 0; y < ib->yres; y++) {
			fp = (float *)ProwPtr(ib,y);
			tmCvLums(luvb.ls, fp, ib->xres);
			for (x = 0; x < ib->xres; x++, fp++) {
				if (luvb.ls[x] < MINBRT)
					continue;
				*fp = tmLuminance(luvb.tms->lumap[luvb.ls[x] -
						luvb.tms->mbrmin]+luva.tms->hbrmin);
			}
		}
	else					/* separate luminance/chroma */
		for (y = 0; y < ib->yres; y++) {
			fp = (float *)ProwPtr(ib,y);
			tmCvColors(luvb.tms, luvb.ls, TM_NOCHROM, (COLOR *)fp, ib->xres);
			for (x = 0; x < ib->xres; x++, fp += 3) {
				const double	den = fp[0] + 15.*fp[1] + 3.*fp[2];
				double	lum;
				int	u, v;
				if (luvb.ls[x] < MINBRT)
					continue;
				u = (int)(UVSCALE*4.*fp[0]/den);
				v = (int)(UVSCALE*9.*fp[1]/den);
				if ((u < 0) | (u > 255) | (v < 0) | (v > 255))
					continue;
				u = mapuv[0][u];
				v = mapuv[1][v];
				lum = tmLuminance(luvb.tms->lumap[luvb.ls[x] -
						luvb.tms->mbrmin]+luva.tms->hbrmin);
				fp[0] = 2.25*u/v * lum;
				fp[1] = lum;
				fp[2] = ((3.*UVSCALE - .75*u)/v - 5.) * lum;
			}
		}
						/* return to original CS */
	ok = PconvertColorSpace(ib, orig_csp, 1.f);
cleanup:
	if (!ok) {
		if (luva.tms->lastError != TM_E_OK)
			DMESG(DMCparameter, tmErrorMessage[luva.tms->lastError]);
		else if (luvb.tms->lastError != TM_E_OK)
			DMESG(DMCparameter, tmErrorMessage[luvb.tms->lastError]);
	}
	free(luva.ls);
	free(luvb.ls);
	tmDone(luva.tms);
	tmDone(luvb.tms);
	return ok;
}

/* Modify image ib to match the histogram of image ia */
int
PmatchHisto(ImgStruct *ib, const ImgStruct *ia)
{
	if ((ia == NULL) | (ib == NULL))
		return false;
	if ((ia->img == NULL) | (ib->img == NULL))
		return false;
	if (ib->csp->dtype != IDTfloat &&
			(ia->csp->dtype != ib->csp->dtype) |
			(ImgPixelLen[ia->csp->format] !=
				ImgPixelLen[ib->csp->format])) {
		DMESG(DMCparameter,
			"Cannot match histograms with different pixel types");
		return false;
	}
	DMESGF(DMCtrace, "Matching histogram of %s image",
				PdescribeCS(ia->csp, NULL));
	switch (ib->csp->dtype) {
	case IDTubyte:
		return byteMapHisto(ib, ia);
	case IDTushort:
		return shortMapHisto(ib, ia);
	case IDTfloat:
		return floatMapHisto(ib, ia);
	}
	DMESG(DMCparameter, "Unsupported histogram match type");
	return false;
}
