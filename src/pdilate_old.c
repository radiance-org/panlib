/*
 *  pdilate.c
 *  panlib
 *
 *  Dilate or erode an image.
 *  Algorithm is roughly O(N*M) where N is #pixels and M is disk radius.
 *  Worst-case performance is O(N*M^2), which can happen in vertical gradients.
 *
 *  Created by Greg Ward on 1/9/16.
 *  Copyright 2016 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"
#include "dmessage.h"
#include "pimage.h"

#ifndef true
#define true	1
#define false	0
#endif

/* Structure to hold gray scanline group and index extreme values */
typedef struct {
	float		r;			/* radius (negative to erode) */
	int		ir;			/* positive radius (pixels) */
	int		psiz;			/* pixel size */
	int		cy;			/* current (central) y-position */
	int		(*cmpf)(const void *p1, const void *p2);
	const ImgStruct	*inp;			/* input image poiinter */
	void		*ccp;			/* color conversion pointer */
	short		(*srcpos)[2];		/* extreme source positions */
	short		(*rcrescent)[2];	/* right-crescent source pos. */
	int		nrcr;			/* #points in right crescent */
	uby8		*sca[1];		/* scanlines (extends struct) */
} ScanBar;

#define SBsrcP(sb,x)		((sb)->sca[(sb)->srcpos[x][1]-(sb)->cy+(sb)->ir] + \
						(sb)->srcpos[x][0]*(sb)->psiz)

#define firstBetter(p1,p2,sb)	((1 - 2*((sb)->r < 0)) * (*(sb)->cmpf)( \
					(sb)->sca[(p1)[1]-(sb)->cy+(sb)->ir] + \
						(p1)[0]*(sb)->psiz, \
					(sb)->sca[(p2)[1]-(sb)->cy+(sb)->ir] + \
						(p2)[0]*(sb)->psiz ) > 0)

/* Compare unsigned byte values */
static int
cmp_byte(const void *p1, const void *p2)
{
	return (int)(*(const uby8 *)p1) - (int)(*(const uby8 *)p2);
}

/* Compare unsigned short values */
static int
cmp_ushort(const void *p1, const void *p2)
{
	return (int)(*(const unsigned short *)p1) - (int)(*(const unsigned short *)p2);
}

/* Compare float values */
static int
cmp_float(const void *p1, const void *p2)
{
	float	v1 = *(const float *)p1;
	float	diff = v1 - *(const float *)p2;

	if (v1 > 1e-10)
		v1 *= 1e-5;
	else if (v1 < -1e-10)
		v1 *= -1e-5;
	else
		v1 = 1e-15;
	if (diff > v1) return 1;
	if (diff < -v1) return -1;
	return 0;
}

/* Free allocated scan bar */
static void
freeScanBar(ScanBar *sbp)
{
	if (sbp->ccp != NULL) {			/* free conversion buffers */
		int	i = 2*sbp->ir+1;
		while (i--)
			free(sbp->sca[i]);
		PfreeColorConv(sbp->ccp);
	}
	free(sbp->rcrescent);
	free(sbp->srcpos);
	free(sbp);
}

/* Allocate new structure for scan bar conversion */
static ScanBar *
newScanBar(const ImgColorSpace *dcsp, const ImgStruct *simp, float r)
{
	int		i = (int)(fabs(r) + .9999);
	ScanBar		*sbp;

	if ((i >= simp->xres) | (i >= simp->yres)) {
		DMESG(DMCparameter, "Radius too large for image");
		return NULL;
	}
	if (dcsp->format != IPFy) {
		DMESG(DMCparameter, "Destination image must be grayscale");
		return NULL;
	}
	sbp = (ScanBar *)malloc(sizeof(ScanBar) + 2*sizeof(uby8 *)*i);
	if (sbp == NULL) {
		DMESG(DMCmemory, "Cannot allocate scan bar");
		return NULL;
	}
	switch (dcsp->dtype) {
	case IDTubyte:
		sbp->cmpf = cmp_byte;
		break;
	case IDTushort:
		sbp->cmpf = cmp_ushort;
		break;
	case IDTfloat:
		sbp->cmpf = cmp_float;
		break;
	default:
		DMESG(DMCparameter, "Unsupported destination pixel type");
		free(sbp);
		return NULL;
	}
	sbp->r = r;
	sbp->ir = i;
	sbp->psiz = ImgPixelSize(dcsp);
	sbp->cy = -1;
	sbp->inp = simp;
	if (PmatchColorSpace(dcsp, simp->csp, PICMall)) {
		sbp->ccp = NULL;
	} else if ((sbp->ccp = PcreateColorConv(dcsp, simp->csp, 1.f, 1)) == NULL) {
		free(sbp);
		return NULL;
	}
	sbp->rcrescent = NULL;
	memset(sbp->sca, 0, (2*sbp->ir+1)*sizeof(uby8 *));
	sbp->srcpos = (short (*)[2])malloc(2*sizeof(short)*simp->xres);
	if (sbp->srcpos == NULL)
		goto memerr;
	sbp->rcrescent = (short (*)[2])malloc(2*sizeof(short)*(int)(3.142*sbp->ir+1.5));
	if (sbp->rcrescent == NULL)
		goto memerr;
	sbp->nrcr = 0;				/* create right crescent */
	for (i = 0; i*i <= r*r; i++) {
		int	x = (int)sqrt(r*r - i*i);
		sbp->rcrescent[sbp->nrcr][0] = x;
		sbp->rcrescent[sbp->nrcr++][1] = i;
		if (!i) continue;
		sbp->rcrescent[sbp->nrcr][0] = x;
		sbp->rcrescent[sbp->nrcr++][1] = -i;
	}
	for (i = 2*sbp->ir+1; i--; )		/* initialize scanlines */
		if (sbp->ccp != NULL) {
			sbp->sca[i] = (uby8 *)malloc(sbp->psiz*simp->xres);
			if (sbp->sca[i] == NULL)
				goto memerr;
			if (i <= sbp->ir)
				continue;
			if (!PmapPixels(sbp->sca[i], ProwPtr(simp,i-sbp->ir-1),
					simp->xres, sbp->ccp)) {
				freeScanBar(sbp);
				return NULL;
			}
		} else if (i > sbp->ir) {
			sbp->sca[i] = ProwPtr(simp,i-sbp->ir-1);
		}
	return sbp;
memerr:
	DMESG(DMCmemory, "malloc() failed in newScanBar()");
	freeScanBar(sbp);
	return NULL;
}

/* Dilate (or erode) the current scanline */
static int
dilateScan(uby8 *dsp, ScanBar *sbp, const int y)
{
	const int	xs = 1 - 2*(y&1);	/* serpentine for gradients */
	uby8		*slp;
	int		n, i;
	short		tpt[2];
						/* advance & check */
	if ((y != ++sbp->cy) | (y >= sbp->inp->yres)) {
		DMESG(DMCparameter, "Bad y-position in call to dilateScan()");
		return false;
	}
	slp = sbp->sca[0];			/* rotate scanlines */
	memmove(sbp->sca, sbp->sca+1, sizeof(uby8 *)*2*sbp->ir);
	sbp->sca[2*sbp->ir] = slp;
	if (y+sbp->ir < sbp->inp->yres) {
		if (sbp->ccp == NULL)
			sbp->sca[2*sbp->ir] = ProwPtr(sbp->inp,y+sbp->ir);
		else if (!PmapPixels(slp, ProwPtr(sbp->inp,y+sbp->ir),
						sbp->inp->xres, sbp->ccp))
			return false;
	}
						/* identify min/max positions */
	for (n = sbp->inp->xres; n--; ) {
		const int	x = (xs<0) ? n : sbp->inp->xres-1 - n;
		int		axisFlags = 0;	/* try neighbors' extrema */
		int		dx, dy;
		if ((0 <= x-xs) & (x-xs < sbp->inp->xres)) {
			dx = x - sbp->srcpos[x-xs][0];
			dy = y - sbp->srcpos[x-xs][1];
			axisFlags |= (dx*dx + dy*dy <= sbp->r*sbp->r);
		}
		if (y > 0) {
			dx = x - sbp->srcpos[x][0];
			dy = y - sbp->srcpos[x][1];
			axisFlags |= (dx*dx + dy*dy <= sbp->r*sbp->r) << 1;;
		}
	switch_again:				/* which neighbor works? */
		switch (axisFlags) {
		case 3:				/* both are usable */
			i = (*sbp->cmpf)(SBsrcP(sbp,x-xs), SBsrcP(sbp,x));
			if (i) {		/* is one of them better? */
				axisFlags = 1 << ((i > 0) ^ (sbp->r > 0));
				goto switch_again;
			}
						/* prefer source further down */
			if (sbp->srcpos[x-xs][1] > sbp->srcpos[x][1]) {
				sbp->srcpos[x][0] = sbp->srcpos[x-xs][0];
				sbp->srcpos[x][1] = sbp->srcpos[x-xs][1];
			}
						/* check quarter-circle */
			for (i = sbp->nrcr; i--; ) {
				if (sbp->rcrescent[i][1] < 0)
					continue;
				tpt[0] = x + xs*sbp->rcrescent[i][0];
				if ((0 > tpt[0]) | (tpt[0] >= sbp->inp->xres))
					continue;
				tpt[1] = y + sbp->rcrescent[i][1];
				if (tpt[1] >= sbp->inp->yres)
					continue;
				if (!firstBetter(sbp->srcpos[x], tpt, sbp)) {
					sbp->srcpos[x][0] = tpt[0];
					sbp->srcpos[x][1] = tpt[1];
				}
			}
			break;
		case 1:				/* horiz. neighbor usable */
			sbp->srcpos[x][0] = sbp->srcpos[x-xs][0];
			sbp->srcpos[x][1] = sbp->srcpos[x-xs][1];
						/* check horiz. half-circle */
			for (i = sbp->nrcr; i--; ) {
				tpt[0] = x + xs*sbp->rcrescent[i][0];
				if ((0 > tpt[0]) | (tpt[0] >= sbp->inp->xres))
					continue;
				tpt[1] = y + sbp->rcrescent[i][1];
				if ((0 > tpt[1]) | (tpt[1] >= sbp->inp->yres))
					continue;
				if (!firstBetter(sbp->srcpos[x], tpt, sbp)) {
					sbp->srcpos[x][0] = tpt[0];
					sbp->srcpos[x][1] = tpt[1];
				}
			}
			break;
		case 2:				/* above neighbor usable */
						/* check lower half-circle */
			for (i = sbp->nrcr; i--; ) {
				tpt[0] = x + sbp->rcrescent[i][1];
				if ((0 > tpt[0]) | (tpt[0] >= sbp->inp->xres))
					continue;
				tpt[1] = y + sbp->rcrescent[i][0];
				if (tpt[1] >= sbp->inp->yres)
					continue;
				if (!firstBetter(sbp->srcpos[x], tpt, sbp)) {
					sbp->srcpos[x][0] = tpt[0];
					sbp->srcpos[x][1] = tpt[1];
				}
			}
			break;
		case 0:				/* perform fresh search */
			sbp->srcpos[x][0] = x;
			sbp->srcpos[x][1] = y;
			for (dy = -sbp->ir; dy <= sbp->ir; dy++) {
				int	w;
				if (y+dy < 0) continue;
				if (y+dy >= sbp->inp->yres) break;
				if (dy*dy > sbp->r*sbp->r)
					continue;
				w = (int)sqrt(sbp->r*sbp->r - dy*dy);
				for (dx = -w; dx <= w; dx++) {
					if (x+xs*dx < 0) continue;
					if (x+xs*dx >= sbp->inp->xres) break;
					if (!(dx|dy))
						continue;
					tpt[0] = x + xs*dx;
					tpt[1] = y + dy;
					if (!firstBetter(sbp->srcpos[x], tpt, sbp)) {
						sbp->srcpos[x][0] = tpt[0];
						sbp->srcpos[x][1] = tpt[1];
					}
				}
			}
			break;
		}
						/* transfer selected pixel */
		memcpy(dsp + x*sbp->psiz, SBsrcP(sbp,x), sbp->psiz);
	}
	return true;
}

/* Dilate (or erode) image using a disk of the given radius */
int
PdilateImage(ImgStruct *ib, const ImgStruct *ia, float rad)
{
	ScanBar		*sbp;
	int		y;
						/* check parameters */
	if ((ib == NULL) | (ia == NULL) ||
			(ia->csp == NULL) | (ia->img == NULL))
		return false;
	if (ib->csp == NULL) {
		ib->csp = ia->csp;
		ib->img = NULL;
	}
	if (ib->img == NULL) {
		ib->xres = ia->xres; ib->yres = ia->yres;
	} else if ((ib->xres != ia->xres) | (ib->yres != ia->yres)) {
		DMESG(DMCparameter, "Mismatched image dimensions");
		return false;
	}
	if (rad*rad < 1)
		return PmapImage(ib, ia, 1.f);
						/* color/buffered output? */
	if (ib->csp->format != IPFy || PimagesOverlap(ib, ia)) {
		ImgStruct	oimg;
		int		c;
		oimg = *ia;
		oimg.img = NULL;
		if (ia->csp->format != IPFy) {	/* multi-channel case */
			ImgColorSpace	grayCS;
			ImgStruct	ogry;
			PcopyCS(&grayCS, ia->csp);
			grayCS.format = IPFy;
			ogry.img = NULL;
			ogry.csp = &grayCS;
			ogry.xres = ia->xres; ogry.yres = ia->yres;
			for (c = ImgPixelLen[ia->csp->format]; c--; ) {
				ImgStruct	igry = ogry;
				igry.img = NULL;
				if (!PcopyComponent(&igry, 0, ia, c) ||
						!PdilateImage(&ogry, &igry, rad)) {
					PfreeImage(&igry);
					PfreeImage(&oimg);
					PfreeImage(ib);
					return false;
				}
				PfreeImage(&igry);
				PcopyComponent(&oimg, c, &ogry, 0);
			}
			PfreeImage(&ogry);
		} else if (!PdilateImage(&oimg, ia, rad)) {
			PfreeImage(ib);
			return false;
		}
		if (ib->img == NULL && PmatchColorSpace(ib->csp, oimg.csp, PICMall))
			c = PlinkImage(ib, &oimg);
		else
			c = PmapImage(ib, &oimg, 1.f);
		PfreeImage(&oimg);
		return c;
	}
						/* grayscale destination */
	sbp = newScanBar(ib->csp, ia, rad);
	if (sbp == NULL)
		return false;
	if (!PnewImage(ib, .0))
		return false;
	for (y = 0; y < ia->yres; y++)
		if (!dilateScan(ProwPtr(ib,y), sbp, y)) {
			freeScanBar(sbp);
			PfreeImage(ib);
			return false;
		}
	freeScanBar(sbp);
	return true;
}
