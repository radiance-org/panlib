/*
 *  pisum.c
 *
 *  Pancing image summing and blending routines.
 *
 *  Created by Greg Ward on 5/29/08.
 *  Copyright 2016 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "system.h"
#include "dmessage.h"
#include "pimage.h"
#include "tiff.h"		/* needed for uint16, etc. */

#ifndef true
#define true	1
#define false	0
#endif

/* Floating point summation callback */
static int
sumFloat(const uby8 *slp, int len, int y, void *udp)
{
	const float *	ps = (const float *)slp;
	ImgStruct *	ib = (ImgStruct *)udp;
	int		n = len * ImgPixelLen[ib->csp->format];
	float *		pd = (float *)ProwPtr(ib, y);

	while (n--)
		*pd++ += *ps++;

	return len;
}

/* Unsigned byte summation callback */
static int
sumUByte(const uby8 *ps, int len, int y, void *udp)
{
	ImgStruct *	ib = (ImgStruct *)udp;
	int		n = len * ImgPixelLen[ib->csp->format];
	uby8 *		pd = ProwPtr(ib, y);
	int		res;

	while (n--) {
		res = *pd + *ps++;
		if (res > 0xff) res = 0xff;
		*pd++ = res;
	}
	return len;
}

/* Unsigned short summation callback */
static int
sumUShort(const uby8 *slp, int len, int y, void *udp)
{
	const uint16 *	ps = (const uint16 *)slp;
	ImgStruct *	ib = (ImgStruct *)udp;
	int		n = len * ImgPixelLen[ib->csp->format];
	uint16 *	pd = (uint16 *)ProwPtr(ib, y);
	int		res;

	while (n--) {
		res = *pd + *ps++;
		if (res > 0xffff) res = 0xffff;
		*pd++ = res;
	}
	return len;
}

/* Sum one image into another of the same size (but possibly different type) */
int
PsumImage(ImgStruct *ib, const ImgStruct *ia, const float sf)
{
	static PscanlineMethod	*spa[] = {
		sumUByte, sumUShort, NULL, sumFloat
	};
						/* argument checks */
	if (ib == NULL || ib->csp == NULL)
		return false;
	if (ia == NULL || (ia->img == NULL) | (ia->csp == NULL))
		return false;
	if (ib->img == NULL)			/* starting fresh */
		return PmapImage(ib, ia, sf);
	if (sf == 0)				/* nothing to add? */
		return true;
	if (ib->csp->format == IPFn ||
			!PmatchColorSpace(ib->csp, &ICS_Y, PICMgamma)) {
		DMESG(DMCparameter, "Unsupported non-linear image summation");
		return false;
	}
	if ((ib->xres != ia->xres) | (ib->yres != ia->yres)) {
		DMESG(DMCparameter, "Mismatched image dimensions");
		return false;
	}
	if (ia->img == ib->img ? (ImgPixelSize(ia->csp) != ImgPixelSize(ib->csp))
				: PimagesOverlap(ib, ia)) {
		DMESG(DMCparameter, "Cannot sum image onto itself");
		return false;
	}
						/* check for direct sum */
	if (ib->csp->dtype == IDTfloat &&
			PmatchColorSpace(ib->csp, ia->csp, PICMall)) {
		int	y;
		for (y = ia->yres; y--; ) {
			const float *	sfp = (const float *)ProwPtr(ia, y);
			float *		dfp = (float *)ProwPtr(ib, y);
			int		n = ia->xres*ImgPixelLen[ia->csp->format];
			
			while (n--)
				*dfp++ += *sfp++ * sf;
		}
		return true;
	}
						/* callback summation */
	return (PconvertImageCB(ia, ib->csp, sf, spa[ib->csp->dtype], ib) >= 0);
}

/* Blend between two images using separate alpha channel */
int
PblendImages(ImgStruct *ib, const ImgStruct *ia0, const ImgStruct *ialpha,
					const ImgStruct *ia1)
{
	float *			bbuf = NULL;
	void *			bccv = NULL;
	float *			a0buf = NULL;
	void *			a0ccv = NULL;
	float *			a1buf = NULL;
	void *			a1ccv = NULL;
	float *			alphabuf = NULL;
	void *			alphaccv = NULL;
	const ImgColorSpace *	wcsp = NULL;
	int			nc, y;
						/* argument checks */
	if (ib == NULL || ib->csp == NULL)
		return false;
	if (ia0 == NULL || (ia0->img == NULL) | (ia0->csp == NULL))
		return false;
	if (ia1 == NULL || (ia1->img == NULL) | (ia1->csp == NULL))
		return false;
	if (ialpha == NULL || (ialpha->img == NULL) | (ialpha->csp == NULL))
		return false;
	if ((ia0->xres != ialpha->xres) | (ia0->yres != ialpha->yres) ||
			(ia1->xres != ialpha->xres) | (ia1->yres != ialpha->yres)) {
		DMESG(DMCparameter, "Mismatched source dimensions");
		return false;
	}
	if (ib->img == NULL) {			/* starting fresh? */
		ib->xres = ialpha->xres; ib->yres = ialpha->yres;
		if (!PnewImage(ib, .0))
			return false;
	} else if ((ib->xres != ialpha->xres) | (ib->yres != ialpha->yres)) {
		DMESG(DMCparameter, "Mismatched destination dimensions");
		return false;
	}
						/* set up conversion buffers */
	nc = ImgPixelLen[ib->csp->format];
	if (!PmatchColorSpace(ib->csp, &ICS_Y, PICMdtype|PICMgamma)) {
		if (ia0->csp->dtype == IDTfloat && ImgPixelLen[ia0->csp->format] == nc)
			wcsp = ia0->csp;
		else if (ia1->csp->dtype == IDTfloat && ImgPixelLen[ia1->csp->format] == nc)
			wcsp = ia1->csp;
		else
			wcsp = (nc==1) ? &ICS_Y : &ICS_XYZ;
		bccv = PcreateColorConv(ib->csp, wcsp, 1.f, 1);
		if (bccv == NULL) return false;
		bbuf = (float *)malloc(ImgPixelSize(wcsp)*ib->xres);
		if (bbuf == NULL) goto memerr;
	} else
		wcsp = ib->csp;
	if (!PmatchColorSpace(ia0->csp, wcsp, PICMall)) {
		a0ccv = PcreateColorConv(wcsp, ia0->csp, 1.f, 1);
		if (a0ccv == NULL) return false;
		a0buf = (float *)malloc(ImgPixelSize(wcsp)*ia0->xres);
		if (a0buf == NULL) goto memerr;
	}
	if (!PmatchColorSpace(ia1->csp, wcsp, PICMall)) {
		a1ccv = PcreateColorConv(wcsp, ia1->csp, 1.f, 1);
		if (a1ccv == NULL) return false;
		a1buf = (float *)malloc(ImgPixelSize(wcsp)*ia1->xres);
		if (a1buf == NULL) goto memerr;
	}
	if (!PmatchColorSpace(ialpha->csp, &ICS_Alpha, PICMall) &&
			!PmatchColorSpace(ialpha->csp, &ICS_Y, PICMall)) {
		alphaccv = PcreateColorConv(&ICS_Alpha, ialpha->csp, 1.f, 1);
		if (alphaccv == NULL) return false;
		alphabuf = (float *)malloc(sizeof(float)*ialpha->xres);
		if (alphabuf == NULL) goto memerr;
		if (ImgPixelLen[ialpha->csp->format] != 1)
			DMESG(DMCwarning, "Converting alpha to single channel");
	}
	strcpy(dmessage_buf, "Blending ");
	PdescribeCS(ia0->csp, strchr(dmessage_buf,'\0'));
	if (PmatchColorSpace(ia0->csp, ia1->csp, PICMptype|PICMgamma)) {
		strcat(dmessage_buf, " images");
	} else {
		strcat(dmessage_buf, " and ");
		PdescribeCS(ia1->csp, strchr(dmessage_buf,'\0'));
	}
	if (!PmatchColorSpace(ia0->csp, ib->csp, PICMptype|PICMgamma)) {
		strcat(dmessage_buf, " to ");
		PdescribeCS(ib->csp, strchr(dmessage_buf,'\0'));
	}
	DMESG(DMCtrace, dmessage_buf);
						/* process each scanline */
	for (y = ib->yres; y--; ) {
		int		n = ib->xres;
		const float *	s0fp;
		const float *	s1fp;
		const float *	alfp;
		float *		dfp;

		if (!a0ccv)
			s0fp = (const float *)ProwPtr(ia0,y);
		else if (PmapPixels((uby8 *)a0buf, ProwPtr(ia0,y), ia0->xres, a0ccv))
			s0fp = a0buf;
		else
			break;
		if (!a1ccv)
			s1fp = (const float *)ProwPtr(ia1,y);
		else if (PmapPixels((uby8 *)a1buf, ProwPtr(ia1,y), ia1->xres, a1ccv))
			s1fp = a1buf;
		else
			break;
		if (!alphaccv)
			alfp = (const float *)ProwPtr(ialpha,y);
		else if (PmapPixels((uby8 *)alphabuf, ProwPtr(ialpha,y), ialpha->xres, alphaccv))
			alfp = alphabuf;
		else
			break;
		if (bccv)
			dfp = bbuf;
		else
			dfp = (float *)ProwPtr(ib,y);

		while (n--) {			/* process float pixels */
			int	c = nc;
			while (c--)
				*dfp++ = (1.f - *alfp) * *s0fp++
					+ *alfp * *s1fp++ ;
			++alfp;
		}
		if (bccv && !PmapPixels(ProwPtr(ib,y), (uby8 *)bbuf, ib->xres, bccv))
			break;
	}
						/* clean up */
	if (bccv) {
		PfreeColorConv(bccv);
		free(bbuf);
	}
	if (a0ccv) {
		PfreeColorConv(a0ccv);
		free(a0buf);
	}
	if (a1ccv) {
		PfreeColorConv(a1ccv);
		free(a1buf);
	}
	if (alphaccv) {
		PfreeColorConv(alphaccv);
		free(alphabuf);
	}
	if (y >= 0) {				/* did we abort? */
		PfreeImage(ib);
		return false;
	}
	return true;				/* else we're good */
memerr:
	DMESG(DMCmemory, "malloc() failed in PblendImages()");
	PfreeImage(ib);
	return false;
}

/* Blend image operation given by callback */
int
PblendImageOpCB(ImgStruct *ib, const ImgStruct *ialpha, PimageOp *op, void *udp)
{
	const int	fmarg = 4;		/* filter margin */
	int		psiz;
	int		ok;
	int		x, y, i;
	ImgStruct	iwork, isub, asub;
	ImgRect		subrect;
	double		frac;
	const uby8	*p;

	if (ib == NULL || (ib->csp == NULL) | (ib->img == NULL))
		return false;
	if (ialpha == NULL || (ialpha->img == NULL) | (ialpha->csp == NULL))
		return false;
	if ((ib->xres != ialpha->xres) | (ib->yres != ialpha->yres)) {
		DMESG(DMCparameter, "Mismatched image dimensions");
		return false;
	}
	if (op == NULL)
		return false;			/* assume this is a mistake */
	psiz = ImgPixelSize(ialpha->csp);
	/*
	 * Optimize operation by trimming zero borders in alpha mask.
	 * This way, a small alpha mask somewhere in the image only
	 * operates on that region.  Presumably, the operation and
	 * subsequent blend are more expensive than finding margins.
	 */
	subrect.xright = subrect.ybottom = 0;
	subrect.xleft = ialpha->xres;
	subrect.ytop = ialpha->yres;
	for (y = 0; y < ialpha->yres; y++) {	/* search from top-left */
		p = ProwPtr(ialpha,y);
		for (x = 0; x < ialpha->xres; x++) {
			for (i = psiz; i--; )
				if (*p++) break;
			if (i >= 0) {
				subrect.xleft = x;
				subrect.ytop = y;
				break;
			}
		}
		if (x < ialpha->xres) break;
	}
	if (y >= ialpha->yres) {		/* alpha zero everywhere? */
		DMESG(DMCtrace, "Zero alpha everywhere => no-op");
		return true;
	}
	for (y = ialpha->yres; --y > 0; ) {	/* search from bottom-right */
		p = ProwPtr(ialpha,y) + ialpha->xres*psiz;
		for (x = ialpha->xres; x-- > 0; ) {
			for (i = psiz; i--; )
				if (*--p) break;
			if (i >= 0) {
				subrect.xright = x+1;
				subrect.ybottom = y+1;
				break;
			}
		}
		if (x >= 0) break;
	}
	for (x = 0; x < subrect.xleft; x++) {	/* search from left */
		for (y = subrect.ytop; y < subrect.ybottom; y++) {
			p = PpixP(ialpha,x,y,psiz);
			for (i = psiz; i--; )
				if (*p++) break;
			if (i >= 0) {
				subrect.xleft = x;
				break;
			}
		}
	}
	for (x = ialpha->xres; --x > subrect.xright; ) {	/* from right */
		for (y = subrect.ytop; y < subrect.ybottom; y++) {
			p = PpixP(ialpha,x,y,psiz);
			for (i = psiz; i--; )
				if (*p++) break;
			if (i >= 0) {
				subrect.xright = x+1;
				break;
			}
		}
	}
	subrect.xleft -= fmarg;			/* add margins for filtering */
	subrect.xleft *= (subrect.xleft > 0);
	if ((subrect.xright += fmarg) > ib->xres)
		subrect.xright = ib->xres;
	subrect.ytop -= fmarg;
	subrect.ytop *= (subrect.ytop > 0);
	if ((subrect.ybottom += fmarg) > ib->yres)
		subrect.ybottom = ib->yres;
	if (!PlinkSubimage(&isub, ib, &subrect))
		return false;			/* should not happen here! */
	PlinkSubimage(&asub, ialpha, &subrect);
	frac = PrectArea(&subrect) / (double)(ib->xres*ib->yres);
	DTESTF((frac < .999), DMCtrace, "Operating on %.1f%% submimage", 100.*frac);
	iwork.img = NULL;			/* clear destination */
	iwork.csp = isub.csp;			/* call-back can reassign CS */
	iwork.xres = isub.xres; iwork.yres = isub.yres;
	ok = (*op)(&iwork, &isub, udp);		/* operate and blend result */
	if (ok > 0)
		ok = PblendImages(&isub, &isub, &asub, &iwork);
	PfreeImage(&iwork);			/* free temps */
	PfreeImage(&isub);
	PfreeImage(&asub);
	return ok;
}
