/*
 *  pblendpan2.c
 *
 *  Created by Greg Ward on 4/26/06.
 *  Copyright 2006 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include "dmessage.h"
#include "pimage.h"

#ifndef true
#define true	1
#define false	0
#endif

#undef MIN
#define MIN(a,b)	((a)<(b) ? (a) : (b))

#ifndef DECOMP_DIV
#define DECOMP_DIV	4	/* resampling at each decomposition layer */
#endif

#ifndef MAX_BLEND
#define MAX_BLEND	400	/* maximum blend width */
#endif

extern int
PblendPan2(ImgStruct *blnd, const ImgRect *anch,
			const ImgStruct *inpa, const ImgStruct *inpb);

/*
 * This is essentially a reimplementation of Burt & Adelson's 1983
 * TOG paper, "A Multiresolution Spline with Application to Image
 * Mosaics."  The main difference is that I use a multiplicative
 * separation for better behavior with HDR sources.
 */

/* Separate image into low and high frequency versions */
static int
sepFreq2(ImgStruct *hfi, ImgStruct *smi, const ImgStruct *src)
{
	ImgStruct	tmpImg;
	int		x, y;
						/* copy/map source */
	hfi->xres = tmpImg.xres = src->xres;
	hfi->yres = tmpImg.yres = src->yres;
	if (!PmapImage(hfi, src, 1.f))
		return false;
						/* resample to get LF */
	tmpImg.csp = hfi->csp;
	tmpImg.img = NULL;
	if (!PsizeImage(smi, hfi, PSgaussian) ||
			!PsizeImage(&tmpImg, smi, PSlinear)) {
		PfreeImage(smi);
		return false;
	}
						/* divide to get HF */
	DASSERT(hfi->csp->dtype == IDTfloat);
	for (y = hfi->yres; y-- > 0; ) {
		const float *	sp = (const float *)ProwPtr(&tmpImg,y);
		float *		dp = (float *)ProwPtr(hfi,y);
		for (x = ImgPixelLen[hfi->csp->format]*hfi->xres; x-- > 0; ) {
			if (*sp == 0) { ++dp; ++sp; continue; }
			*dp++ /= *sp++;
		}
	}
	PfreeImage(&tmpImg);
	return true;
}

/* Splice two images at midline of indicated region */
static int
spliceImg(ImgStruct *retImg, const ImgStruct *luImg, const ImgStruct *rbImg,
			const ImgRect *blndRgn)
{
	ImgRect		rbRgn;
	ImgStruct	rbSrc;

	if (!PmapImage(retImg, luImg, 1.f))
		return false;
	if (!PlegalRect(blndRgn, retImg->xres, retImg->yres)) {
		DMESG(DMCparameter, "Illegal rectangle in spliceImg");
		return false;
	}
	rbRgn.xleft = rbRgn.ytop = 0;
	rbRgn.xright = retImg->xres;
	rbRgn.ybottom = retImg->yres;
	if ((blndRgn->xleft > 0) | (blndRgn->xright < retImg->xres))
		rbRgn.xleft = (blndRgn->xleft + blndRgn->xright) / 2;
	else
		rbRgn.ytop = (blndRgn->ytop + blndRgn->ybottom) / 2;
	if (!PlinkSubimage(&rbSrc, rbImg, &rbRgn)) {
		DMESG(DMCparameter, "PlinkSubimage failed in spliceImg");
		return false;
	}
	if (!PcopyImage(retImg, &rbSrc, rbRgn.xleft, rbRgn.ytop))
		return false;
	PfreeImage(&rbSrc);
	return true;
}

/* Multiply src floating-point image into dst */
static void
multImg(ImgStruct *dst, const ImgStruct *src)
{
	int	x, y;

	DASSERT(src->csp->dtype == IDTfloat);
	DASSERT(dst->csp->dtype == IDTfloat);
	DASSERT((dst->xres == src->xres) & (dst->yres == src->yres));
	DASSERT(ImgPixelLen[dst->csp->format] == ImgPixelLen[src->csp->format]);
	for (y = 0; y < dst->yres; y++) {
		const float *	sp = (const float *)ProwPtr(src,y);
		float *		dp = (float *)ProwPtr(dst,y);
		for (x = ImgPixelLen[dst->csp->format]*dst->xres; x--; )
			*dp++ *= *sp++;
	}
}

/* Blend two image regions using a frequency decomposition */
static int
joinDecomp(ImgStruct *retImg, const ImgStruct *luImg, const ImgStruct *rbImg,
			const ImgRect *blndRgn)
{
	ImgStruct	smlLU, smlRB, smlRet;
	ImgStruct	hfLU, hfRB;
	ImgRect		smlRgn;
	ImgStruct	rsmRet;

	if (retImg == NULL || retImg->csp == NULL)
		return false;
					/* check recursion limit */
	if ((blndRgn->xright - blndRgn->xleft <= 2*DECOMP_DIV-1) |
			(blndRgn->ybottom - blndRgn->ytop <= 2*DECOMP_DIV-1))
		return spliceImg(retImg, luImg, rbImg, blndRgn);
		
	if (retImg->img == NULL) {
		retImg->xres = luImg->xres;
		retImg->yres = luImg->yres;
	}				/* else should verify size? */
	smlRgn.xleft = blndRgn->xleft / DECOMP_DIV;
	smlRgn.xright = blndRgn->xright / DECOMP_DIV;
	smlRgn.ytop = blndRgn->ytop / DECOMP_DIV;
	smlRgn.ybottom = blndRgn->ybottom / DECOMP_DIV;
	smlLU.csp = hfLU.csp = luImg->csp;
	smlLU.xres = luImg->xres / DECOMP_DIV;
	smlLU.yres = luImg->yres / DECOMP_DIV;
	smlLU.img = hfLU.img = NULL;
	if (!sepFreq2(&hfLU, &smlLU, luImg))
		return false;
	smlRB.csp = hfRB.csp = rbImg->csp;
	smlRB.xres = rbImg->xres / DECOMP_DIV;
	smlRB.yres = rbImg->yres / DECOMP_DIV;
	smlRB.img = hfRB.img = NULL;
	if (!sepFreq2(&hfRB, &smlRB, rbImg))
		return false;
	smlRet.csp = retImg->csp;
	smlRet.img = NULL;
	if (!joinDecomp(&smlRet, &smlLU, &smlRB, &smlRgn))
		return false;
	PfreeImage(&smlLU); PfreeImage(&smlRB);
	rsmRet.csp = retImg->csp;
	rsmRet.xres = retImg->xres;
	rsmRet.yres = retImg->yres;
	rsmRet.img = NULL;
	if (!PsizeImage(&rsmRet, &smlRet, PSlinear))
		return false;
	PfreeImage(&smlRet);
	if (!spliceImg(retImg, &hfLU, &hfRB, blndRgn))
		return false;
	PfreeImage(&hfLU); PfreeImage(&hfRB);
	multImg(retImg, &rsmRet);
	PfreeImage(&rsmRet);
	return true;
}

/* Blend sections of a panorama:
 *	Input is two overlapping images in matching color spaces.
 *	The rectangle "anch" has its left-top coordinate
 *	set to a feature in "inpa" that matches "inpb" when the
 *	right-bottom corner of "inpb" is at the right-bottom corner of
 *	"anch".  (This is not necessarily a legal rectangle.)
 *	On completion, PblendPan2() creates the image "blnd" to cover
 *	the overlapping region between "inpa" and "inpb".
 */
int
PblendPan2(ImgStruct *blnd, const ImgRect *anch,
			const ImgStruct *inpa, const ImgStruct *inpb)
{
	ImgStruct	srca, srcb;
	ImgStruct	subImg, iblnd;
	ImgRect		rgn;
	int		xmatch, ymatch;
						/* check arguments */
	if (blnd == NULL || blnd->csp == NULL)
		return false;
	if (anch == NULL)
		return false;
	if (inpa == NULL || (inpa->img == NULL) | (inpa->csp == NULL))
		return false;
	if (inpb == NULL || (inpb->img == NULL) | (inpb->csp == NULL))
		return false;
						/* create source images */
	if (blnd->csp->dtype == IDTfloat)
		iblnd.csp = srca.csp = srcb.csp = blnd->csp;
	else if (blnd->csp->format == IPFy)
		iblnd.csp = srca.csp = srcb.csp = &ICS_Y;
	else
		iblnd.csp = srca.csp = srcb.csp = &ICS_RGB709;
	iblnd.img = srca.img = srcb.img = subImg.img = NULL;
	rgn.xleft = anch->xright - inpb->xres;
	if (rgn.xleft < 0) rgn.xleft = 0;
	rgn.xright = anch->xright;
	if (rgn.xright > inpa->xres) rgn.xright = inpa->xres;
	rgn.ytop = anch->ybottom - inpb->yres;
	if (rgn.ytop < 0) rgn.ytop = 0;
	rgn.ybottom = anch->ybottom;
	if (rgn.ybottom > inpa->yres) rgn.ybottom = inpa->yres;
	xmatch = anch->xleft - rgn.xleft;
	ymatch = anch->ytop - rgn.ytop;
	if (!PlinkSubimage(&subImg, inpa, &rgn) ||
			!PmapImage(&srca, &subImg, 1.f))
		goto fail;
	PfreeImage(&subImg);
	rgn.xleft += inpb->xres - anch->xright;
	rgn.xright += inpb->xres - anch->xright;
	rgn.ytop += inpb->yres - anch->ybottom;
	rgn.ybottom += inpb->yres - anch->ybottom;
	if (!PlinkSubimage(&subImg, inpb, &rgn) ||
			!PmapImage(&srcb, &subImg, 1.f))
		goto fail;
	PfreeImage(&subImg);
						/* blend scanlines */
	if (abs(inpb->xres - anch->xright)*(inpa->yres + inpb->yres) >
			abs(inpb->yres - anch->ybottom)*(inpa->xres + inpb->xres)) {
		int	width = MIN(xmatch,srca.xres-xmatch);
		if (width > MAX_BLEND) width = MAX_BLEND;
		rgn.ytop = 0; rgn.ybottom = srca.yres;
		rgn.xleft = xmatch - width/2;
		rgn.xright = xmatch + width/2;
		if (inpb->xres < anch->xright) {
			if (!joinDecomp(&iblnd, &srca, &srcb, &rgn))
				goto fail;
		} else {
			if (!joinDecomp(&iblnd, &srcb, &srca, &rgn))
				goto fail;
		}
	} else {				/* vertical panorama */
		int	height = 2*MIN(ymatch,srca.yres-ymatch);
		if (height > MAX_BLEND) height = MAX_BLEND;
		rgn.xleft = 0; rgn.xright = srca.xres;
		rgn.ytop = ymatch - height/2;
		rgn.ybottom = ymatch + height/2;
		if (inpb->yres < anch->ybottom) {
			if (!joinDecomp(&iblnd, &srca, &srcb, &rgn))
				goto fail;
		} else {
			if (!joinDecomp(&iblnd, &srcb, &srca, &rgn))
				goto fail;
		}
	}
	PfreeImage(&srca); PfreeImage(&srcb);
/*
float *matchP = (float *)PpixPtr(&iblnd,xmatch,ymatch);
matchP[0] = 0;
matchP[1] = 1000.f;
matchP[2] = 0;
*/						/* copy/convert result */
	if (!PmapImage(blnd, &iblnd, 1.f))
		goto fail;
	PfreeImage(&iblnd);
	return true;
fail:
	PfreeImage(&subImg); PfreeImage(&iblnd);
	PfreeImage(&srca); PfreeImage(&srcb);
	PfreeImage(blnd);
	return false;
}
