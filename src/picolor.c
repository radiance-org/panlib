/*
 *  picolor.c
 *  panlib
 *
 *  Pancine color conversion routines.
 *
 *  Created by Greg Ward on 2/16/06.
 *  Copyright 2006 Anyhere Software. All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "system.h"
#include "dmessage.h"
#include "pimage.h"
#include "rtmath.h"		/* needed for types & normal conversion */
#include "color.h"
#include "random.h"

#ifndef M_LN10
#define M_LN10		2.30258509299404568402
#endif

#define gamLINEAR(g)	((0.98f <= (g)) & ((g) <= 1.02f))

#define mkPOS(x)	((x) < 1e-20f ? 1e-20f : (x))

#ifndef true
#define true	1
#define false	0
#endif

#if defined(_WIN32) || defined(_WIN64)
#undef powf
#define powf(x,e)	(float)pow((double)(x), (double)(e))
#undef log10f
#define log10f(x)	(float)log10((double)(x))
#undef expf
#define expf(x)		(float)exp((double)(x))
#endif

#define CMAGIC		0x51F8571E
#define	GTSIZ		1024	/* size of reverse gamma lookup table */
#define GTPART		0.012f	/* lower value for gamma lookup */

int	PrandCQuant = 0;	/* random color quantization? */

#define applyGamma8(v,c)	((v) <= 0 ? 0 : (v) >= 1.f ? 255 : \
					PrandCQuant ? (int)( 255.*powf(v,(c)->outgam) \
								+ frandom() ) : \
					((v) < GTPART) & ((c)->outgam < .9f) ? \
					(int)(256.*powf(v,(c)->outgam)) : \
					(c)->toByte[(int)((v)*(float)GTSIZ)])

/* Conversion structure for color transformation */
typedef struct _ColorCv {
	uint32		magic;		/* initialized magic number */
	short		inpsiz, outsiz;	/* input & output pixel sizes */
	float		inplogo;	/* input offset for log encoding */
	float		inpgam;		/* input gamma value or ln(base) */
	float		fromByte[256];	/* 8-bit de-gamma/exp table */
	int		usemat;		/* need matrix? */
	void		(*getf)(float *drgb, const uby8 *sp, int n,
					const struct _ColorCv *);
	COLORMAT	cmat;		/* color conversion matrix */
	float		outlogo;	/* output offset for log encoding */
	float		outgam;		/* 1/output gamma */
	uby8		toByte[GTSIZ];	/* 8-bit re-gamma table */
	float		offset;		/* amount to add before output */
	void		(*putf)(uby8 *dp, const float *srgb, int n,
					const struct _ColorCv *);
} ColorCv;

static void
setGammas(ColorCv *cvp)
{
	int	i;

	if (cvp->inplogo > 0)		/* log decoding */
		for (i = 256; i--; )
			cvp->fromByte[i] = cvp->inplogo *
				expf(cvp->inpgam*(i+.5f)*(float)(M_LN10/256.));
	else if (cvp->inpgam == 1.f)	/* special case for speed */
		for (i = 256; i--; )
			cvp->fromByte[i] = (i+.5f)*(1.f/256.f);
	else				/* gamma decoding */
		for (i = 256; i--; )
			cvp->fromByte[i] = powf((i+.5f)/256.f, cvp->inpgam);
	
	if (cvp->outlogo > 0)
		return;
					/* gamma encoding */
	if (cvp->outgam == 1.f)
		for (i = GTSIZ; i--; )
			cvp->toByte[i] = (int)((i+.5f)*(256.f/GTSIZ));
	else
		for (i = GTSIZ; i--; )
			cvp->toByte[i] = (int)(256.f*powf((i+.5f)/(float)GTSIZ,
								cvp->outgam));
}

static void
applyMat(float *sco, const float *sci, int n, const COLORMAT cm)
{
					/* modified from colortrans() */
	for ( ; n-- > 0; sci += 3) {
		*sco++ = cm[0][0]*sci[0] + cm[0][1]*sci[1] + cm[0][2]*sci[2];
		*sco++ = cm[1][0]*sci[0] + cm[1][1]*sci[1] + cm[1][2]*sci[2];
		*sco++ = cm[2][0]*sci[0] + cm[2][1]*sci[1] + cm[2][2]*sci[2];
	}
}

static void
Gfrom3float(float *drgb, const uby8 *sp, int n, const ColorCv *cvp)
{
	const float	*sfp = (const float *)sp;
	
	if (cvp->inplogo > 0) {		/* log decoding */
		const float	mult = cvp->inpgam*M_LN10;
		for (n *= 3; n-- > 0; )
			*drgb++ = cvp->inplogo * expf(*sfp++ * mult);
		return;
	}
	if (gamLINEAR(cvp->inpgam)) {	/* linear decoding */
		memcpy(drgb, sfp, sizeof(float)*3*n);
		return;
	}
					/* gamma decoding */
	for (n *= 3; n-- > 0; sfp++)
		*drgb++ = powf(mkPOS(*sfp), cvp->inpgam);
}

static void
Gfrom1float(float *drgb, const uby8 *sp, int n, const ColorCv *cvp)
{
	const float *	sfp = (const float *)sp;
	
	if (cvp->inplogo > 0) {		/* log decoding */
		const float	mult = cvp->inpgam*M_LN10;
		for ( ; n-- > 0; drgb += 3)
			drgb[0] = drgb[1] = drgb[2] =
				cvp->inplogo * expf(*sfp++ * mult);
		return;
	}
	if (gamLINEAR(cvp->inpgam)) {	/* linear decoding */
		while (n-- > 0) {
			*drgb++ = *sfp;
			*drgb++ = *sfp;
			*drgb++ = *sfp++;
		}
		return;
	}
					/* gamma decoding */
	for ( ; n-- > 0; drgb += 3, sfp++)
		drgb[0] = drgb[1] = drgb[2] = powf(mkPOS(*sfp), cvp->inpgam);
}

static void
Gfrom3short(float *drgb, const uby8 *sp, int n, const ColorCv *cvp)
{
	const uint16 *	ssp = (const uint16 *)sp;
	int		i;
	
	if (cvp->inplogo > 0) {		/* log decoding */
		const float	mult = cvp->inpgam * (M_LN10/0x10000);
		for (i = 3*n; i-- > 0; )
			*drgb++ = cvp->inplogo * expf((*ssp++ + .5f) * mult);
		return;
	}
	for (i = 3*n; i-- > 0; )
		*drgb++ = (*ssp++ + .5f) * (1.f/0x10000);
	if (gamLINEAR(cvp->inpgam))	/* linear decoding */
		return;
					/* gamma decoding */
	for (i = 3*n; i-- > 0; ) {
		--drgb;
		*drgb = powf(*drgb, cvp->inpgam);
	}
}

static void
Gfrom1short(float *drgb, const uby8 *sp, int n, const ColorCv *cvp)
{
	const uint16 *	ssp = (const uint16 *)sp;
	
	if (cvp->inplogo > 0) {		/* log decoding */
		const float	mult = cvp->inpgam * (M_LN10/0x10000);
		for ( ; n-- > 0; drgb += 3)
			drgb[0] = drgb[1] = drgb[2] =
				cvp->inplogo * expf((*ssp++ + .5f) * mult);
		return;
	}
	if (gamLINEAR(cvp->inpgam)) {	/* linear decoding */
		for ( ; n-- > 0; drgb += 3)
			drgb[0] = drgb[1] = drgb[2] =
					(*ssp++ + .5f) * (1.f/0x10000);
		return;
	}
					/* gamma decoding */
	for ( ; n-- > 0; drgb += 3)
		drgb[0] = drgb[1] = drgb[2] =
				powf((*ssp++ + .5f) * (1.f/0x10000), cvp->inpgam);
}

static void
Gfrom3byte(float *drgb, const uby8 *sp, int n, const ColorCv *cvp)
{
	for (n *= 3; n-- > 0; )
		*drgb++ = cvp->fromByte[*sp++];
}

static void
Gfrom1byte(float *drgb, const uby8 *sp, int n, const ColorCv *cvp)
{
	for ( ; n-- > 0; drgb += 3)
		drgb[0] = drgb[1] = drgb[2] = cvp->fromByte[*sp++];
}

static void
GfromNormal(float *dvec, const uby8 *sp, int n, const ColorCv *cvp)
{
	const int32 *	snp = (const int32 *)sp;
	FVECT		vec;

	for ( ; n-- > 0; dvec += 3) {
		decodedir(vec, *snp++);
		VCOPY(dvec, vec);
	}
}

static void
Pto3float(uby8 *dp, const float *srgb, int n, const ColorCv *cvp)
{
	float *	dfp = (float *)dp;
	
	if (cvp->outlogo > 0) {		/* log encoding */
		const float	mult = 1.f/cvp->outlogo;
		for (n *= 3; n-- > 0; srgb++)
			if (*srgb > 0)
				*dfp++ = log10f(*srgb * mult) * cvp->outgam;
			else
				*dfp++ = -30.f * cvp->outgam;
		return;
	}
	if (gamLINEAR(cvp->outgam)) {	/* linear encoding */
		memcpy(dfp, srgb, sizeof(float)*3*n);
		return;
	}
					/* gamma encoding */
	for (n *= 3; n-- > 0; srgb++)
		*dfp++ = powf(mkPOS(*srgb), cvp->outgam);
}

static void
Pto1float(uby8 *dp, const float *srgb, int n, const ColorCv *cvp)
{
	float *	dfp = (float *)dp;
	
	++srgb;			/* Y channel */
	if (cvp->outlogo > 0) {		/* log encoding */
		const float	mult = 1.f/cvp->outlogo;
		for ( ; n-- > 0; srgb += 3)
			if (*srgb > 0)
				*dfp++ = log10f(*srgb * mult) * cvp->outgam;
			else
				*dfp++ = -30.f * cvp->outgam;
		return;
	}
	if (gamLINEAR(cvp->outgam)) {	/* linear encoding */
		for ( ; n-- > 0; srgb += 3)
			*dfp++ = *srgb;
		return;
	}
					/* gamma encoding */
	for ( ; n-- > 0; srgb += 3)
		*dfp++ = powf(mkPOS(*srgb), cvp->outgam);
}

static void
Pto3short(uby8 *dp, const float *srgb, int n, const ColorCv *cvp)
{
	uint16 *	dsp = (uint16 *)dp;
	
	if (cvp->outlogo > 0) {		/* log encoding */
		const float	mult1 = 1.f/cvp->outlogo;
		const float	mult2 = cvp->outgam*0x10000;
		int		iv;
		for (n *= 3; n-- > 0; srgb++)
			if (*srgb > cvp->outlogo) {
				iv = log10f(*srgb * mult1) * mult2;
				*dsp++ = (iv >= 0xffff) ? 0xffff : iv*(iv>0);
			} else
				*dsp++ = 0;
		return;
	}
	if (gamLINEAR(cvp->outgam)) {	/* linear encoding */
		for (n *= 3; n-- > 0; srgb++)
			*dsp++ = (*srgb <= 0) ? 0 : (*srgb >= 1.f) ? 0xffff :
					(int)(*srgb * 65536.f);
		return;
	}
					/* gamma encoding */
	for (n *= 3; n-- > 0; srgb++)
		*dsp++ = (*srgb <= 0) ? 0 : (*srgb >= 1.f) ? 0xffff :
				(int)(powf(*srgb, cvp->outgam) * 65536.f);
}

static void
Pto1short(uby8 *dp, const float *srgb, int n, const ColorCv *cvp)
{
	uint16 *	dsp = (uint16 *)dp;
	
	++srgb;			/* Y channel */
	if (cvp->outlogo > 0) {		/* log encoding */
		const float	mult1 = 1.f/cvp->outlogo;
		const float	mult2 = cvp->outgam*0x10000;
		int		iv;
		for ( ; n-- > 0; srgb += 3)
			if (*srgb > cvp->outlogo) {
				iv = log10f(*srgb * mult1) * mult2;
				*dsp++ = (iv >= 0xffff) ? 0xffff : iv*(iv>0);
			} else
				*dsp++ = 0;
		return;
	}
	if (gamLINEAR(cvp->outgam)) {	/* linear encoding */
		for ( ; n-- > 0; srgb += 3)
			*dsp++ = (*srgb <= 0) ? 0 : (*srgb >= 1.f) ? 0xffff :
					(int)(*srgb * 65536.f);
		return;
	}
					/* gamma encoding */
	for ( ; n-- > 0; srgb += 3)
		*dsp++ = (*srgb <= 0) ? 0 : (*srgb >= 1.f) ? 0xffff :
				(int)(powf(*srgb, cvp->outgam) * 65536.f);
}

static void
Pto3byte(uby8 *dp, const float *srgb, int n, const ColorCv *cvp)
{
	if (cvp->outlogo > 0) {		/* log encoding */
		const float	mult1 = 1.f/cvp->outlogo;
		const float	mult2 = cvp->outgam*256;
		int		iv;
		for (n *= 3; n-- > 0; srgb++)
			if (*srgb > cvp->outlogo) {
				iv = log10f(*srgb * mult1) * mult2;
				*dp++ = (iv >= 255) ? 255 : iv*(iv>0);
			} else
				*dp++ = 0;
		return;
	}
					/* gamma encoding */
	for (n *= 3; n-- > 0; srgb++)
		*dp++ = applyGamma8(*srgb, cvp);
}

static void
Pto1byte(uby8 *dp, const float *srgb, int n, const ColorCv *cvp)
{
	++srgb;			/* Y channel */
	if (cvp->outlogo > 0) {		/* log encoding */
		const float	mult1 = 1.f/cvp->outlogo;
		const float	mult2 = cvp->outgam*256;
		int		iv;
		for ( ; n-- > 0; srgb += 3)
			if (*srgb > cvp->outlogo) {
				iv = log10f(*srgb * mult1) * mult2;
				*dp++ = (iv >= 255) ? 255 : iv*(iv>0);
			} else
				*dp++ = 0;
		return;
	}
					/* gamma encoding */
	for ( ; n-- > 0; srgb += 3)
		*dp++ = applyGamma8(*srgb, cvp);
}

static void
PtoNormal(uby8 *dp, const float *svec, int n, const ColorCv *cvp)
{
	int32 *	dsp = (int32 *)dp;
	FVECT	vec;
	
	for ( ; n-- > 0; svec += 3) {
		VCOPY(vec, svec);
		*dsp++ = (normalize(vec) == 0) ? 0 : encodedir(vec);
	}
}

/* Set up a color conversion */
static int
setColorConv(ColorCv *cvp, const ImgColorSpace *cdp,
			const ImgColorSpace *csp, float sf, int cv_wp)
{
	COLORMAT	tmat;

	cvp->magic = ~CMAGIC;
	if ((sf <= 0) & (cdp->logorig > 0) || (sf < 0 &&
			(cdp->dtype != IDTfloat) & (cdp->dtype != IDTint))) {
		DMESGF(DMCparameter, "Illegal scale factor (%e)", sf);
		return false;
	}
	if ((csp->gamma <= 0) | (cdp->gamma <= 0)) {
		DMESG(DMCparameter, "Illegal gamma");
		return false;
	}
	cvp->inpsiz = ImgPixelSize(csp);
	cvp->outsiz = ImgPixelSize(cdp);
	cvp->inplogo = csp->logorig;
	cvp->inpgam = csp->gamma;
	cvp->getf = NULL;
	switch (csp->dtype) {
	case IDTfloat:
		if (ImgPixelLen[csp->format] == 3)
			cvp->getf = Gfrom3float;
		else if (ImgPixelLen[csp->format] == 1)
			cvp->getf = Gfrom1float;
		break;
	case IDTushort:
		if (ImgPixelLen[csp->format] == 3)
			cvp->getf = Gfrom3short;
		else if (ImgPixelLen[csp->format] == 1)
			cvp->getf = Gfrom1short;
		break;
	case IDTubyte:
		if (ImgPixelLen[csp->format] == 3)
			cvp->getf = Gfrom3byte;
		else if (ImgPixelLen[csp->format] == 1)
			cvp->getf = Gfrom1byte;
		break;
	case IDTint:
		if (csp->format == IPFn)
			cvp->getf = GfromNormal;
		break;
	}
	if (cvp->getf == NULL)
		goto nodice;
	if (cv_wp)
		cv_wp = !PmatchColorSpace(cdp, csp, PICMwhite);
	memset(cvp->cmat, 0, sizeof(cvp->cmat));
	cvp->cmat[0][0] = cvp->cmat[1][1] = cvp->cmat[2][2] = 1.f;
	cvp->usemat = 0;
	switch (csp->format) {
	case IPFy:
	case IPFa:
	case IPFd:
		break;
	case IPFxyz:
		break;		/* XXX ignoring white point as well */
	case IPFrgb:
		if (!colorprimsOK((RGBPRIMP)csp->chroma) ||
			!(cv_wp ? comprgb2xyzWBmat(cvp->cmat, (RGBPRIMP)csp->chroma)
				: comprgb2xyzmat(cvp->cmat, (RGBPRIMP)csp->chroma))) {
			DMESG(DMCwarning, "Bad color conversion from source");
			comprgb2xyzmat(cvp->cmat, stdprims);
		}
		cvp->usemat = 1;
		break;
	case IPFn:
	case IPFvec3:
		if ((cdp->format != IPFvec3) & (cdp->format != IPFn))
			cvp->usemat = -1;
		break;
	default:
		cvp->usemat = -1;
		break;
	}
	if (cvp->usemat < 0)
		goto nodice;
	switch (cdp->format) {
	case IPFy:
	case IPFa:		/* treat alpha, Y, depth as interchangeable */
	case IPFd:
		break;
	case IPFxyz:
		break;		/* XXX ignoring white point as well */
	case IPFrgb:
		if (!colorprimsOK((RGBPRIMP)cdp->chroma) ||
			!(cv_wp ? compxyz2rgbWBmat(tmat, (RGBPRIMP)cdp->chroma)
				: compxyz2rgbmat(tmat, (RGBPRIMP)cdp->chroma))) {
			DMESG(DMCwarning, "Bad color conversion to destination");
			compxyz2rgbmat(tmat, stdprims);
		}
		multcolormat(cvp->cmat, cvp->cmat, tmat);
		cvp->usemat = !(cvp->usemat && PmatchColorSpace(cdp, csp, PICMprims));
		break;
	case IPFn:
	case IPFvec3:
		if ((csp->format != IPFvec3) & (csp->format != IPFn))
			cvp->usemat = -1;
		break;
	default:
		cvp->usemat = -1;
		break;
	}
	if (cvp->usemat < 0)
		goto nodice;
	cvp->outlogo = cdp->logorig;
	cvp->outgam = 1./cdp->gamma;
	cvp->putf = NULL;
	switch (cdp->dtype) {
	case IDTfloat:
		if (ImgPixelLen[cdp->format] == 3)
			cvp->putf = Pto3float;
		else if (ImgPixelLen[cdp->format] == 1)
			cvp->putf = Pto1float;
		break;
	case IDTushort:
		if (ImgPixelLen[cdp->format] == 3)
			cvp->putf = Pto3short;
		else if (ImgPixelLen[cdp->format] == 1)
			cvp->putf = Pto1short;
		break;
	case IDTubyte:
		if (ImgPixelLen[cdp->format] == 3)
			cvp->putf = Pto3byte;
		else if (ImgPixelLen[cdp->format] == 1)
			cvp->putf = Pto1byte;
		break;
	case IDTint:
		if (cdp->format == IPFn)
			cvp->putf = PtoNormal;
		break;
	}
	if (cvp->putf == NULL)
		goto nodice;
	if (!cvp->usemat & (cvp->inplogo <= 0) & (cvp->outlogo <= 0)) {
						/* combine gammas */
		const float	com_gamma = cvp->inpgam * cvp->outgam;
		if (com_gamma >= 1.f) {		/* decode only */
			cvp->inpgam = com_gamma;
			cvp->outgam = 1.f;
			if (sf > 0)
				sf = powf(sf, 1.f/cdp->gamma);
			else if (sf < 0)
				sf = -powf(-sf, 1.f/cdp->gamma);
		} else {			/* encode only */
			cvp->inpgam = 1.f;
			cvp->outgam = com_gamma;
			if (sf > 0)
				sf = powf(sf, 1.f/csp->gamma);
			else if (sf < 0)
				sf = -powf(-sf, 1.f/csp->gamma);
		}
	}
	cvp->offset = 0;
	if (cvp->inplogo <= 0) {		/* apply scale factor */
		int	i, j;
		for (j = 3; j--; )
			for (i = 3; i--; )
				cvp->cmat[j][i] *= sf;
		cvp->usemat |= (0.999f > sf) | (sf > 1.001f);
	} else if (!cvp->usemat & (cvp->outlogo > 0)) {
						/* side-step log transcoding */
		cvp->inplogo = cvp->outlogo = 0;
		cvp->inpgam = cvp->outgam = 1.f;
		cvp->offset = log10f(sf*csp->logorig/cdp->logorig) / cdp->gamma;
		cvp->cmat[0][0] = cvp->cmat[1][1] =
				cvp->cmat[2][2] = csp->gamma/cdp->gamma;
		cvp->usemat |= (fabs(1. - cvp->cmat[0][0]) > 0.001);
	} else
		cvp->inplogo *= sf;		/* a little right lie */
	setGammas(cvp);
	cvp->magic = CMAGIC;			/* mark as initialized */
	return true;
nodice:
	DMESG(DMCparameter, "Unsupported color space conversion");
	return false;
}

/* Check for unmanageable image overlap */
static int
badOverlap(const ImgStruct *ib, const ImgStruct *ia)
{
	if (ib->img == ia->img)
		return false;		/* this overlap we can handle */
	if (!PimagesOverlap(ib, ia))
		return false;		/* no memory overlap -- OK */
	if ((ib->rowsize < 0) ^ (ia->rowsize < 0))
		return false;		/* code doesn't handle this */

					/* check for overstepping */
	return (ib->img < ia->img) ^ (ib->img + (ssize_t)ib->yres*ib->rowsize <
					ia->img + (ssize_t)ia->yres*ia->rowsize);
}

/* Apply color conversion to an array of pixel values */
int
PmapPixels(uby8 *pout, const uby8 *pinp, int n, const void *cv)
{
#define	CSCANLEN	4096
	const ColorCv *	cvp = (const ColorCv *)cv;
	float		scanbuf[3*(CSCANLEN+1)];

	if (n <= 0)
		return (n == 0);
	if ((pout == NULL) | (pinp == NULL) | (cvp == NULL))
		return false;
	if (cvp->magic != CMAGIC) {
		DMESG(DMCparameter, "Invalid color conversion");
		return false;
	}
	if (pout == pinp && (cvp->outsiz > cvp->inpsiz) & (n > CSCANLEN)) {
		DMESG(DMCparameter, "Cannot map pixel buffer onto itself");
		return false;
	}
	while (n > 0) {
		int	len = n;
		if (len > CSCANLEN)
			len = CSCANLEN;
		if (cvp->usemat) {
			(*cvp->getf)(scanbuf+3, pinp, len, cvp);
			applyMat(scanbuf, scanbuf+3, len,
					(const float (*)[3])cvp->cmat);
		} else
			(*cvp->getf)(scanbuf, pinp, len, cvp);
		if (cvp->offset) {
			float *	scp = scanbuf + len*3;
			while (scp > scanbuf)
				*--scp += cvp->offset;
		}
		(*cvp->putf)(pout, scanbuf, len, cvp);
		pinp += len*cvp->inpsiz;
		pout += len*cvp->outsiz;
		n -= len;
	}
	return true;
#undef CSCANLEN
}

/* Scanline color conversion */
static int
convertColor(ImgStruct *ib, const ImgStruct *ia, const ColorCv *cvp)
{
	int		n = ib->yres;
	const uby8 *	sa;
	uby8 *		sb;
					/* check for failure case */
	if (badOverlap(ib, ia)) {
		DMESG(DMCparameter, "Bad image overlap for color conversion");
		return false;
	}
					/* allow overlaps */
	if (ib->img < ia->img || (ib->img == ia->img &&
				  abs(ib->rowsize) <= abs(ia->rowsize))) {
		sb = ib->img;
		sa = ia->img;
		while (n-- > 0) {
			if (!PmapPixels(sb, sa, ib->xres, cvp))
				return false;
			sb += ib->rowsize;
			sa += ia->rowsize;
		}
	} else {
		sb = ib->img + (ssize_t)ib->yres*ib->rowsize;
		sa = ia->img + (ssize_t)ia->yres*ia->rowsize;
		while (n-- > 0) {
			sb -= ib->rowsize;
			sa -= ia->rowsize;
			if (!PmapPixels(sb, sa, ib->xres, cvp))
				return false;
		}
	}
	return true;			/* all done */
}

/* Direct 1-1 copy routine, no color conversion */
static int
straightCopy(ImgStruct *ib, const ImgStruct *ia)
{
	const int	sz = ib->xres*ImgPixelSize(ib->csp);
	int		n = ib->yres;
	const uby8 *	sa;
	uby8 *		sb;
					/* check for failure case */
	if (badOverlap(ib, ia)) {
		DMESG(DMCparameter, "Bad image overlap for copy");
		return false;
	}
					/* allow overlaps */
	if (ib->img < ia->img || (ib->img == ia->img &&
					abs(ib->rowsize) < abs(ia->rowsize))) {
		sb = ib->img;
		sa = ia->img;
		while (n-- > 0) {
			memmove(sb, sa, sz);
			sb += ib->rowsize;
			sa += ia->rowsize;
		}
	} else if (ib->img > ia->img || (ib->img == ia->img &&
					abs(ib->rowsize) > abs(ia->rowsize))) {
		sb = ib->img + (ssize_t)ib->yres*ib->rowsize;
		sa = ia->img + (ssize_t)ia->yres*ia->rowsize;
		while (n-- > 0) {
			sb -= ib->rowsize;
			sa -= ia->rowsize;
			memmove(sb, sa, sz);
		}
	}
	return true;			/* all done */
}

/* General color space conversion/copy routine */
int
PmapImage(ImgStruct *ib, const ImgStruct *ia, float sf)
{
	ColorCv		ccvt;
						/* check format */
	if (ia == NULL || (ia->img == NULL) | (ia->csp == NULL))
		return false;
	if (ib == NULL || ib->csp == NULL)
		return false;
	if (ib->img != NULL) {			/* already allocated? */
		if ((ib->xres != ia->xres) | (ib->yres != ia->yres)) {
			DMESG(DMCparameter, "Mismatched image dimensions");
			return false;
		}
		if (ImgPixelSize(ib->csp)*ib->xres > abs(ib->rowsize)) {
			DMESG(DMCparameter, "New pixel size does not fit");
			return false;
		}
	} else {				/* else need to allocate */
		ib->xres = ia->xres; ib->yres = ia->yres;
		if (!PnewImage(ib, .0))
			return false;
	}
	if (ib->csp->format == IPFn)		/* scaling irrelevant */
		sf = 1.f;
	else if (sf == 0) {			/* optimization */
		DMESGF((ib->csp->logorig > 0) ? DMCwarning : DMCtrace,
			"Scaling %s image by zero", PdescribeCS(ib->csp, NULL));
		PclearImage(ib, NULL);
		return true;
	}
						/* check for direct copy */
	if ((0.999f <= sf) & (sf <= 1.001f) &&
			PmatchColorSpace(ib->csp, ia->csp, PICMequiv)) {
		if (straightCopy(ib, ia))
			return true;
		PfreeImage(ib);
		return false;
	}
						/* perform conversion */
	strcpy(dmessage_buf, "Converting ");
	PdescribeCS(ia->csp, strchr(dmessage_buf,'\0'));
	strcat(dmessage_buf, " to ");
	PdescribeCS(ib->csp, strchr(dmessage_buf,'\0'));
	DMESG(DMCtrace, dmessage_buf);
	if (setColorConv(&ccvt, ib->csp, ia->csp, sf, true) &&
			convertColor(ib, ia, &ccvt))
		return true;
						/* free dest. on failure */
	PfreeImage(ib);
	return false;
}

/* In situ color space conversion (will delink image) */
int
PconvertColorSpace(ImgStruct *ib, const ImgColorSpace *dcs, float sf)
{
	ImgColorSpace	*mbase_csp = NULL;
	ImgColorSpace	csTemp;
	ImgStruct	iorig;
	int		ok;
	
	if (ib == NULL || (ib->img == NULL) | (ib->csp == NULL) | (dcs == NULL))
		return false;
	if (ib->mbase != NULL)
		mbase_csp = (ImgColorSpace *)MobjMem(ib->mbase);
	iorig.img = NULL;
	/* The following line allowed images to shrink (or expand)
	 * with varying pixel size, but it's generally preferable
	 * to fit the image into memory, so I replaced it with
	 * a tight constraint on pixel size (but not scanline size):
	if (ImgPixelSize(dcs)*ib->xres > ib->rowsize ||
	 */
	if (ImgPixelSize(dcs) != ImgPixelSize(ib->csp) ||
			(ib->mbase != NULL && (ib->mbase->nref > 1) |
					(ib->mbase->freo != SysFree))) {
		iorig = *ib;
		ib->img = NULL;			/* need to reallocate */
	} else if (!dcs->cstatic & (dcs != mbase_csp)) {
		if (ib->csp == mbase_csp) {
				/* tricky business to avoid reallocation */
			if (!PlinkImage(&iorig, ib))
				return false;
			PcopyCS(&csTemp, mbase_csp);
			iorig.csp = &csTemp;	/* remember current CS */
			PcopyCS(mbase_csp, dcs);
			dcs = mbase_csp;	/* ib now preset to new CS */
		} else {
			iorig = *ib;
			ib->img = NULL;		/* need to reallocate */
		}
	}
	if (iorig.img == NULL && !PlinkImage(&iorig, ib))
		return false;
	ib->csp = dcs;				/* apply conversion */
	ok = PmapImage(ib, &iorig, sf);
	PfreeImage(&iorig);			/* unlink original */
	return ok;
}

/* Convert image color space, sending scanlines to given callback function */
int
PconvertImageCB(const ImgStruct *ia, const ImgColorSpace *dcs, float sf,
					PscanlineMethod *scb, void *udp)
{
	int		rv, rvsum = 0;
	ColorCv		ccvt;
	uby8 *		destbuf;
	int		y;

	if (ia == NULL || (ia->csp == NULL) | (dcs == NULL) | (scb == NULL))
		return -1;
	if (dcs->format == IPFn)		/* scaling irrelevant */
		sf = 1.f;
						/* short-cut direct case */
	if ((0.999f <= sf) & (sf <= 1.001f) &&
			PmatchColorSpace(dcs, ia->csp, PICMequiv)) {
		for (y = 0; y < ia->yres; y++) {
			rv = (*scb)(ProwPtr(ia, y), ia->xres, y, udp);
			if (rv < 0)
				return rv;
			rvsum += rv;
		}
		return rvsum;
	}
						/* set up conversion */
	if (!setColorConv(&ccvt, dcs, ia->csp, sf, true))
		return -1;
						/* allocate scanline */
	destbuf = (uby8 *)malloc(ccvt.outsiz*ia->xres);
	if (destbuf == NULL) {
		DMESG(DMCmemory, "malloc() failed in PconvertImageCB");
		return -1;
	}
						/* perform scan conversion */
	for (y = 0; y < ia->yres; y++) {
		if (PmapPixels(destbuf, ProwPtr(ia, y), ia->xres, &ccvt))
			rv = (*scb)(destbuf, ia->xres, y, udp);
		else
			rv = -1;
		if (rv < 0) {
			rvsum = rv;
			break;
		}
		rvsum += rv;
	}
	free(destbuf);				/* clean up and return */
	return rvsum;
}

/* Allocate an opaque color conversion structure */
void *
PcreateColorConv(const ImgColorSpace *cdp, const ImgColorSpace *csp,
			float sf, int cv_wp)
{
	ColorCv *	ccbuf;
	
	if ((cdp == NULL) | (csp == NULL))
		return NULL;
	ccbuf = (ColorCv *)malloc(sizeof(ColorCv));
	if (ccbuf == NULL) {
		DMESG(DMCmemory, "malloc() failed in PcreateColorConv");
		return NULL;
	}
	if (setColorConv(ccbuf, cdp, csp, sf, cv_wp))
		return (void *)ccbuf;
	free(ccbuf);
	return NULL;
}

/* Get/add cached color conversion (no free routine) */
const void *
PgetColorConv(const ImgColorSpace *cdp, const ImgColorSpace *csp)
{
	static struct _ColorConvEnt {
		struct _ColorConvEnt	*next;
		const ImgColorSpace	*cdp;
		const ImgColorSpace	*csp;
		ColorCv			ccvt;
	}			*ccList = NULL;
	struct _ColorConvEnt	*cent, *cdmatch = NULL;
	ImgColorSpace		*newCS;

	if ((cdp == NULL) | (csp == NULL))
		return NULL;
					/* look for entry in cache */
	for (cent = ccList; cent != NULL; cent = cent->next) {
		if (cdmatch == NULL && PmatchColorSpace(cdp, cent->cdp, PICMequiv)) {
			cdmatch = cent;
			if (cdp->cstatic && !cdmatch->cdp->cstatic) {
				for (cent = cdmatch->next; cent != NULL &&
						cent->cdp == cdmatch->cdp;
							cent = cent->next)
					cent->cdp = cdp;
				free((ImgColorSpace *)cdmatch->cdp);
				(cent = cdmatch)->cdp = cdp;
			}
		}
		if (cdmatch != NULL) {
			if (cent->cdp != cdmatch->cdp)
				break;
			if (PmatchColorSpace(csp, cent->csp, PICMequiv)) {
				if (csp->cstatic && !cent->csp->cstatic) {
					free((ImgColorSpace *)cent->csp);
					cent->csp = csp;
				}
				return (const void *)&cent->ccvt;
			}
		}
	}
					/* need to create new entry */
	cent = (struct _ColorConvEnt *)malloc(sizeof(struct _ColorConvEnt));
	if (cent == NULL)
		goto memerr;
	if (!setColorConv(&cent->ccvt, cdp, csp, 1.f, true)) {
		free(cent);
		return NULL;
	}
	if (cdmatch != NULL) {		/* insert in list next to siblings */
		cent->cdp = cdmatch->cdp;
		cent->next = cdmatch->next;
		cdmatch->next = cent;
	} else {			/* no siblings, so push on top */
		if (cdp->cstatic) {
			cent->cdp = cdp;
		} else {
			newCS = (ImgColorSpace *)malloc(sizeof(ImgColorSpace));
			if (newCS == NULL)
				goto memerr;
			PcopyCS(newCS, cdp);
			cent->cdp = newCS;
		}
		cent->next = ccList;
		ccList = cent;
	}
	if (csp->cstatic) {
		cent->csp = csp;
	} else {
		newCS = (ImgColorSpace *)malloc(sizeof(ImgColorSpace));
		if (newCS == NULL)
			goto memerr;
		PcopyCS(newCS, csp);
		cent->csp = newCS;
	}
	return (const void *)&cent->ccvt;
memerr:
	DMESG(DMCmemory, "malloc() failed in PgetColorConv");
	return NULL;
}

/* PconvertPixel:
 *  Convert color for a single pixel.
 *  This is a rather expensive operation since our routines are
 *  all optimized for converting a large numbr of pixels using
 *  tables, etc.  However, we cache conversion structures to
 *  avoid set-up costs for repeat operations.
 */
PixelVal
PconvertPixel(PixelVal pval, const ImgColorSpace *cdp)
{
	if ((cdp == NULL) | (pval.csp == NULL))
		return Pblack;
	if (PmatchColorSpace(cdp, pval.csp, PICMequiv)) {
		pval.csp = cdp;		/* expected side-effect */
		return pval;
	}
	if (!PmapPixels(pval.v.b, pval.v.b, 1, PgetColorConv(cdp, pval.csp)))
		return Pblack;
	pval.csp = cdp;
	return pval;
}

/* compute distance-squared to black-body curve (approx. as line) */
static float
blackBodyWeighting(float U, float V)
{
	const float	blueU=0.2871f, blueV=0.3828f,
				redU=0.4517f, redV=0.2034f;
	const float	d = (blueU-redU)*(blueU-redU)
					+ (blueV-redV)*(blueV-redV);
	const float	halfDist2 = d/(6.f*6.f);
	float		d1, d2, dist2 = -1.f;
	
	d1 = (U-blueU)*(U-blueU) + (V-blueV)*(V-blueV);
	d2 = (U-redU)*(U-redU) + (V-redV)*(V-redV);
	
	if (d2 > d1) {			/* check if past endpoints */
		if (d2 - d1 > d)
			dist2 = d1;
	} else {
		if (d1 - d2 > d)
			dist2 = d2;
	}
	if (dist2 <= 0) {
		d1 = d + d2 - d1;
		dist2 = d2 - 0.25f*d1*d1/d;
	}
	return halfDist2 / (halfDist2 + dist2);
}

/* In situ, modified gray-world white balance */
int
PautoWhiteBal(ImgStruct *ib, float frac, int orig)
{
	const float	minWeight = 0.1f;
	extern COLORMAT	vkmat, ivkmat;
	int		ok;
	ImgStruct	iTmp;
	ColorCv		ccvt;
	float *		fpp;
	COLOR		rgb_corr, rgb_whtorig, xyz_val;
	ImgColorSpace	adjSharpCS;
	float		U, V, w, U_sum, V_sum, w_sum;
	int		x, y;
	
	if (ib == NULL || (ib->img == NULL) | (ib->csp == NULL))
		return false;
	if (frac <= 0.01f)
		return true;
	if (ImgPixelLen[ib->csp->format] != 3) {
		DMESG(DMCparameter, "Input to PautoWhiteBalance not in color");
		return false;
	}
	if (!setColorConv(&ccvt, &ICS_SharpRGB, ib->csp, 1.f, false))
		return false;
	iTmp.img = NULL;
	iTmp.csp = &ICS_SharpRGB;
	if (!PnewImage(&iTmp, .0))
		return false;
	U_sum = V_sum = w_sum = 0;		/* BB-biased white avg. */
	for (y = 0; y < iTmp.yres; y++) {
		fpp = (float *)ProwPtr(&iTmp,y);
		if (!PmapPixels((uby8 *)fpp, ProwPtr(ib,y), ib->xres, &ccvt)) {
			w_sum = 0;
			break;
		}
		for (x = 0; x < iTmp.xres; x++, fpp += 3) {
			w = fpp[0] + fpp[1] + fpp[2];
			if (w <= 0)
				continue;
			U = fpp[0]/w;
			V = fpp[1]/w;
			w = blackBodyWeighting(U, V);
			if (w < minWeight)
				continue;
			U_sum += w*U; V_sum += w*V; w_sum += w;
		}
	}
	ok = (w_sum > minWeight*iTmp.xres*iTmp.yres/1000.f);
	if (ok) {				/* shift white */
		U = U_sum/w_sum; V = V_sum/w_sum;
		rgb_corr[0] = U/V; rgb_corr[1] = 1.f;
		rgb_corr[2] = (1.f - U)/V - 1.f;
		colortrans(xyz_val, ivkmat, rgb_corr);
		xyz_val[0] /= xyz_val[1]; xyz_val[2] /= xyz_val[1];
		xyz_val[1] = 1.f;
		sprintf(dmessage_buf,
	"Shifting white from (%.4f,%.4f) %.0f%% to estimated WPT (%.4f,%.4f)",
				ib->csp->chroma[3][0], ib->csp->chroma[3][1],
				100.f*frac,
				xyz_val[0]/(xyz_val[0]+1.f+xyz_val[2]),
				1.f/(xyz_val[0]+1.f+xyz_val[2]));
		DMESG(DMCtrace, dmessage_buf);
		colortrans(rgb_corr, vkmat, xyz_val);
		xyz_val[0] = ib->csp->chroma[3][0]/ib->csp->chroma[3][1];
		xyz_val[1] = 1.f;
		xyz_val[2] = (1.f - ib->csp->chroma[3][0])/ib->csp->chroma[3][1] - 1.f;
		colortrans(rgb_whtorig, vkmat, xyz_val);
		rgb_corr[0] = frac*rgb_corr[0] + (1.f-frac)*rgb_whtorig[0];
		rgb_corr[1] = frac*rgb_corr[1] + (1.f-frac)*rgb_whtorig[1];
		rgb_corr[2] = frac*rgb_corr[2] + (1.f-frac)*rgb_whtorig[2];
		colortrans(xyz_val, ivkmat, rgb_corr);
						/* transform it back */
		PcopyCS(&adjSharpCS, &ICS_SharpRGB);
		w = xyz_val[0] + xyz_val[1] + xyz_val[2];
		adjSharpCS.chroma[3][0] = xyz_val[0]/w;
		adjSharpCS.chroma[3][1] = xyz_val[1]/w;
		iTmp.csp = &adjSharpCS;		/* we make this true below */
		if (!orig) {			/* use new white point */
			ImgColorSpace	finalCS;
			PcopyCS(&finalCS, ib->csp);
			finalCS.chroma[3][0] = adjSharpCS.chroma[3][0];
			finalCS.chroma[3][1] = adjSharpCS.chroma[3][1];
			PfreeImage(ib);		/* need to store new CS */
			ib->csp = &finalCS;
			ok = PnewImage(ib, .0);
		}
		if (ok)
			ok = setColorConv(&ccvt, ib->csp, iTmp.csp, 1.f, orig);
		if (ok) {
			COLORMAT	rat;	/* prepend white scaling */
			memset(rat, 0, sizeof(rat));
			rat[0][0] = 1.f/rgb_corr[0];
			rat[1][1] = 1.f/rgb_corr[1];
			rat[2][2] = 1.f/rgb_corr[2];
			multcolormat(ccvt.cmat, rat, ccvt.cmat);
			ccvt.usemat = 1;
			ok = convertColor(ib, &iTmp, &ccvt);
		}
	}
	PfreeImage(&iTmp);
	return ok;
}
