/*
 *  jhtonemap2.c
 *
 *  Created by Greg Ward on 9/7/10.
 *  Copyright 2010 Dolby Laboratories. All rights reserved.
 *
 *  Advanced tone-mapping operator for JPEG-HDR compression.
 */

#define JPEG_INTERNALS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "pimage.h"
#include "jpeghdr.h"

#ifndef M_LN2
#define M_LN2           0.693147180559945309417232121458176568
#endif
#ifndef M_LN10
#define M_LN10		2.30258509299404568401799145468436421
#endif

#ifndef DISP_LN_RNG
#define	DISP_LN_RNG	3.8f			/* log(display_range) */
#endif

static const float	gamtGoal = 99.5f;	/* HDR in-gamut goal */
static const float	rendGoal = 98.f;	/* render percentile goal */
static const JHSAMPLE	rendGamt[3] = {1,1,1};	/* upper rendering gamut */

#if defined(_WIN32) || defined(_WIN64)

#ifdef _MSC_VER
#include <BaseTsd.h>
typedef SSIZE_T	ssize_t;
#endif
#undef powf
#define powf(x,e)	(float)pow((double)(x),(double)(e))
#undef logf
#define logf(x)		(float)log((double)(x))
#undef expf
#define expf(x)		(float)exp((double)(x))
#undef log10f
#define log10f(x)	(float)log10((double)(x))
#undef exp10f
#define exp10f(x)	(float)exp(M_LN10*(double)(x))

#else

#define exp10f(x)	expf((float)M_LN10*(float)(x))

#endif

/* Compute log luminance image and histogram */
LOCAL(void)
get_log_image (jh_compress_ptr jhinf, ImgStruct *limg)
{
	static const ImgColorSpace	logSpace = {
		IDTushort, IPFy, TRUE, 20.f, JH_LUM_MIN,
		{{1./3., 1./3.}}
	};
	const float	logmin = log10f(JH_LUM_MIN);
	unsigned short	*slp;
	int		l16min, l16max;
	JHSAMPLE	pv[3];
	float		lv;
	int		i, j, ystep;
					/* get log luminance image */
	limg->csp = &logSpace;
	limg->xres = jhinf->cinfo.image_width;
	limg->yres = jhinf->cinfo.image_height;
	limg->img = NULL;
	if (!PnewImage(limg, .0))
		ERREXIT1(&jhinf->cinfo, JERR_OUT_OF_MEMORY, 217);
	l16min = 0xffff; l16max = 0;
	for (j = 0; j < limg->yres; j++) {
		slp = (unsigned short *)ProwPtr(limg, j);
		for (i = 0; i < limg->xres; i++, slp++) {
			hdr_pval(pv, &jhinf->hdr, i, j);
			lv = hdr_lum(pv);
			if (lv <= JH_LUM_MIN) {
				*slp = 0;
				continue;
			}
			lv = log10f(lv) - logmin;
			if (lv >= logSpace.gamma) {
				*slp = 0xffff;
				continue;
			}
			*slp = (float)(0xffff/logSpace.gamma) * lv;
			if (*slp < l16min)
				l16min = *slp;
			if (*slp > l16max)
				l16max = *slp;
		}
	}
					/* create histogram */
	if (jhinf->hist_len > 0)
		free((void *)jhinf->histo);
	jhinf->hist_len = (JH_HIST_SIZ * (l16max - l16min)) >> 16;
	jhinf->histo = (long *)calloc(jhinf->hist_len, sizeof(long));
	if (jhinf->histo == NULL)
		ERREXIT1(&jhinf->cinfo, JERR_OUT_OF_MEMORY, 218);
	jhinf->hist_lminmax[0] = (l16min*logSpace.gamma/0xffff + logmin)*M_LN10;
	jhinf->hist_lminmax[1] = (l16max*logSpace.gamma/0xffff + logmin)*M_LN10;
	ystep = 1;
	while ((unsigned long long)limg->yres/ystep*limg->xres !=
			(unsigned long)(limg->yres/ystep*limg->xres))
		ystep++;
	jhinf->histot = 0;
	for (j = 0; j < limg->yres; j += ystep) {
		slp = (unsigned short *)ProwPtr(limg, j);
		for (i = 0; i < limg->xres; i++, slp++) {
			if (!*slp) continue;
			if (*slp == 0xffff) continue;
			jhinf->histo[(*slp - l16min) * jhinf->hist_len
						/ (l16max+1 - l16min)]++;
			jhinf->histot++;
		}
	}
	j = jhinf->histot / 200;	/* raise floor above noise */
	for (i = 0; j >= 0; i++)
		j -= jhinf->histo[i];
	l16min += (i-1)*(l16max - l16min)/jhinf->hist_len;
	for (j = 0; j < limg->yres; j += ystep) {
		slp = (unsigned short *)ProwPtr(limg, j);
		for (i = 0; i < limg->xres; i++, slp++)
			if (*slp < l16min)
				*slp = l16min;
	}
}

/* Compute histogram adjustment tone-mapping operator */
LOCAL(boolean)
comp_tmap (jh_compress_ptr jhinf, unsigned short tmap[],
				const ImgColorSpace *csp)
{
	unsigned long	*histo;
	unsigned long	histot;
	unsigned long	trimmed, tol, ceilc;
	int		dispmin, dispmax;
	int		i, j;
	double		mult, offs, bx, nf;
					/* check for within display range */
	if (jhinf->hist_lminmax[1] - jhinf->hist_lminmax[0] <= DISP_LN_RNG)
		return FALSE;
					/* copy luminance histogram */
	if (jhinf->histo == NULL)
		return FALSE;
	histo = (unsigned long *)malloc(sizeof(unsigned long)*(jhinf->hist_len + 1));
	if (histo == NULL)
		return FALSE;
	*histo++ = 0;			/* sneaky [-1] entry */
	memcpy((void *)histo, (void *)jhinf->histo,
			sizeof(unsigned long)*jhinf->hist_len);
	histot = jhinf->histot;
					/* trim super-linear bin counts */
	tol = histot / 100;
	do {
		ceilc = 0.5 + (double)histot *
			(jhinf->hist_lminmax[1] - jhinf->hist_lminmax[0]) /
			((double)jhinf->hist_len * DISP_LN_RNG);
		trimmed = 0;
		for (i = 0; i < jhinf->hist_len; i++)
			if (histo[i] > ceilc) {
				trimmed += histo[i] - ceilc;
				histo[i] = ceilc;
			}
		histot -= trimmed;
	} while (trimmed > tol);
					/* convert to cumulative */
	for (i = 1; i < jhinf->hist_len; i++)
		histo[i] += histo[i-1];
					/* equalization mapping */
	mult = csp->gamma*(M_LN10/(1<<16)) * jhinf->hist_len /
			(jhinf->hist_lminmax[1] - jhinf->hist_lminmax[0]);
	offs = (log(csp->logorig) - jhinf->hist_lminmax[0]) * jhinf->hist_len /
			(jhinf->hist_lminmax[1] - jhinf->hist_lminmax[0]);
	dispmax = -log10f(csp->logorig)/csp->gamma * (1<<16);
	dispmin = dispmax - DISP_LN_RNG/M_LN10*(1<<16)/csp->gamma;
	nf = (double)(dispmax - dispmin) / histot;
	for (j = 0; j < (1<<16); j++) {
		i = bx = (j+.5)*mult + offs;
		if (bx <= 0) {
			tmap[j] = dispmin;
			continue;
		}
		if (i >= jhinf->hist_len) {
			tmap[j] = dispmax;
			continue;
		}
		bx -= (double)i;
		tmap[j] = ((1.-bx)*histo[i-1] + bx*histo[i])*nf + dispmin + .5;
	}
	--histo;			/* done with modified histogram */
	free((void *)histo);
	return TRUE;
}

/* Get gamma for this level in the ratio image pyramid */
/* The sum from i=1 to N should equal 1.0 */
LOCAL(double)
get_level_gamma (int i, int N)
{
	return i * 2 / (double)(N*(N+1));
}

/* Compute multiscale tone adjustment on ratio between global TM and HDR */
LOCAL(boolean)
comp_multiscale (ImgStruct *ratIm, int i, const int N)
{
	const int		l16one = -log10f(ratIm->csp->logorig) /
						ratIm->csp->gamma * (1<<16);
	const double		gamma = get_level_gamma(i, N);
	ImgStruct		ratSmall, ratResamp;
	const unsigned short	*ssp;
	unsigned short		*dsp;
	int			x, y;
					/* check end to recursion */
	if (i <= 1) {
		for (y = 0; y < ratIm->yres; y++) {
			dsp = (unsigned short *)ProwPtr(ratIm,y);
			for (x = 0; x < ratIm->xres; x++, dsp++)
				*dsp = ((int)*dsp - l16one)*gamma +
						l16one + 0.5;
		}
		return TRUE;
	}
	ratSmall = *ratIm;		/* downsample image */
	ratSmall.xres = ratIm->xres >> 1;
	ratSmall.yres = ratIm->yres >> 1;
	ratSmall.img = NULL;
	if (!PsizeImage(&ratSmall, ratIm, PSgaussian))
		return FALSE;
					/* apply recursively */
	if (!comp_multiscale(&ratSmall, i-1, N)) {
		PfreeImage(&ratSmall);
		return FALSE;
	}
	ratResamp = *ratIm;		/* upsample result */
	ratResamp.img = NULL;
	if (!PsizeImage(&ratResamp, &ratSmall, PSlinear))	
		return FALSE;
	PfreeImage(&ratSmall);		/* multiply into original */
	for (y = 0; y < ratIm->yres; y++) {
		ssp = (const unsigned short *)ProwPtr(&ratResamp,y);
		dsp = (unsigned short *)ProwPtr(ratIm,y);
		for (x = 0; x < ratIm->xres; x++, dsp++)
			*dsp = ((int)*dsp - l16one)*gamma +
					(double)*ssp++ + 0.5;
	}
	PfreeImage(&ratResamp);		/* that's all there is to it! */
	return TRUE;
}

/* Apply multiscale local tone-mapping operator (allocates tmi if NULL) */
GLOBAL(void)
jpeghdr_tonemap_multiscale (jh_compress_ptr jhinf)
{
	static const ImgColorSpace	bitmapSpace = {
		IDTubyte, IPFn, TRUE, 1, 0
	};
	float		lumap[1<<16];
	UINT8		dmap[1<<16];
	unsigned short	*tmap = (unsigned short *)lumap;
	ImgStruct	logImg, ratImg, okMap;
	UINT8		*tmpp;
	uby8		*oksp;
	float		gammaRedu, dmin, dmul, pctls[2];
	unsigned short	*slp, minmax16[2];
	int		l16min, l16max, l16one;
	unsigned long	npix2check, npixOK;
	JHSAMPLE	pv[3], sf;
	int		i, j, n;
					/* sanity check */
	if ((jhinf->cinfo.image_width <= 0) | (jhinf->cinfo.image_height <= 0))
		return;
					/* create log Y image & histogram */
	get_log_image(jhinf, &logImg);
					/* get log scale origin */
	l16one = -log10f(logImg.csp->logorig)/logImg.csp->gamma * (1<<16);
					/* compute tone mapping */
	if (!comp_tmap(jhinf, tmap, logImg.csp)) {
		PfreeImage(&logImg);
		jpeghdr_tonemap_default(jhinf);
		return;
	}
					/* create TM ratio image */
	ratImg = logImg;
	ratImg.img = NULL;
	if (!PrenderImageI(&ratImg, &logImg))
		ERREXIT1(&jhinf->cinfo, JERR_OUT_OF_MEMORY, 219);
	for (j = 0; j < ratImg.yres; j++) {
		slp = (unsigned short *)ProwPtr(&ratImg, j);
		for (i = 0; i < ratImg.xres; i++, slp++)
			*slp = (int)tmap[*slp] - (int)*slp + l16one;
	}
					/* determine pyramid height */
	i = ratImg.xres; j = ratImg.yres;
	for (n = 1; (i > 2) & (j > 2); n++) {
		i >>= 1; j >>= 1;
	}
					/* apply multiscale adjustment */
	if (!comp_multiscale(&ratImg, n, n))
		ERREXIT1(&jhinf->cinfo, JERR_OUT_OF_MEMORY, 220);

					/* multiply into luminance */
	for (j = 0; j < ratImg.yres; j++) {
		const unsigned short	*rp =
				(const unsigned short *)ProwPtr(&ratImg, j);
		slp = (unsigned short *)ProwPtr(&logImg, j);
		for (i = 0; i < ratImg.xres; i++)
			*slp++ += (int)*rp++ - l16one;
	}
					/* done with ratio image */
	PfreeImage(&ratImg);
					/* recompute min & max */
	pctls[0] = 0;
	pctls[1] = rendGoal;
	PcomputePercentiles(minmax16, pctls, 2, &logImg, PHCluminance);
	l16min = (minmax16[0] < tmap[0]) ? minmax16[0] : tmap[0];
	l16max = minmax16[1];
					/* create gamut sampling bitmap */
	okMap.csp = &bitmapSpace;
	okMap.xres = (logImg.xres + 7)>>3;
	okMap.yres = logImg.yres;
	okMap.img = NULL;
	if (!PnewImage(&okMap, .0))
		ERREXIT1(&jhinf->cinfo, JERR_OUT_OF_MEMORY, 221);
					/* make lookup to speed conversions */
	gammaRedu = logImg.csp->gamma / (float)(1<<16);
	for (i = 1<<16; i--; )
		lumap[i] = logImg.csp->logorig*exp10f((i+.5f)*gammaRedu);
	
    for (n = 0; n < 2; n++) {		/* check HDR gamut, then rendered */
	const int	goal = n ? rendGoal : gamtGoal;
	PclearImage(&okMap, NULL);
#if JH_PSMP_MAX > 0			/* reduce gamut-checking time */
	if ((double)logImg.xres*logImg.yres > JH_PSMP_MAX) {
		const int	samp_prob = (double)RAND_MAX*JH_PSMP_MAX /
					((double)logImg.xres*logImg.yres) + .5;
		npix2check = 0;
		for (j = 0; j < logImg.yres; j++) {
			oksp = ProwPtr(&okMap, j);
			for (i = 0; i < logImg.xres; oksp += !(++i&7))
				if (rand() <= samp_prob)
					++npix2check;
				else
					*oksp |= 1 << (i&7);
		}
	} else
#endif
		npix2check = logImg.xres*logImg.yres;

	npixOK = 0;
	do {				/* meet in-gamut percentage goal */
		if (npixOK)
			l16max += (goal - 99.9f*npixOK/npix2check) *
					(0.01f*(1<<16))/logImg.csp->gamma + 1;
		for (j = 0; j < logImg.yres; j++) {
			oksp = ProwPtr(&okMap, j);
			slp = (unsigned short *)ProwPtr(&logImg, j);
			for (i = 0; i < logImg.xres; oksp += !(++i&7), slp++) {
				int	bitv;
				while (*oksp == 0xff && i+8 < logImg.xres) {
					bitv = (i+8) & ~07;
					slp += bitv-i;
					i = bitv;
					++oksp;
				}
				bitv = 1 << (i&7);
				if (*oksp & bitv)
					continue;
				sf = lumap[(int)*slp - l16max + l16one];
				if (sf <= .072f) {
					*oksp |= bitv;
					++npixOK;
					continue;
				}
				hdr_pval(pv, &jhinf->hdr, i, j);
				sf /= hdr_lum(pv);
				pv[0] *= sf; pv[1] *= sf; pv[2] *= sf;
				if (n ?		pv[0] <= rendGamt[0] &&
						pv[1] <= rendGamt[1] &&
						pv[2] <= rendGamt[2]
					      : jpeghdr_in_gamut(jhinf, pv)) {
					*oksp |= bitv;
					++npixOK;
				}
			}
		}
	} while (npixOK < goal*0.01f*npix2check);
    }					/* done with gamut checks */
	PfreeImage(&okMap);
					/* compute display mapping */
	dmin = lumap[l16min - l16max + l16one];
	dmul = 1.f/(1.f - dmin);
	gammaRedu = 1.f/jhinf->gamma;
	for (i = 1<<16; i-- > l16max; )
		dmap[i] = 255;
	while (i > l16min)
		dmap[i--] = 256. * powf(
				dmul*(lumap[i - l16max + l16one] - dmin),
				gammaRedu );
	while (i >= 0)
		dmap[i--] = 0;

	jpeghdr_alloc_tmi(jhinf);	/* final tone-mapped 8-bit image */
	tmpp = jhinf->tmi;
	for (j = 0; j < jhinf->cinfo.image_height; j++) {
		slp = (unsigned short *)ProwPtr(&logImg, j);
		for (i = jhinf->cinfo.image_width; i--; )
			*tmpp++ = dmap[*slp++];
	}
	PfreeImage(&logImg);		/* all done -- clean up */
}
