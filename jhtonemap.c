/*
 *  jhtonemap.c
 *
 *  Created by Greg Ward on 10/6/04.
 *  Copyright 2004 Sunnybrook Technologies, <www.sunnybrooktech.com>. 
 *  All rights reserved.
 *
 *  Basic routines for tone-mapping & resampling HDR JPEG images.
 */

#define JPEG_INTERNALS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "jpeghdr.h"

#if defined(_WIN32) || defined(_WIN64)
#undef powf
#define powf(x,e)	(float)pow((double)(x),(double)(e))
#undef logf
#define logf(x)		(float)log((double)(x))
#undef expf
#define expf(x)		(float)exp((double)(x))
#endif

#ifndef M_LN2
#define M_LN2           0.69314718055994530942  /* log e2 */
#endif

#define set_sampling(xr,yr)	const int samp_prob = JH_PSMP_MAX > 0 && \
					JH_PSMP_MAX < (xr)*(yr) ? \
					(int)((double)RAND_MAX*JH_PSMP_MAX / \
						((double)(xr)*(yr))) : 0

#define	skip_sample()		(samp_prob && rand() > samp_prob)

/* Compute natural log histogram */
GLOBAL(void)
jpeghdr_comp_histo (jh_compress_ptr jhinf, int len, float minv, float maxv)
{
	float		hist_sca;
	JHSAMPLE	pv[3];
	float		lv;
	int		i, j;
					/* set up sub-sampling */
	set_sampling(jhinf->cinfo.image_width, jhinf->cinfo.image_height);
					/* reset histogram */
	jhinf->histot = 0;
	if (jhinf->hist_len > 0) {
		free((void *)jhinf->histo);
		jhinf->histo = NULL;
		jhinf->hist_len = 0;
	}
	if ((jhinf->cinfo.image_width <= 0) | (jhinf->cinfo.image_height <= 0))
		return;
	if (minv < JH_LUM_MIN)
		minv = JH_LUM_MIN+1e-20f;
	if (len <= 0 || minv >= maxv)
		return;
	jhinf->histo = (long *)calloc(len, sizeof(long));
	if (jhinf->histo == NULL)
		ERREXIT1(&jhinf->cinfo, JERR_OUT_OF_MEMORY, 199);
	jhinf->hist_len = len;
	jhinf->hist_lminmax[0] = logf(minv);
	jhinf->hist_lminmax[1] = logf(maxv);
	hist_sca = (float)len/(jhinf->hist_lminmax[1] - jhinf->hist_lminmax[0]);
					/* sample image */
	for (j = 0; j < jhinf->cinfo.image_height; j++)
		for (i = 0; i < jhinf->cinfo.image_width; i++) {
			if (skip_sample())
				continue;
			hdr_pval(pv, &jhinf->hdr, i, j);
			lv = hdr_lum(pv);
			if (lv < minv || lv >= maxv)
				continue;
			lv = logf(lv) - jhinf->hist_lminmax[0];
			jhinf->histo[(int)(hist_sca*lv)]++;
			jhinf->histot++;
		}
}

/* Allocate tone-mapped image holder */
GLOBAL(void)
jpeghdr_alloc_tmi (jh_compress_ptr jhinf)
{
	if (jhinf->tmi != NULL)
		return;			/* assume it's OK */
	if ((jhinf->cinfo.image_width <= 0) |
			(jhinf->cinfo.image_height <= 0))
		return;			/* nothing to allocate */
			/* XXX should we call jpeg_get_large()? */
	jhinf->tmi = (UINT8 *)malloc(sizeof(UINT8) *
			jhinf->cinfo.image_width * jhinf->cinfo.image_height);
	if (jhinf->tmi == NULL)
		ERREXIT1(&jhinf->cinfo, JERR_OUT_OF_MEMORY, 200);
}

/* Free allocated tone-mapped image */
GLOBAL(void)
jpeghdr_free_tmi (jh_compress_ptr jhinf)
{
	if (jhinf->tmi == NULL)
		return;
	free((void *)jhinf->tmi);
	jhinf->tmi = NULL;
}

/* Apply default tone-mapping operator (allocates tmi if NULL) */
GLOBAL(void)
jpeghdr_tonemap_default (jh_compress_ptr jhinf)
{
	const float	invgam = 1.f/jhinf->gamma;
	double		lavg = .0;
	unsigned long	nv = 0;
	long		cnt;
	JHSAMPLE	pv[3];
	float		lv, tm;
	UINT8		*tmpp;
	float		lmax, lmin;
	double		alph;
	float		lum_max;
	float		gsca, tsca;
	int		i, j, k;
					/* set up sub-sampling */
	set_sampling(jhinf->cinfo.image_width, jhinf->cinfo.image_height);

					/* sanity check */
	if ((jhinf->cinfo.image_width <= 0) | (jhinf->cinfo.image_height <= 0))
		return;
					/* compute histogram */
	if (!jhinf->hist_len)
		jpeghdr_comp_histo(jhinf, JH_HIST_SIZ, JH_LUM_MIN, JH_LUM_MAX);
					/* compute log min., avg., max. */
	for (k = jhinf->hist_len; k--; ) {
		if (!jhinf->histo[k])
			continue;
		lavg += (double)jhinf->histo[k] * jpeghdr_histo_lval(jhinf, k);
		nv += jhinf->histo[k];
	}
	jpeghdr_alloc_tmi(jhinf);
	if (!nv) {			/* black or out of range image */
		lavg = jpeghdr_histo_lval(jhinf, 0);
		nv = 1;
		lmin = lmax = lavg;
	} else {
		lavg /= (double)nv;
					/* drop bottom 0.5%, top 0.05% */
		cnt = nv/200 + 1;
		for (k = 0; k < jhinf->hist_len; k++)
			if ((cnt -= jhinf->histo[k]) <= 0)
				break;
		lmin = jpeghdr_histo_lval(jhinf, k);
		cnt = nv/2000 + 1;
		for (k = jhinf->hist_len; k--; )
			if ((cnt -= jhinf->histo[k]) <= 0)
				break;
		lmax = jpeghdr_histo_lval(jhinf, k);
	}
	lum_max = expf(lmax);
	if (lmax - lmin <= 4.f) {	/* image is LDR already */
		alph = -1.;
		gsca = 0.7/lum_max;
	} else {			/* compute Reinhard parameters */
		alph = 0.18*pow(4.0, (2.0*lavg-lmin-lmax)/(lmax-lmin));
		gsca = alph/exp(lavg);
	}
	tsca = 1.f;			/* TMO scaling to bring in gamut */
	for (j = 0; j < jhinf->cinfo.image_height && tsca > .42f; j++)
		for (i = 0; i < jhinf->cinfo.image_width; i++) {
			float		ttsca = tsca;
			boolean		changing;
			float		mul;
			if (skip_sample())
				continue;
			hdr_pval(pv, &jhinf->hdr, i, j);
			lv = hdr_lum(pv);
			if ((lv <= JH_LUM_MIN) | (lv > lum_max))
				continue;
			do {
				tm = lv*gsca;
				if (alph < 0)	/* linear TMO? */
					tm *= ttsca;
				else
					tm *= ttsca/(1.f + tm);
				if (tm <= .072f)
					break;
				mul = tm/lv;	/* gamut check */
				pv[0] *= mul;
				pv[1] *= mul;
				pv[2] *= mul;
				if ((changing = !jpeghdr_in_gamut(jhinf, pv)))
					ttsca *= .97f;
			} while (changing);
			if (ttsca >= .4f)	/* within acceptable range? */
				tsca = ttsca;
		}
	tmpp = jhinf->tmi;		/* tone-map image */
	for (j = 0; j < jhinf->cinfo.image_height; j++)
		for (i = 0; i < jhinf->cinfo.image_width; i++, tmpp++) {
			hdr_pval(pv, &jhinf->hdr, i, j);
			lv = hdr_lum(pv);
			if (lv <= JH_LUM_MIN) {
				*tmpp = 0;
				continue;
			}
			tm = lv*gsca;
			if (alph < 0)		/* linear TMO? */
				tm *= tsca;
			else
				tm *= tsca/(1.f + tm);

			*tmpp = (int)(256.f*powf(tm, invgam));
		}
}
