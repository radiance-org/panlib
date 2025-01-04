/*
 *  jhcomp.c
 *
 *  Created by Greg Ward on 10/5/04.
 *  Copyright 2004 Sunnybrook Technologies, <www.sunnybrooktech.com>. 
 *  All rights reserved.
 *
 *  Routines for compressing JPEG files, maintaining high dynamic range
 *  information in subband marker.
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

#define rjitter()	(.5f - (1.f/RAND_MAX)*rand())

#define FEQ(a,b)	(((a) <= (b)+1e-4f) & ((a) >= (b)-1e-4f))

#ifndef TM_MIN
#define TM_MIN		20		/* minimum allowed tone-mapped value */
#endif

/* Initialize HDR JPEG compression object */
GLOBAL(void)
jpeghdr_CreateCompress (jh_compress_ptr jhinf, int jpvers, int jhvers,
				      size_t jpssiz, size_t jhssiz)
{
	/* clear struct in case of jump to jpeghdr_destroy_compress */
	memset((void *)((char *)jhinf+jpssiz), 0, jhssiz-jpssiz);
	if (jhvers != JH_LIB_VERSION)
		ERREXIT2(&jhinf->cinfo, JERR_BAD_LIB_VERSION,
				JH_LIB_VERSION, jhvers);
	if (jhssiz != (int)sizeof(jpeghdr_compress_struct))
		ERREXIT2(&jhinf->cinfo, JERR_BAD_STRUCT_SIZE,
			(int)sizeof(jpeghdr_compress_struct),
				(int)jhssiz);
	jhinf->gamma = 2.2f;
	jhinf->quality = 90;
	jhinf->correction = JHprecorr;
	jhinf->alpha = 1.f;
	jhinf->beta = 1.f;
	jhinf->samp2nits = JH_SAMP2NITS_UNK;
	jpeg_CreateCompress(&jhinf->cinfo, jpvers, jpssiz);
}

/* Destruction of JPEG HDR compression objects */
GLOBAL(void)
jpeghdr_destroy_compress (jh_compress_ptr jhinf)
{
	memset((void *)&jhinf->hdr, 0, sizeof(hdr_image_struct));
	jpeghdr_comp_histo(jhinf, 0, 0, 0);
	jpeghdr_free_tmi(jhinf);
	jpeg_destroy_compress(&jhinf->cinfo);
	memset((void *)&jhinf->cinfo, 0, sizeof(struct jpeg_compress_struct));
}

/* Callback to transfer 32-bit floating-point RGB pixels */
METHODDEF(void)
cvt_floatRGB (JHSAMPLE *rv, const hdr_image_struct *im, int h, int v)
{
	const float	*fp = (const float *)hdr_iptr(im, h, v);
	rv[0] = fp[0]; rv[1] = fp[1]; rv[2] = fp[2];
}

/* Assign 32-bit IEEE float HDR source framebuffer */
GLOBAL(void)
jpeghdr_src_floatRGB (jh_compress_ptr jhinf, const float *img,
		int width, int height)
{
	if (img == NULL || width <= 0 || height <= 0)
		ERREXIT(&jhinf->cinfo, JERR_EMPTY_IMAGE);
	jhinf->cinfo.image_width = width;
	jhinf->cinfo.image_height = height;
	jhinf->hdr.im_base = (const void *)img;
	jhinf->hdr.h_stride = sizeof(float) * 3;
	jhinf->hdr.v_stride = jhinf->hdr.h_stride * width;
	jhinf->hdr.getval = &cvt_floatRGB;
	jhinf->hdr.c_data = NULL;
}
					
/* Callback to transfer 64-bit floating-point RGB pixels */
METHODDEF(void)
cvt_doubleRGB (JHSAMPLE *rv, const hdr_image_struct *im, int h, int v)
{
	const double	*dp = (const double *)hdr_iptr(im, h, v);
	rv[0] = dp[0]; rv[1] = dp[1]; rv[2] = dp[2];
}

/* Assign 64-bit IEEE double HDR source framebuffer */
GLOBAL(void)
jpeghdr_src_doubleRGB (jh_compress_ptr jhinf, const double *img,
		int width, int height)
{
	if (img == NULL || width <= 0 || height <= 0)
		ERREXIT(&jhinf->cinfo, JERR_EMPTY_IMAGE);
	jhinf->cinfo.image_width = width;
	jhinf->cinfo.image_height = height;
	jhinf->hdr.im_base = (const void *)img;
	jhinf->hdr.h_stride = sizeof(double) * 3;
	jhinf->hdr.v_stride = jhinf->hdr.h_stride * width;
	jhinf->hdr.getval = &cvt_doubleRGB;
	jhinf->hdr.c_data = NULL;
}

/* Check tone-mapping and fix any problematic values */
LOCAL(boolean)
check_tonemap (jh_compress_ptr jhinf)
{
	boolean		tmOK = FALSE;
	UINT8		*tmpp = jhinf->tmi;
	long		cnt = (long)jhinf->cinfo.image_width *
					jhinf->cinfo.image_height;
	
	if (tmpp == NULL)
		return FALSE;
	for ( ; cnt-- > 0; tmpp++)
		if (*tmpp < TM_MIN)
			*tmpp = TM_MIN;
		else if (*tmpp > TM_MIN)
			tmOK = TRUE;
	return tmOK;
}

/* Compute proper range of log luminance values for subband */
LOCAL(void)
log_range (float loglim[2], jh_compress_ptr jhinf, const float gamma_recip[256])
{
	float		rmin = 1.f;
	float		rmax = 1.f;
	float		lmin, lmax;
	JHSAMPLE	pv[3];
	float		ratio;
	UINT8		*tmpp;
	float		v;
	int		i, j;
	long		cnt;
					/* find ratio value range */
	if (!jhinf->hist_len)
		jpeghdr_comp_histo(jhinf, JH_HIST_SIZ, JH_LUM_MIN, JH_LUM_MAX);
					/* drop top & bottom 0.02% of values */
	cnt = jhinf->histot / 5000;
	if (!cnt) cnt = 1;
	for (j = 0; j < jhinf->hist_len; j++)
		if ((cnt -= jhinf->histo[j]) <= 0)
			break;
	lmin = expf(jpeghdr_histo_lval(jhinf, j));
	cnt = jhinf->histot / 5000;
	if (!cnt) cnt = 1;
	for (j = jhinf->hist_len; j--; )
		if ((cnt -= jhinf->histo[j]) <= 0)
			break;
	lmax = expf(jpeghdr_histo_lval(jhinf, j));
	tmpp = jhinf->tmi;
	for (j = 0; j < jhinf->cinfo.image_height; j++)
		for (i = 0; i < jhinf->cinfo.image_width; i++, tmpp++) {
			hdr_pval(pv, &jhinf->hdr, i, j);
			v = hdr_lum(pv);
			if ((lmin > v) | (v > lmax))
				continue;
			ratio = v * gamma_recip[*tmpp];
			if (ratio < rmin)
				rmin = ratio;
			if (ratio > rmax)
				rmax = ratio;
		}
	loglim[0] = logf(rmin);
	loglim[1] = logf(rmax);		/* check subband range */
	if (loglim[1] - loglim[0] < JH_RNG_MIN) {
		ratio = .5f*(JH_RNG_MIN - loglim[1] + loglim[0]);
		loglim[0] -= ratio;
		loglim[1] += ratio;
	} else if (loglim[1] - loglim[0] > JH_RNG_MAX) {
		loglim[0] = loglim[1] - JH_RNG_MAX;
	}
}

/* Compute ratio image (& take care of zeroes on input) */
LOCAL(UINT8 *)
compute_subband (jh_compress_ptr jhinf, float loglim[2],
			const float gamma_lookup[256])
{
	float		gamma_recip[256];
	float		logscale;
	UINT8		*sbi_fullres;
	JHSAMPLE	pv[3];
	float		ratio;
	UINT8		*tmpp;
	UINT8		*sbpp;
	float		v;
	int		i, j;
					/* reciprocate gamma table */
	for (i = 256; i--; )
		gamma_recip[i] = 1.f/gamma_lookup[i];
					/* find ratio value range */
	log_range(loglim, jhinf, gamma_recip);
	logscale = 256.f/(loglim[1] - loglim[0]);

					/* allocate subband image */
	sbi_fullres = (UINT8 *)malloc(sizeof(UINT8) *
			jhinf->cinfo.image_width * jhinf->cinfo.image_height);
	if (sbi_fullres == NULL)
		ERREXIT1(&jhinf->cinfo, JERR_OUT_OF_MEMORY, 202);
	tmpp = jhinf->tmi;
	sbpp = sbi_fullres;		/* compute subband */
	for (j = 0; j < jhinf->cinfo.image_height; j++)
		for (i = 0; i < jhinf->cinfo.image_width; i++) {
			hdr_pval(pv, &jhinf->hdr, i, j);
			ratio = hdr_lum(pv);
			if (ratio <= JH_LUM_MIN) {
				*tmpp++ = 0;
				*sbpp++ = 0;
				continue;
			}
			ratio *= gamma_recip[*tmpp++];
			v = logscale*(logf(ratio) - loglim[0]);
			*sbpp++ = (v < 0.5f) ? 0 : (v >= 255.5f) ? 255 :
					(int)(v + rjitter()); 
		}
	return sbi_fullres;
}

/* Resample subband image (& precorrect tone-mapped image) */
LOCAL(UINT8 *)
resample_subband (jh_compress_ptr jhinf, UINT8 *sbi_fullres,
			int sbi_width, int sbi_height, const float loglim[2])
{
	const float	invgam = 1.f/jhinf->gamma;
	const float	logscale = 256.f/(loglim[1] - loglim[0]);
	float		ri_lookup[256];
	UINT8		*sbi_final;
	JHSAMPLE	pv[3];
	float		tmtarget;
	UINT8		*tmpp;
	UINT8		*sbpp;
	float		v;
	int		i, j;
					/* allocate reduced subband */
	sbi_final = (UINT8 *)malloc(sizeof(UINT8) * sbi_width * sbi_height);
	if (sbi_final == NULL)
		ERREXIT1(&jhinf->cinfo, JERR_OUT_OF_MEMORY, 203);
	jpeghdr_resample8(sbi_fullres, jhinf->cinfo.image_width,
				jhinf->cinfo.image_height,
				sbi_final, sbi_width, sbi_height);
	if (jhinf->correction != JHprecorr)
		return sbi_final;
					/* precorrect tmi */
	jpeghdr_resample8(sbi_final, sbi_width, sbi_height,
			sbi_fullres, jhinf->cinfo.image_width,
			jhinf->cinfo.image_height);
	for (i = 256; i--; )
		ri_lookup[i] = expf(-((i+.5f)/logscale + loglim[0]));
	tmpp = jhinf->tmi;
	sbpp = sbi_fullres;
	for (j = 0; j < jhinf->cinfo.image_height; j++)
		for (i = 0; i < jhinf->cinfo.image_width;
					i++, tmpp++, sbpp++) {
			if (!*tmpp)
				continue;
			hdr_pval(pv, &jhinf->hdr, i, j);
			tmtarget = hdr_lum(pv)*ri_lookup[*sbpp];
			v = 256.f*powf(tmtarget, invgam);
			*tmpp = (v <= TM_MIN) ? TM_MIN :
					(v >= 255.f) ? 255 :
					(int)v;
		}
	return sbi_final;
}

/* Compress subband image to temporary file using JPEG encoding */
LOCAL(FILE *)
compress_subband (UINT8 *sbi, int width, int height,
			int qual, struct jpeg_error_mgr *errh)
{
	FILE				*fp = tmpfile();
	struct jpeg_compress_struct	cinfo;

	if (fp == NULL)
		return NULL;
	cinfo.err = errh;
	jpeg_create_compress(&cinfo);
	jpeg_stdio_dest(&cinfo, fp);
	cinfo.image_width = width;
	cinfo.image_height = height;
	cinfo.input_components = 1;
	cinfo.in_color_space = JCS_GRAYSCALE;
	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, qual, TRUE);
	jpeg_start_compress(&cinfo, TRUE);
	while (cinfo.next_scanline < cinfo.image_height) {
		JSAMPROW	rptr = (JSAMPLE *)sbi +
						cinfo.next_scanline*width;
		jpeg_write_scanlines(&cinfo, &rptr, 1);
	}
	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);
	return fp;				/* do not rewind */
}

/* Write out final tone-mapped image */
LOCAL(void)
ldr_write (jh_compress_ptr jhinf, const float gamma_lookup[256])
{
	JSAMPLE		*ycc_scanline;
	JHSAMPLE	pv[3];
	float		ratio;
	UINT8		*tmpp;
	int		i, j;
					/* allocate scanline */
	ycc_scanline = (JSAMPLE *)malloc(sizeof(JSAMPLE)*3 *
					jhinf->cinfo.image_width);
	if (ycc_scanline == NULL)
		ERREXIT1(&jhinf->cinfo, JERR_OUT_OF_MEMORY, 204);
	tmpp = jhinf->tmi;		/* write out LDR image */
	for (j = 0; j < jhinf->cinfo.image_height; j++) {
		for (i = 0; i < jhinf->cinfo.image_width; i++, tmpp++) {
			if (!*tmpp) {
				ycc_scanline[i*3] = 0;
				ycc_scanline[i*3 + 1] =
				ycc_scanline[i*3 + 2] = 128;
				continue;
			}
			hdr_pval(pv, &jhinf->hdr, i, j);
			ratio = gamma_lookup[*tmpp] / hdr_lum(pv);
			pv[0] *= ratio; pv[1] *= ratio; pv[2] *= ratio;
			jpeghdr_desaturate(jhinf, pv);
			jpeghdr_rgb2ycc(pv, ycc_scanline + i*3);
		}
		jpeg_write_scanlines(&jhinf->cinfo, &ycc_scanline, 1);
	}
	free((void *)ycc_scanline);	/* cleanup */
}

/* Desaturate an RGB floating point value */
GLOBAL(void)
jpeghdr_desaturate (jh_compress_ptr jhinf, JHSAMPLE pv[3])
{
	JHSAMPLE	lum;
	float		satrat;

	if (FEQ(jhinf->alpha, 1.f) && FEQ(jhinf->beta, 1.f))
		return;
	lum = hdr_lum(pv);
	if (lum <= JH_LUM_MIN)
		return;
	if (!FEQ(jhinf->beta, 1.f)) {
		JHSAMPLE	vmin = (pv[0] < pv[1]) ? pv[0] : pv[1];
		if (pv[2] < vmin) vmin = pv[2];
		if (lum - vmin <= 1e-4f)
			return;
		satrat = jhinf->alpha * powf(1.f - vmin/lum, jhinf->beta - 1.f);
	} else
		satrat = jhinf->alpha;
	lum *= 1.f - satrat;
	pv[0] = lum + satrat*pv[0];
	pv[1] = lum + satrat*pv[1];
	pv[2] = lum + satrat*pv[2];
}

/* Convert a floating point RGB value to a 24-bit YCC pixel */
GLOBAL(void)
jpeghdr_rgb2ycc (JHSAMPLE rgb[3], JSAMPLE ycc[3])
{
	float	rgbgam[3];
	float	yccgam[3];
	int	i;
					/* convert linear RGB to sRGB */
	for (i = 3; i--; )
		if (rgb[i] < -0.0031308f)
			rgbgam[i] = -1.055f*powf(-rgb[i], 1.f/2.4f) + .055f;
		else if (rgb[i] <= 0.0031308f)
			rgbgam[i] = 12.92f*rgb[i];
		else
			rgbgam[i] = 1.055f*powf(rgb[i], 1.f/2.4f) - .055f;

	yccgam[0] = 0.299f*rgbgam[0] + 0.587f*rgbgam[1] + 0.114f*rgbgam[2];
	yccgam[1] = -0.16874f*rgbgam[0] + -0.33126f*rgbgam[1] + 0.5f*rgbgam[2] + 0.5f;
	yccgam[2] = 0.5f*rgbgam[0] + -0.41869f*rgbgam[1] + -0.08131f*rgbgam[2] + 0.5f;
	ycc[0] = (yccgam[0] <= 0) ? 0 : (yccgam[0] >= 1.f) ? 255 :
						(int)(256.f*yccgam[0]);
	ycc[1] = (yccgam[1] <= 0) ? 0 : (yccgam[1] >= 1.f) ? 255 :
						(int)(256.f*yccgam[1]);
	ycc[2] = (yccgam[2] <= 0) ? 0 : (yccgam[2] >= 1.f) ? 255 :
						(int)(256.f*yccgam[2]);
}

/* Is the given RGB color going to fit inside the YCC gamut? */
GLOBAL(boolean)
jpeghdr_in_gamut (jh_compress_ptr jhinf, const JHSAMPLE rgb_t[3])
{
	JHSAMPLE	rgb[3];
	JSAMPLE		ycc[3];
	
	rgb[0] = rgb_t[0]; rgb[1] = rgb_t[1]; rgb[2] = rgb_t[2];
	jpeghdr_desaturate(jhinf, rgb);
	jpeghdr_rgb2ycc(rgb, ycc);
	return ( (0 < GETJSAMPLE(ycc[0])) & (GETJSAMPLE(ycc[0]) < 255) &&
		(0 < GETJSAMPLE(ycc[2])) & (GETJSAMPLE(ycc[2]) < 255) &&
		(0 < GETJSAMPLE(ycc[1])) & (GETJSAMPLE(ycc[1]) < 255) );
}

/* Perform HDR Subband JPEG compression */
GLOBAL(void)
jpeghdr_do_compress (jh_compress_ptr jhinf)
{
	float		loglim[2];
	float		gamma_lookup[256];
	UINT8		*sbi_fullres;
	UINT8		*sbi_final;
	int		sbi_width, sbi_height;
	FILE		*sbc_fp;
	long		sbc_length;
	float		sbi_scale;
	int		i, n, mi;
					/* fix tone-mapping */
	if (!check_tonemap(jhinf))	/* XXX fails for blank images! */
		ERREXIT(&jhinf->cinfo, JERR_CONVERSION_NOTIMPL);

					/* fill gamma lookup table */
	for (i = 256; i--; )
		gamma_lookup[i] = powf((1.f/256.f)*(i + .5f), jhinf->gamma);

					/* compute subband image */
	sbi_fullres = compute_subband(jhinf, loglim, gamma_lookup);

					/* check downsampling setting */
	if (jhinf->quality >= 95)
		jhinf->correction = JHfullsamp;
	if (jhinf->correction == JHfullsamp) {	/* no downsampling */
		sbi_width = jhinf->cinfo.image_width;
		sbi_height = jhinf->cinfo.image_height;
		sbi_final = sbi_fullres;
	} else {				/* downsample ratio image */
		/* XXX this is highly heuristic... */
		sbi_scale = (jhinf->quality - 25)/70.f *
				powf(512.f*512.f /
					(jhinf->cinfo.image_width *
						jhinf->cinfo.image_height),
					0.25f);
		if (sbi_scale > 0.8f)
			sbi_scale = 0.8f;
		else if (sbi_scale < 0.25f)
			sbi_scale = 0.25f;
		sbi_width = (int)(sbi_scale*jhinf->cinfo.image_width + 0.5f);
		sbi_height = (int)(sbi_scale*jhinf->cinfo.image_height + 0.5f);
						/* downsample subband */
		sbi_final = resample_subband(jhinf, sbi_fullres,
						sbi_width, sbi_height, loglim);
		free((void *)sbi_fullres);
	}
					/* compress log ratio image */
	i = (3*jhinf->quality +
			jhinf->quality*jhinf->cinfo.image_width/sbi_width +
			jhinf->quality*(int)(loglim[1]-loglim[0]+.5f)/5)/5;
	if (i > 100) i = 100;
	sbc_fp = compress_subband(sbi_final, sbi_width, sbi_height,
					i, jhinf->cinfo.err);
	if (sbc_fp == NULL)
		ERREXIT(&jhinf->cinfo, JERR_FILE_WRITE);
	sbc_length = ftell(sbc_fp);
	rewind(sbc_fp);
					/* start JPEG-HDR header */
	jhinf->cinfo.input_components = 3;
	jhinf->cinfo.in_color_space = JCS_YCbCr;
	jpeg_set_defaults(&jhinf->cinfo);
	jpeg_set_quality(&jhinf->cinfo, jhinf->quality, TRUE);
	jpeg_start_compress(&jhinf->cinfo, TRUE);
	mi = 0;				/* add compressed subband marker(s) */
	strcpy((char *)sbi_final, JH_MARKER_MAGIC);
	sprintf((char *)sbi_final + JH_MARKER_MAGLN, JH_MARKER_FMT0,
			JH_LIB_VERSION, loglim[0], loglim[1], jhinf->samp2nits,
			jhinf->alpha, jhinf->beta, jhinf->correction);
	while (sbc_length > 0) {
		i = strlen((char *)sbi_final + JH_MARKER_MAGLN) +
				JH_MARKER_MAGLN + 1;
		n = fread((void *)(sbi_final+i), 1, 0xFF00-i, sbc_fp);
		if (n <= 0)
			ERREXIT(&jhinf->cinfo, JERR_FILE_READ);
		jpeg_write_marker(&jhinf->cinfo, JH_APPM,
				(JOCTET *)sbi_final, i+n);
		sbc_length -= n;
		sprintf((char *)sbi_final + JH_MARKER_MAGLN,
				JH_MARKER_FMT1, ++mi);
	}
	fclose(sbc_fp);			/* done with ratio image */
	free((void *)sbi_final);
					/* write tone-mapped image */
	ldr_write(jhinf, gamma_lookup);
					/* clean up */
	jpeg_finish_compress(&jhinf->cinfo);
}

/* Simplified HDR compression interface */
GLOBAL(long)
jpeghdr_write_memory (const char *fname,
			const JHSAMPLE *img, int width, int height,
			short qual, JHCorrMethod corr, float alpha, float beta,
			float samp2nits)
{
	FILE				*fp;
	long				nbytes;
	jpeghdr_compress_struct		jhinf;
	struct jpeg_error_mgr		jerr;

	if (img == NULL || (width <= 0) | (height <= 0))
		return 0;
	if (fname == NULL)
		fp = stdout;
	else if ((fp = fopen(fname, "wb")) == NULL)
		return 0;
	jhinf.cinfo.err = jpeg_std_error(&jerr);
	jpeghdr_create_compress(&jhinf);
	jpeg_stdio_dest(&jhinf.cinfo, fp);
	jpeghdr_src_floatRGB(&jhinf, img, width, height);
	jhinf.samp2nits = samp2nits;
	jhinf.quality = qual;
	jhinf.correction = corr;
	jhinf.alpha = alpha; jhinf.beta = beta;
	jpeghdr_tonemap_default(&jhinf);
	jpeghdr_do_compress(&jhinf);
	jpeghdr_destroy_compress(&jhinf);
	nbytes = (fflush(fp) == 0) ? ftell(fp) : 0;
	if (fp != stdout)
		fclose(fp);
	return nbytes;
}
