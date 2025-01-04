/*
 *  jhdecomp.c
 *
 *  Created by Greg Ward on 12/20/04.
 *  Copyright 2004 Sunnybrook Technologies, <www.sunnybrooktech.com>. 
 *  All rights reserved.
 *
 *  Routines for decompressing JPEG files, reconstructing high dynamic range
 *  from information in subband marker.
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

#define FEQ(a,b)	(((a) <= (b)+1e-4f) & ((a) >= (b)-1e-4f))

/* Initialize HDR JPEG compression object */
GLOBAL(void)
jpeghdr_CreateDecompress (jh_decompress_ptr jhinf, int jpvers, int jhvers,
				      size_t jpssiz, size_t jhssiz)
{
	/* clear struct in case of jump to jpeghdr_destroy_decompress */
	memset((void *)((char *)jhinf+jpssiz), 0, jhssiz-jpssiz);
	if (jhvers != JH_LIB_VERSION)
		ERREXIT2(&jhinf->cinfo, JERR_BAD_LIB_VERSION,
				JH_LIB_VERSION, jhvers);
	if (jhssiz != (int)sizeof(jpeghdr_decompress_struct))
		ERREXIT2(&jhinf->cinfo, JERR_BAD_STRUCT_SIZE,
			(int)sizeof(jpeghdr_decompress_struct),
				(int)jhssiz);
	jhinf->correction = JHprecorr;	/* reset by input file... */
	jhinf->alpha = 1.f;
	jhinf->beta = 1.f;
	jhinf->samp2nits = JH_SAMP2NITS_UNK;
	jpeg_CreateDecompress(&jhinf->cinfo, jpvers, jpssiz);
}

/* Free HDR decompress temp space */
LOCAL(void)
jh_free_decompress(jh_decompress_ptr jhinf)
{
	if (jhinf->tmi_ycc[0] != NULL) {
		free((void *)jhinf->tmi_ycc[0]);
		free((void *)jhinf->tmi_ycc[1]);
		free((void *)jhinf->tmi_ycc[2]);
	}
	jhinf->tmi_ycc[0] = jhinf->tmi_ycc[1] = jhinf->tmi_ycc[2] = NULL;
	if (jhinf->ycc_scan != NULL) {
		free((void *)jhinf->ycc_scan);
		jhinf->ycc_scan = NULL;
	}
	if (jhinf->sb_fp != NULL) {
		fclose(jhinf->sb_fp);
		jhinf->sb_fp = NULL;
	}
	if (jhinf->sbi != NULL) {
		free((void *)jhinf->sbi);
		jhinf->sbi = NULL;
	}
	jhinf->sbm_read = 0;
}

/* Finish HDR decompression cycle (cleans up temp buffers) */
GLOBAL(boolean)
jpeghdr_finish_decompress (jh_decompress_ptr jhinf)
{
	jh_free_decompress(jhinf);
	return jpeg_finish_decompress(&jhinf->cinfo);
}

/* Abort decompression cycle (free up memory resources) */
GLOBAL(void)
jpeghdr_abort_decompress (jh_decompress_ptr jhinf)
{
	jh_free_decompress(jhinf);
	jpeg_abort_decompress(&jhinf->cinfo);
}

/* Destruction of JPEG HDR compression objects */
GLOBAL(void)
jpeghdr_destroy_decompress (jh_decompress_ptr jhinf)
{
	jh_free_decompress(jhinf);
	jpeg_destroy_decompress(&jhinf->cinfo);
	memset((void *)&jhinf->cinfo, 0, sizeof(struct jpeg_decompress_struct));
}

/* Read a byte from JPEG data source */
LOCAL(int)
read_next_byte (j_decompress_ptr cinfo)
{
	if (cinfo->src->bytes_in_buffer == 0 &&
			!(*cinfo->src->fill_input_buffer)(cinfo))
		return EOF;
	cinfo->src->bytes_in_buffer--;
	return *cinfo->src->next_input_byte++;
}

/* Read specified number of bytes from JPEG data source */
LOCAL(size_t)
read_bytes (j_decompress_ptr cinfo, char *buf, size_t siz)
{
	int	nbread = 0;

	while (cinfo->src->bytes_in_buffer < siz) {
		const size_t	n = cinfo->src->bytes_in_buffer;
		if (n) {
			memcpy((void *)buf, 
				(const void *)cinfo->src->next_input_byte, n);
			buf += n;
			siz -= n;
			nbread += n;
			cinfo->src->next_input_byte += n;
			cinfo->src->bytes_in_buffer = 0;
		}
		if (!(*cinfo->src->fill_input_buffer)(cinfo))
			return nbread;
	}
	memcpy((void *)buf, (const void *)cinfo->src->next_input_byte, siz);
	nbread += siz;
	cinfo->src->next_input_byte += siz;
	cinfo->src->bytes_in_buffer -= siz;
	return nbread;
}

/* Skip (ignore) specified number of bytes from JPEG data source */
LOCAL(boolean)
skip_bytes (j_decompress_ptr cinfo, size_t len)
{
	(*cinfo->src->skip_input_data)(cinfo, len);
	return TRUE;
}

/* Load subband marker (WARNING: punts on suspended input) */
METHODDEF(boolean)
subband_parser (j_decompress_ptr cinfo)
{
	jh_decompress_ptr	jhinf = (jh_decompress_ptr)cinfo;
	char			buf[512];
	size_t			length, n;
	int			i;

	i = read_next_byte(cinfo) << 8;
	i |= read_next_byte(cinfo);
	if ((i < 0) | (i > 0xffff))
		return FALSE;
	length = i - 2;
	n = sizeof(buf) - 1;
	if ((int)n > length)
		n = length;
	if (read_bytes(cinfo, buf, n) != n)
		return FALSE;
	buf[n] = '\0';
	length -= n;
	if (jhinf->sb_fp == NULL) {		/* first subband marker */
		float	logmin, logmax;
		int	corr;
		if (strncmp(buf, JH_MARKER_MAGIC, JH_MARKER_MAGLN-1) ||
				sscanf(buf+JH_MARKER_MAGLN, JH_MARKER_FMT0,
				&jhinf->cvers,
				&logmin, &logmax,
				&jhinf->samp2nits,
				&jhinf->alpha, &jhinf->beta,
				&corr) != 7 ||
				jhinf->cvers/10 > JH_LIB_VERSION/10 ||
				logmin >= logmax)
			return skip_bytes(cinfo, length);
						/* prepare decode table */
		jhinf->correction = corr;
		for (i = 256; i--; )
			jhinf->sb_dec[i] = expf( logmin +
					(1.f/256.f)*(i+.5f)*(logmax-logmin) );
						/* open temporary file */
		jhinf->sb_fp = tmpfile();
		if (jhinf->sb_fp == NULL)
			ERREXIT(cinfo, JERR_FILE_WRITE);
	} else if (strncmp(buf, JH_MARKER_MAGIC, JH_MARKER_MAGLN-1) ||
			sscanf(buf+JH_MARKER_MAGLN, JH_MARKER_FMT1, &i) != 1 ||
			i != jhinf->sbm_read) {
		jhinf->sbm_read = 0;
		return skip_bytes(cinfo, length);
	}
						/* copy after marker header */
	i = strlen(buf+JH_MARKER_MAGLN) + JH_MARKER_MAGLN + 1;
	if (fwrite(buf+i, 1, n-i, jhinf->sb_fp) != n-i)
		ERREXIT(cinfo, JERR_FILE_WRITE);
	while (length > 0) {
		n = sizeof(buf);
		if (n > length)
			n = length;
		if (read_bytes(cinfo, buf, n) != n)
			return FALSE;
		length -= n;
		if (fwrite(buf, 1, n, jhinf->sb_fp) != n)
			ERREXIT(cinfo, JERR_FILE_WRITE);
	}
	jhinf->sbm_read++;			/* increment marker count */
	return TRUE;
}

/* Read JPEG header and determine if HDR subband is available */
GLOBAL(int)
jpeghdr_read_header (jh_decompress_ptr jhinf)
{
	int	rv;
						/* look for our marker */
	jpeg_set_marker_processor(&jhinf->cinfo, JH_APPM, subband_parser);
						/* read header */
	rv = jpeg_read_header(&jhinf->cinfo, TRUE);
						/* reset marker handling */
	jpeg_save_markers(&jhinf->cinfo, JH_APPM, 0);
						/* is it a standard JPEG? */
	if ((jhinf->sbm_read == 0) | (rv != JPEG_HEADER_OK)) {
		if (jhinf->sb_fp != NULL)
			fclose(jhinf->sb_fp);
		jhinf->sb_fp = NULL;
		return rv;			/* not an HDR JPEG */
	}
	return JPEG_HEADER_HDR;			/* else give the good news */
}

/* Decompress subband image from stream, assuming 8-bit JPEG encoding */
LOCAL(boolean)
decompress_subband (jh_decompress_ptr jhinf, int *xres, int *yres)
{
	struct jpeg_decompress_struct	cinfo;

	if (jhinf->sb_fp == NULL)
		return FALSE;
	rewind(jhinf->sb_fp);
	cinfo.err = jhinf->cinfo.err;
	jpeg_create_decompress(&cinfo);
	jpeg_stdio_src(&cinfo, jhinf->sb_fp);
	if (jpeg_read_header(&cinfo, TRUE) != JPEG_HEADER_OK)
		return FALSE;
	cinfo.out_color_space = JCS_GRAYSCALE;
	if (*xres <= cinfo.image_width/8)
		cinfo.scale_denom = 8;
	else if (*xres <= cinfo.image_width/4)
		cinfo.scale_denom = 4;
	else if (*xres <= cinfo.image_width/2)
		cinfo.scale_denom = 2;
	jpeg_start_decompress(&cinfo);
	jhinf->sbi = (UINT8 *)malloc(sizeof(UINT8) *
			cinfo.output_width * cinfo.output_height );
	if (jhinf->sbi == NULL)
		ERREXIT1(&jhinf->cinfo, JERR_OUT_OF_MEMORY, 210);
	*xres = cinfo.output_width;
	*yres = cinfo.output_height;
	while (cinfo.output_scanline < cinfo.output_height) {
		JSAMPROW	rptr[4];
		rptr[0] = jhinf->sbi + cinfo.output_scanline*cinfo.output_width;
		rptr[1] = rptr[0] + cinfo.output_width;
		rptr[2] = rptr[1] + cinfo.output_width;
		rptr[3] = rptr[2] + cinfo.output_width;
		if (!jpeg_read_scanlines(&cinfo, rptr, 4)) {
			jpeg_destroy_decompress(&cinfo);
			jh_free_decompress(jhinf);
			return FALSE;
		}
	}
	jpeg_destroy_decompress(&cinfo);
	fclose(jhinf->sb_fp);
	jhinf->sb_fp = NULL;
	return TRUE;
}

/* Preload LDR (tone-mapped) image */
LOCAL(boolean)
ldr_preload (jh_decompress_ptr jhinf)
{
	const int	out_width = jhinf->cinfo.output_width;
	const int	out_height = jhinf->cinfo.output_height;
	int		i;
	JSAMPLE		*samprow;
					/* preload tone-mapped image */
	for (i = 3; i--; ) {
		jhinf->tmi_ycc[i] = (JSAMPLE *)malloc(sizeof(JSAMPLE) *
							out_width * out_height);
		if (jhinf->tmi_ycc[i] == NULL)
			goto memerr;
	}
	samprow = (JSAMPLE *)malloc(sizeof(JSAMPLE)*3 * out_width);
	if (samprow == NULL)
		goto memerr;
	while (jhinf->cinfo.output_scanline < out_height) {
		JSAMPLE		*jsp;
		while (!jpeg_read_scanlines(&jhinf->cinfo, &samprow, 1)) {
			free((void *)samprow);
			return FALSE;	/* XXX what can we do?? */
		}
		jsp = samprow + out_width*3;
		i = jhinf->cinfo.output_scanline * out_width;
		while (jsp > samprow) {
			jhinf->tmi_ycc[2][--i] = *--jsp;
			jhinf->tmi_ycc[1][i] = *--jsp;
			jhinf->tmi_ycc[0][i] = *--jsp;
		}
	}
	free((void *)samprow);
	jpeg_finish_decompress(&jhinf->cinfo);
	return TRUE;
memerr:
	ERREXIT1(&jhinf->cinfo, JERR_OUT_OF_MEMORY, 213);
	return FALSE;	/* pro forma return */
}

/* Compute post-correction on ratio image */
LOCAL(boolean)
ri_postcorr (jh_decompress_ptr jhinf, int sb_width, int sb_height)
{
	const int	out_width = jhinf->cinfo.output_width;
	const int	out_height = jhinf->cinfo.output_height;
	UINT8		*sbi_corr;
	UINT8		*tmi_sml;
	UINT8		*tmi_rs;
	int		samp_rad;
	int		x, y, i, j;

	/******* Post-correct ratio image *******
	 * Note that all post-correction calculations are done in 8 bits.
	 * This is not exactly to spec, but since an 8-bit log ratio value
	 * is roughly equivalent an 8-bit gamma-compressed Y encoding
	 * (very roughly), it works well enough and saves us having to
	 * perform a lot of log-linear conversions along the way.
	 * Post-correction is a hack to begin with, after all.
	 */
	if (sizeof(JSAMPLE) != sizeof(UINT8))
		ERREXIT1(&jhinf->cinfo, JERR_BAD_PRECISION, sizeof(JSAMPLE)*8);
	tmi_sml = (UINT8 *)malloc(sizeof(UINT8) * sb_width * sb_height);
	tmi_rs = (UINT8 *)malloc(sizeof(UINT8) * out_width * out_height);
	if ((tmi_sml == NULL) | (tmi_rs == NULL))
		goto memerr;
	jpeghdr_resample8(jhinf->tmi_ycc[0], out_width, out_height,
			tmi_sml, sb_width, sb_height);
	jpeghdr_resample8(tmi_sml, sb_width, sb_height,
			tmi_rs, out_width, out_height);
	free((void *)tmi_sml);
	sbi_corr = (UINT8 *)malloc(sizeof(UINT8) * out_width * out_height);
	if (sbi_corr == NULL)
		goto memerr;
	samp_rad = out_width / sb_width;
						/* correct ratio image */
	memcpy((void *)sbi_corr, (const void *)jhinf->sbi,
			sizeof(UINT8)*out_width*out_height);
	for (y = out_height - samp_rad; y-- > samp_rad; )
	    for (x = out_width - samp_rad; x-- > samp_rad; ) {
		float	sigma, mult;
		int	sb_min=255, tm_min=255;
		int	sb_max=0, tm_max=0;
						/* synthesize high-freq. */
		if (jhinf->sbi[y*out_width + x] == 0)
			continue;
		for (j = y-samp_rad; j <= y+samp_rad; j += samp_rad)
		    for (i = x-samp_rad; i <= x+samp_rad; i += samp_rad) {
			const int	k = j*out_width + i;
			if (tmi_rs[k] < tm_min)
				tm_min = tmi_rs[k];
			if (tmi_rs[k] > tm_max)
				tm_max = tmi_rs[k];
			if (jhinf->sbi[k] < sb_min)
				sb_min = jhinf->sbi[k];
			if (jhinf->sbi[k] > sb_max)
				sb_max = jhinf->sbi[k];
		    }
		if (tm_max - tm_min <= 3)	/* XXX assumed JPEG epsilon */
			continue;
		i = y*out_width + x;
		sigma = (float)((sb_max - sb_min)*tmi_rs[i]) /
				((tm_max - tm_min)*jhinf->sbi[i]);
		mult = (float)jhinf->tmi_ycc[0][i] / tmi_rs[i];
		if (sigma < 1.f)
			mult = powf(mult, sigma);
		sb_max = (int)(mult*jhinf->sbi[i] + .5f);
		sbi_corr[i] = (sb_max > 255) ? 255 : sb_max;
	    }
	free((void *)tmi_rs);
	free((void *)jhinf->sbi);		/* save corrected ratio image */
	jhinf->sbi = sbi_corr;
	return TRUE;
memerr:
	ERREXIT1(&jhinf->cinfo, JERR_OUT_OF_MEMORY, 214);
	return FALSE;	/* pro forma return */
}

/* Start HDR decompression cycle, returning FALSE if input is suspended */
GLOBAL(boolean)
jpeghdr_start_decompress (jh_decompress_ptr jhinf)
{
	int		out_width, out_height;
	boolean		do_postcorrect;
	int		sb_width, sb_height;
	UINT8		*sbi_us;

	if (jhinf->sb_fp == NULL)
		ERREXIT1(&jhinf->cinfo, JERR_BAD_STATE,
				jhinf->cinfo.global_state);
	jhinf->output_scanline = 0;
	if ((jhinf->cinfo.out_color_space = jhinf->cinfo.jpeg_color_space)
			!= JCS_YCbCr)
		ERREXIT(&jhinf->cinfo, JERR_BAD_J_COLORSPACE);
	if (!jpeg_start_decompress(&jhinf->cinfo))
		return FALSE;
	out_width = jhinf->cinfo.output_width;
	out_height = jhinf->cinfo.output_height;
						/* decompress subband */
	sb_width = out_width; sb_height = out_height;
	if (!decompress_subband(jhinf, &sb_width, &sb_height))
		return FALSE;
	do_postcorrect = (jhinf->correction == JHpostcorr) &
				(sb_width < out_width);
	if (!do_postcorrect) {
		jhinf->ycc_scan = (JSAMPLE *)malloc(sizeof(JSAMPLE)*3 *
							out_width);
		if (jhinf->ycc_scan == NULL)
			goto memerr;
	}
	if ((sb_width == out_width) & (sb_height == out_height))
		return TRUE;			/* no correction needed */
						/* preload for post-correction */
	if (do_postcorrect && !ldr_preload(jhinf))
		return FALSE;
						/* upsample ratio image */
	sbi_us = (UINT8 *)malloc(sizeof(UINT8) * out_width * out_height);
	if (sbi_us == NULL)
		goto memerr;
	if ((jhinf->cvers == 10) & !do_postcorrect)
		jpeghdr_resample8orig(jhinf->sbi, sb_width, sb_height,
				sbi_us, out_width, out_height);
	else
		jpeghdr_resample8(jhinf->sbi, sb_width, sb_height,
				sbi_us, out_width, out_height);
	free((void *)jhinf->sbi);		/* save upsampled ratio image */
	jhinf->sbi = sbi_us;
	if (do_postcorrect)			/* post-correct ratio image? */
		return ri_postcorr(jhinf, sb_width, sb_height);
	return TRUE;
memerr:
	ERREXIT1(&jhinf->cinfo, JERR_OUT_OF_MEMORY, 215);
	return FALSE;	/* pro forma return */
}

/* Resaturate an RGB floating point value */
GLOBAL(void)
jpeghdr_resaturate (jh_decompress_ptr jhinf, JHSAMPLE pv[3])
{
	int		imin, i1, i2;
	JHSAMPLE	lum;
	float		satrat;

	if (FEQ(jhinf->alpha, 1.f) && FEQ(jhinf->beta, 1.f))
		return;
	lum = hdr_lum(pv);
	if (lum <= JH_LUM_MIN)
		return;
	if (FEQ(jhinf->beta, 1.f)) {
		satrat = 1.f/jhinf->alpha;
		lum *= 1.f - satrat;
		pv[0] = lum + satrat*pv[0];
		pv[1] = lum + satrat*pv[1];
		pv[2] = lum + satrat*pv[2];
		return;
	}
	imin = (pv[1] < pv[0]);
	if (pv[2] < pv[imin]) imin = 2;
	if (lum - pv[imin] <= 1e-4f)
		return;
	if ((i1 = imin+1) >= 3) i1 -= 3;
	if ((i2 = imin+2) >= 3) i2 -= 3;
	satrat = (lum - pv[imin])/(jhinf->alpha*lum);
	pv[imin] = lum*(1.f - powf(satrat, 1.f/jhinf->beta));
	satrat = powf(1.f - pv[imin]/lum, 1.f-jhinf->beta)/jhinf->alpha;
	lum *= 1.f - satrat;
	pv[i1] = lum + satrat*pv[i1];
	pv[i2] = lum + satrat*pv[i2];
}

/* Get tone-mapped YCC pixel from last decoded scanline */
GLOBAL(void)
jpeghdr_scan_ycc (jh_decompress_ptr jhinf, JSAMPLE ycc[3], int x)
{
	int	i;

	if (jhinf->output_scanline <= 0)
		return;
	if (jhinf->ycc_scan != NULL) {
		i = x*3;
		ycc[0] = jhinf->ycc_scan[i++];
		ycc[1] = jhinf->ycc_scan[i++];
		ycc[2] = jhinf->ycc_scan[i];
		return;
	}
	if (jhinf->tmi_ycc[0] == NULL)
		return;
	i = (jhinf->output_scanline - 1)*jhinf->cinfo.output_width + x;
	ycc[0] = jhinf->tmi_ycc[0][i];
	ycc[1] = jhinf->tmi_ycc[1][i];
	ycc[2] = jhinf->tmi_ycc[2][i];
}

/* Convert YCC to integer sRGB without clamping */
LOCAL(void)
ycc2srgb (const JSAMPLE ycc[3], int srgb[3])
{
	int	i;
	SHIFT_TEMPS

	i = 91881*(GETJSAMPLE(ycc[2])-128);
	srgb[0] = GETJSAMPLE(ycc[0]) + RIGHT_SHIFT(i,16);
	i = -22553*(GETJSAMPLE(ycc[1])-128) + -46818*(GETJSAMPLE(ycc[2])-128);
	srgb[1] = GETJSAMPLE(ycc[0]) + RIGHT_SHIFT(i,16);
	i = 116130*(GETJSAMPLE(ycc[1])-128);
	srgb[2] = GETJSAMPLE(ycc[0]) + RIGHT_SHIFT(i,16);
}

/* Get tone-mapped sRGB pixel from last decoded scanline */
GLOBAL(void)
jpeghdr_scan_rgb24 (jh_decompress_ptr jhinf, JSAMPLE rgb[3], int x)
{
	JSAMPLE	ycc[3];
	int	srgb[3];
	int	i;
						/* read as YCC */
	jpeghdr_scan_ycc(jhinf, ycc, x);
						/* convert & clamp */
	ycc2srgb(ycc, srgb);
	for (i = 3; i--; )
		rgb[i] = srgb[i] <= 0 ? 0 : srgb[i] >= 255 ? 255 : srgb[i];
}

/* Convert a 24-bit YCC pixel to a floating point RGB value */
GLOBAL(void)
jpeghdr_ycc2rgb (const JSAMPLE ycc[3], JHSAMPLE rgb[3])
{
#define below		(-230)
#define	above		440
	static float	gamma_lookup[above-below];
	int		srgb[3];
						/* convert to sRGB */
	ycc2srgb(ycc, srgb);
						/* gamma lookup */
	if (gamma_lookup[0] == 0) {
		int	i = 0;			/* initialize table */
		for ( ; i <= -11-below; i++)
			gamma_lookup[i] = -powf(1.f/1.055f *
					(-1.f/256.f*(below+i) + .055f), 2.4f);
		for ( ; i < 11-below; i++)
			gamma_lookup[i] = 1.f/256.f/12.92f * (below + i);
		for ( ; i < above-below; i++)
			gamma_lookup[i] = powf(1.f/1.055f *
					(1.f/256.f*(below+i) + .055f), 2.4f);
	}
	rgb[0] = gamma_lookup[srgb[0]-below];
	rgb[1] = gamma_lookup[srgb[1]-below];
	rgb[2] = gamma_lookup[srgb[2]-below];
#undef below
#undef above
}

/* Decompress the next HDR scanline on the input */
GLOBAL(boolean)
jpeghdr_read_scanline (jh_decompress_ptr jhinf, JHSAMPLE *sl)
{
	JSAMPLE	ycc[3];
	int	x, y;

	if (jhinf->output_scanline >= jhinf->cinfo.output_height ||
			jhinf->sbi == NULL ||
			(jhinf->ycc_scan != NULL && jhinf->output_scanline !=
						jhinf->cinfo.output_scanline))
		ERREXIT1(&jhinf->cinfo, JERR_BAD_STATE,
				jhinf->cinfo.global_state);
	 if (jhinf->ycc_scan != NULL &&
			!jpeg_read_scanlines(&jhinf->cinfo, &jhinf->ycc_scan, 1))
		return FALSE;
	y = jhinf->output_scanline++;
	for (x = jhinf->cinfo.output_width; x-- > 0; ) {
		const float	ratio = jhinf->sb_dec[
				 jhinf->sbi[y*jhinf->cinfo.output_width + x] ];
		JHSAMPLE	*rgb = sl + x*3;
		jpeghdr_scan_ycc(jhinf, ycc, x);
		jpeghdr_ycc2rgb(ycc, rgb);
		jpeghdr_resaturate(jhinf, rgb);
		rgb[0] *= ratio;
		rgb[1] *= ratio;
		rgb[2] *= ratio;
	}
	return TRUE;
}

/* Simplified HDR decompression interface */
GLOBAL(int)
jpeghdr_load_memory (const char *fname,
			JHSAMPLE **fimg, JSAMPLE **bimg, int xyres[2],
			float *samp2nits)
{
	FILE				*fp;
	int				hdr_ok;
	jpeghdr_decompress_struct       jhinf;
	struct jpeg_error_mgr		jerr;
	int				fscanlen=0, bscanlen=0;
	JSAMPLE				*bptr;

	if ((fimg == NULL) & (bimg == NULL)
			|| xyres == NULL || (xyres[0] < 0) | (xyres[1] < 0) ||
			(xyres[0] > 0) ^ (xyres[1] > 0))
		return -1;
	if (fname == NULL)
		fp = stdin;
	else if ((fp = fopen(fname, "rb")) == NULL)
		return -1;
	jhinf.cinfo.err = jpeg_std_error(&jerr);
	jpeghdr_create_decompress(&jhinf);
	jpeg_stdio_src(&jhinf.cinfo, fp);
	if (fimg == NULL)		/* LDR only? */
		hdr_ok = jpeg_read_header(&jhinf.cinfo, TRUE);
	else
		hdr_ok = jpeghdr_read_header(&jhinf);
	if (hdr_ok == JPEG_HEADER_OK)
		jhinf.cinfo.out_color_space = JCS_RGB;
	else if (hdr_ok != JPEG_HEADER_HDR)
		goto cleanup;
	if (xyres[0]) {			/* fit image to maximum dimensions */
		jhinf.cinfo.scale_denom = (jhinf.cinfo.image_width +
						xyres[0] - 1) / xyres[0];
		if (jhinf.cinfo.scale_denom < (jhinf.cinfo.image_height +
						xyres[1] - 1) / xyres[1])
			jhinf.cinfo.scale_denom = (jhinf.cinfo.image_height +
						xyres[1] - 1) / xyres[1];
		if (jhinf.cinfo.scale_denom > 8) {
			/* provided buffer too small? */
			if (fimg != NULL && *fimg != NULL) {
				hdr_ok = -1;
				goto cleanup;
			}
			if (bimg != NULL && *bimg != NULL) {
				hdr_ok = -1;
				goto cleanup;
			}
			jhinf.cinfo.scale_denom = 8;
		} else if (jhinf.cinfo.scale_denom > 4)
			jhinf.cinfo.scale_denom = 8;
		else if (jhinf.cinfo.scale_denom > 2)
			jhinf.cinfo.scale_denom = 4;
		else if (!jhinf.cinfo.scale_denom)
			jhinf.cinfo.scale_denom = 1;
	}
					/* set up LDR buffer (if any) */
	jpeg_calc_output_dimensions(&jhinf.cinfo);
	if (bimg != NULL) {
		if (*bimg == NULL) {
			*bimg = (JSAMPLE *)malloc( sizeof(JSAMPLE)*3 *
					jhinf.cinfo.output_width *
					jhinf.cinfo.output_height );
			if (*bimg == NULL)
				goto memerr;
			bscanlen = 3*jhinf.cinfo.output_width;
		} else if (!(bscanlen = 3*xyres[0])) {
			hdr_ok = -1;
			goto cleanup;
		} else if ((xyres[0] > jhinf.cinfo.output_width) |
				(xyres[1] > jhinf.cinfo.output_height)) {
			memset((void *)*bimg, 0, sizeof(JSAMPLE) *
						bscanlen * xyres[1]);
		}
	}
	if (hdr_ok == JPEG_HEADER_HDR) {
					/* set up HDR buffer */
		if (*fimg == NULL) {
			*fimg = (JHSAMPLE *)malloc( sizeof(JHSAMPLE)*3 *
					jhinf.cinfo.output_width *
					jhinf.cinfo.output_height );
			if (*fimg == NULL)
				goto memerr;
			fscanlen = 3*jhinf.cinfo.output_width;
		} else if (!(fscanlen = 3*xyres[0])) {
			hdr_ok = -1;
			goto cleanup;
		} else if ((xyres[0] > jhinf.cinfo.output_width) |
				(xyres[1] > jhinf.cinfo.output_height)) {
			memset((void *)*fimg, 0, sizeof(JHSAMPLE) *
						fscanlen * xyres[1]);
		}
		jpeghdr_start_decompress(&jhinf);
		while (jhinf.output_scanline < jhinf.cinfo.output_height) {
			if (bimg != NULL)
				bptr = *bimg + bscanlen*jhinf.output_scanline;
			jpeghdr_read_scanline(&jhinf, *fimg +
						fscanlen*jhinf.output_scanline);
			if (bimg != NULL) {
				int	x;
				for (x = 0; x < jhinf.cinfo.output_width;
								x++, bptr += 3)
					jpeghdr_scan_rgb24(&jhinf, bptr, x);
			}
		}
	} else /* hdr_ok == JPEG_HEADER_OK */ {
		if (bimg == NULL) {
			hdr_ok = -1;
			goto cleanup;
		}
		jpeg_start_decompress(&jhinf.cinfo);
		while (jhinf.cinfo.output_scanline < jhinf.cinfo.output_height) {
			bptr = *bimg + bscanlen*jhinf.cinfo.output_scanline;
			jpeg_read_scanlines(&jhinf.cinfo, &bptr, 1);
		}
	}
	xyres[0] = jhinf.cinfo.output_width;
	xyres[1] = jhinf.cinfo.output_height;
	if (samp2nits != NULL)
		*samp2nits = jhinf.samp2nits;
cleanup:
	jpeghdr_destroy_decompress(&jhinf);
	if (fp != stdin)
		fclose(fp);
	return hdr_ok;
memerr:
	ERREXIT1(&jhinf.cinfo, JERR_OUT_OF_MEMORY, 216);
	return -1;	/* pro forma return */
}
