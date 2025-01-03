/*
 * jpeghdr.h
 *
 * Copyright (C) 2005 Sunnybrook Technologies <www.sunnybrooktech.com>
 * All Rights Reserved.
 *
 * Extensions to Thomas Lane's JPEG library to encode high dynamic range
 * images using a backwards-compatible subband in application marker 11.
 * Where elegance and compatibility with libjpeg were at odds, we chose
 * compatibility in hopes that it would make code changes easier.
 *
 * A simplified interface for writing from and reading into memory is
 * provided by the jpeghdr_write_memory() and jpeghdr_load_memory()
 * routines.
 *
 * See <http://www.anyhere.com/gward/papers/cic05.pdf> for details, and
 * example library usage at the end of this header.
 *
 * Greg Ward, initial implementation January 2005.
 *
 */

#ifndef JPEGHDR_H
#define JPEGHDR_H

#include "jpeglib.h"

#define JH_LIB_VERSION		11		/* library version 1.1 */

#define JH_APPM			(JPEG_APP0+11)	/* HDR application marker # */
#define JH_MARKER_MAGIC		"HDR_RI "
#define JH_MARKER_MAGLN		7
#define JH_MARKER_FMT0 \
		"ver=%d\n ln0=%f, ln1=%f, s2n=%e alp=%f bet=%f cor=%d\n~"
#define JH_MARKER_FMT1		"ext=%d\n~"

#ifndef JH_LUM_MIN
#define JH_LUM_MIN		1e-10f		/* smallest luminance */
#endif
#ifndef JH_LUM_MAX
#define	JH_LUM_MAX		(1.f/JH_LUM_MIN)
#endif
#ifndef JH_HIST_SIZ
#define JH_HIST_SIZ		2048		/* default histogram size */
#endif
#ifndef JH_RNG_MIN
#define JH_RNG_MIN		2.5f		/* minimum ln subband range */
#endif
#ifndef JH_RNG_MAX
#define JH_RNG_MAX		16.f		/* maximum ln subband range */
#endif
#ifndef JH_PSMP_MAX
#define JH_PSMP_MAX		400000		/* maximum pixels to sample */
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef float   JHSAMPLE;			/* HDR scanline data type */

/* Subband Correction Method (full sampling, precorrection, postcorrection) */
typedef enum { JHfullsamp=0, JHprecorr=1, JHpostcorr=2 } JHCorrMethod;

#define JH_SAMP2NITS_UNK	0.0		/* Unknown sample-to-Nits */

/* Structure for retrieving in-core HDR source image data */
/* Pixels should be returned from getval as linear, CCIR 709 RGB values */
typedef struct hdr_image_s {
	const void *    im_base;		/* base of HDR image */
	int		h_stride;		/* bytes between pixels */
	int		v_stride;		/* bytes between scanlines */
						/* get pixel callback */
	JMETHOD(void, getval, (JHSAMPLE *rv, const struct hdr_image_s *im,
					int h, int v));
	void		*c_data;		/* client data pointer */
} hdr_image_struct;

typedef hdr_image_struct		*hdr_image_ptr;

/* Macro to point at pixel location (h,v) in image im */
#define hdr_iptr(im,h,v) \
	(const void *)( (const char *)(im)->im_base + \
			(h)*(im)->h_stride + (v)*(im)->v_stride )

/* Macro to get HDR pixel value from image im at location (h,v) */
#define hdr_pval(rv,im,h,v) \
	(*(im)->getval)(rv, im, h, v)

/* Macro to compute linear luminance from CCIR 709 RGB pixel */
#define hdr_lum(pv) \
	(0.2126f*(pv)[0] + 0.7152f*(pv)[1] + 0.0722f*(pv)[2])

/* Structure for encoding JPEG file with HDR subband */
typedef struct {
	struct jpeg_compress_struct
			cinfo;			/* JPEG compression struct */
	void		*c_data;		/* available to application */
	hdr_image_struct
			hdr;			/* HDR input image */
	long		*histo;			/* computed ln histogram */
	int		hist_len;		/* histogram length */
	float		hist_lminmax[2];	/* ln histogram bounds */
	unsigned long	histot;			/* total number of samples */
	UINT8		*tmi;			/* tone-mapped image data */
	float		gamma;			/* gamma of tone-mapped data */
	short		quality;		/* quality level (0-100) */
	JHCorrMethod    correction;		/* subband correction method */
	float		alpha, beta;		/* saturation parameters */
	float		samp2nits;		/* sample-to-nits multiplier */
} jpeghdr_compress_struct;

/* Macro to point at 8-bit gray pixel in tone-mapped image array */
#define hdr_tptr(jhinf,h,v) \
	((jhinf)->tmi + (jhinf)->cinfo.image_width*(v) + (h))

/* Macro to compute i'th histogram bin value (natural log) */
#define jpeghdr_histo_lval(jhinf,i) \
	((jhinf)->hist_lminmax[0] + \
		((i)+.5f)/(jhinf)->hist_len * ((jhinf)->hist_lminmax[1] - \
						(jhinf)->hist_lminmax[0]))

typedef jpeghdr_compress_struct		*jh_compress_ptr;

/* Structure for decoding JPEG file with HDR subband */
typedef struct {
	struct jpeg_decompress_struct
			cinfo;			/* JPEG decompression struct */
	int		cvers;			/* compressor version used */
	void		*c_data;		/* available to application */
	int		sbm_read;		/* subband markers read */
	FILE		*sb_fp;			/* subband image temp file */
	UINT8		*sbi;			/* upsampled subband image */
	float		sb_dec[256];		/* subband decoding table */
	JSAMPLE		*tmi_ycc[3];		/* tone-mapped YCbCr image */
	JSAMPLE		*ycc_scan;		/* YCbCr scanline holder */
	int		output_scanline;	/* current output scanline */
	JHCorrMethod    correction;		/* subband correction method */
	float		alpha, beta;		/* saturation parameters */
	float		samp2nits;		/* sample-to-nits multiplier */
} jpeghdr_decompress_struct;

typedef jpeghdr_decompress_struct	*jh_decompress_ptr;

/* Initialization of JPEG HDR compression objects.
 * jpeghdr_create_compress() and jpeghdr_create_decompress() are the
 * exported names that applications should call.  These expand to calls
 * on jpeghdr_CreateCompress and jpeghdr_CreateDecompress with additional
 * information passed for version mismatch checking.
 * NB: you must set up the error-manager BEFORE calling jpeg_hdr_create_xxx.
 */
#define jpeghdr_create_compress(jhinf) \
    jpeghdr_CreateCompress(jhinf, JPEG_LIB_VERSION, JH_LIB_VERSION, \
			sizeof(struct jpeg_compress_struct), \
			sizeof(jpeghdr_compress_struct))
#define jpeghdr_create_decompress(jhinf) \
    jpeghdr_CreateDecompress(jhinf, JPEG_LIB_VERSION, JH_LIB_VERSION, \
			sizeof(struct jpeg_decompress_struct), \
			sizeof(jpeghdr_decompress_struct))
EXTERN(void) jpeghdr_CreateCompress JPP((jh_compress_ptr jhinf,
				      int jpvers, int jhvers,
				      size_t jpssiz, size_t jhssiz));
EXTERN(void) jpeghdr_CreateDecompress JPP((jh_decompress_ptr jhinf,
					int jpvers, int jhvers,
					size_t jpssiz, size_t jhssiz));

/* Assign HDR source framebuffer */
EXTERN(void) jpeghdr_src_floatRGB JPP((jh_compress_ptr jhinf,
			const float *img, int width, int height));
EXTERN(void) jpeghdr_src_doubleRGB JPP((jh_compress_ptr jhinf,
			const double *img, int width, int height));

/* Compute natural log histogram */
EXTERN(void) jpeghdr_comp_histo JPP((jh_compress_ptr jhinf,
					int len, float minv, float maxv));

/* Allocate tone-mapped image holder */
EXTERN(void) jpeghdr_alloc_tmi JPP((jh_compress_ptr jhinf));

/* Is the given RGB color inside our YCC gamut after desaturation? */
EXTERN(boolean) jpeghdr_in_gamut JPP((jh_compress_ptr jhinf,
					const JHSAMPLE rgb[3]));

/* Free allocated tone-mapped image */
EXTERN(void) jpeghdr_free_tmi JPP((jh_compress_ptr jhinf));

/* Apply default tone-mapping operator (allocates tmi if NULL) */
EXTERN(void) jpeghdr_tonemap_default JPP((jh_compress_ptr jhinf));

/* Apply multiscale local tone-mapping operator (optional interface) */
EXTERN(void) jpeghdr_tonemap_multiscale JPP((jh_compress_ptr jhinf));

/* Perform HDR Subband JPEG compression */
EXTERN(void) jpeghdr_do_compress JPP((jh_compress_ptr jhinf));

/* Destruction of JPEG HDR compression object */
EXTERN(void) jpeghdr_destroy_compress JPP((jh_compress_ptr jhinf));

/* Simplified HDR compression interface */
EXTERN(long) jpeghdr_write_memory JPP((const char *fname,
			const JHSAMPLE *img, int width, int height,
			short qual, JHCorrMethod corr, float alpha, float beta,
			float samp2nits));

/* Read JPEG header and determine if HDR subband is available */
EXTERN(int) jpeghdr_read_header JPP((jh_decompress_ptr jhinf));
/* Return value is one of: */
     /* JPEG_SUSPENDED             Suspended due to lack of input data */
     /* JPEG_HEADER_OK             Found LDR image datastream */
     /* JPEG_HEADER_TABLES_ONLY    Found valid table-specs-only datastream */
#define JPEG_HEADER_HDR 3       /* Found valid HDR image */

/* Is jpeghdr_decompress_struct loading an HDR image? */
#define jpeghdr_decompress_ishdr(jh) \
		((jh)->sb_fp != NULL || (jh)->sbi != NULL)

/* Start HDR decompression cycle, returning FALSE if input is suspended */
EXTERN(boolean) jpeghdr_start_decompress JPP((jh_decompress_ptr jhinf));

/* Decompress the next HDR scanline on the input */
EXTERN(boolean) jpeghdr_read_scanline JPP((jh_decompress_ptr jhinf, JHSAMPLE *sl));

/* Abort decompression cycle */
EXTERN(void) jpeghdr_abort_decompress JPP((jh_decompress_ptr jhinf));

/* Finish HDR decompression cycle */
EXTERN(boolean) jpeghdr_finish_decompress JPP((jh_decompress_ptr jhinf));

/* Destruction of JPEG HDR decompression object */
EXTERN(void) jpeghdr_destroy_decompress JPP((jh_decompress_ptr jhinf));

/* Simplified HDR decompression interface */
EXTERN(int) jpeghdr_load_memory JPP((const char *fname,
			JHSAMPLE **fimg, JSAMPLE **bimg, int xyres[2],
			float *samp2nits));

/****************** Compression Example ********************
 * Because we need to make as many as 7 passes over the HDR
 * input data, the library requires that the calling application
 * keep this data in memory and provides a reasonably
 * general method for accessing it via the hdr_image_struct.
 * If the data is stored in a contiguous float or double RGB array,
 * the jpeghdr_src_floatRGB() and jpeghdr_src_doubleRGB() calls set
 * up access directly.  Alternate scanline ordering is handled easily
 * by changing the h_stride and v_stride members of hdr_image_struct.
 * If the data is stored as 16-bit data or is not in CCIR 709 RGB
 * colorspace, the application must create its own accessor callback
 * to perform translation.

	FILE				*fp = fopen("output.jpg", "wb");
	jpeghdr_compress_struct		jhinf;
	struct jpeg_error_mgr		jerr;

	jhinf.cinfo.err = jpeg_std_error(&jerr);
	// Reassign error handling functions as desired
	jpeghdr_create_compress(&jhinf);
	jpeg_stdio_dest(&jhinf.cinfo, fp);
	jpeghdr_src_floatRGB(&jhinf, my_img, my_img_width, my_img_height);
	// Override jhinf.hdr.h_stride and jhinf.hdr.v_stride if non-standard.
	// Assign jhinf.samp2nits appropriately.
	// Reset jhinf.quality and jhinf.correction if desired.
	// Assign jhinf.alpha and jhinf.beta to modify saturation.
	// Change other JPEG encoding defaults as desired.
	// Assign gamma for tone-mapping if non-standard.
	jpeghdr_tonemap_default(&jhinf);
	// Or, assign jhinf.tmi 8-bit grayscale values in scanline order
	jpeghdr_do_compress(&jhinf);
	jpeghdr_destroy_compress(&jhinf);
	fclose(fp);
	
************************************************************/

/****************** Simplified Compression *****************
 * In the case where the output image is being written from
 * a floating-point RGB frame buffer using a CCIR 709 color
 * space, and no special control is desired, the simplified
 * interface may be used.  The jpeghdr_write_memory() routine
 * takes the target file name, a pointer to the floating-
 * point data, and the size of the image, and returns the
 * number of bytes written to the output, or zero if the
 * file could not be open or a write error was encountered.
 * Scanlines are assumed contiguous in the frame buffer,
 * ordered from top to bottom (RGB pixels left to right).
 * Output will go to stdout if the first parameter is NULL.
 
	float	f_vga_buf[480][640][3];
	long	nbytes_written;

	// fill f_vga_buf with CCIR 709 RGB values
	nbytes_written = jpeghdr_write_memory("output.jpg",
				(const float *)f_vga_buf, 640, 480,
				90, JHprecorr, 1.f, 1.f,
				JH_SAMP2NITS_UNK);

************************************************************/

/**************** Decompression Example ********************
 * Because we may not know in advance whether a given JPEG
 * file contains HDR data or not, the caller checks the
 * return value of jpeghdr_read_header() to decide whether
 * to read it as a standard JPEG or an HDR JPEG.  In the
 * latter case, the tone-mapped 24-bit YCC scanline is also
 * provided by the library in a separate buffer, which
 * may be used for rapid display of the image if desired.
 * Note that you may choose to read image in LDR mode with
 * jpeg_read_scanlines() even after calling jpeghdr_read_header(),
 * but you must clean up afterwards with a call to
 * jpeghdr_finish_decompress() or jpeghdr_destroy_decompress(),
 * otherwise some resources will not be properly freed.  Also,
 * it is unsafe to call any standard JPEG library routines
 * directly after calling jpeghdr_start_decompress(), since
 * jpeghdr_start_decompress() may preload the entire file
 * under some circumstances (JHpostcorr).

	FILE				*fp = fopen("input.jpg", "rb");
	JHSAMPLE			*hdrscan;
	JSAMPLE				*ldrscan;
	jpeghdr_decompress_struct       jhinf;
	struct jpeg_error_mgr		jerr;

	jhinf.cinfo.err = jpeg_std_error(&jerr);
	// Reassign error handling functions as desired
	jpeghdr_create_decompress(&jhinf);
	jpeg_stdio_src(&jhinf.cinfo, fp);
	switch (jpeghdr_read_header(&jhinf)) {
	case JPEG_HEADER_OK:		// LDR image
		jpeg_start_decompress(&jhinf.cinfo);
		ldrscan = (JSAMPLE *)malloc(jhinf.cinfo.output_width *
						sizeof(JSAMPLE)*3)
		while (jhinf.cinfo.output_scanline < jhinf.cinfo.output_height) {
			jpeg_read_scanlines(&jhinf.cinfo, &ldrscan, 1);
			// do something with returned LDR scanline
		}
		free((void *)ldrscan);
		break;
	case JPEG_HEADER_HDR:		// HDR image
		jpeghdr_start_decompress(&jhinf);
		hdrscan = (JHSAMPLE *)malloc(jhinf.cinfo.output_width *
						sizeof(JHSAMPLE)*3)
		// Important: test jhinf.output_scanline, not jhinf.cinfo
		while (jhinf.output_scanline < jhinf.cinfo.output_height) {
			jpeghdr_read_scanline(&jhinf, hdrscan);
			// do something with returned HDR scanline
		}
		free((void *)hdrscan);
		break;
	case JPEG_SUSPENDED:
	case JPEG_HEADER_TABLES_ONLY:
	default:			// unhandled return value
		// report error and exit or return
	}
	// jpeghdr_finish_decompress(&jhinf);
	jpeghdr_destroy_decompress(&jhinf);
	fclose(fp);
	
************************************************************/

/**************** Simplified Decompression *****************
 * If the application does not need such tight control
 * over the decompression cycle, the simplified interface
 * jpeghdr_load_memory() may be called instead.  This
 * routine loads a JPEG image into memory that may be
 * preallocated, or the routine will allocate memory as
 * needed.  In either case, a maximum width and height
 * may be specified to restrict memory use.  The routine
 * returns the value returned by jpeghdr_read_header(),
 * or -1 if the file cannot be opened or the provided
 * memory buffer is insufficient.  The simplest form of
 * the call is given in the example below.
 
	int	xyres[2];
	float	*rgb_buf;
	float	samp2nits;
	int	header_ok;

	rgb_buf = NULL;			// calls malloc()
	xyres[0] = xyres[1] = 0;	// load any size
	header_ok = jpeghdr_load_memory("input.jpg",
				&rgb_buf, NULL, xyres,
				&samp2nits);
	// do something with returned buffer
	free((void *)rgb_buf);

 ************************ options *************************
 * If the caller wishes to restrict the loaded image
 * resolution, the xyres[] array above may be set to the
 * maximum width and height desired, and the image will
 * be downsampled as necessary so as not to exceed these
 * dimensions.  If a non-NULL pointer is passed in rgb_buf,
 * its dimensions must be given in xyres[], and the image
 * will be downsampled as necessary to avoid overflowing
 * the caller's frame buffer.  Furthermore, the image will
 * be placed in the upper-left corner of this frame and
 * zero-filled.  The returned dimensions in xyres[] give
 * the actual area of the loaded image in each case, but
 * the buffer size is guaranteed to match these dimensions
 * only if jpeghdr_load_memory() is passed a NULL pointer
 * in rgb_buf as in the above example.  In this example,
 * the library allocates the buffer to the actual image
 * size using malloc(), and the caller is responsible
 * for freeing the returned memory.  Lastly, it is
 * possible to use jpeghdr_load_memory() to simultaneously
 * load a tone-mapped sRGB version of the image by passing
 * a non-NULL pointer for the third parameter.  The
 * allocation and sizing of this sRGB buffer follows the
 * same logic just described.  Conversely, the second
 * argument may be passed NULL if the caller needs a simple
 * interface to read a standard JPEG image.  For example,
 * if a first call to jpeghdr_load_memory() with a non-NULL
 * second argument and NULL third argument returns -1,
 * indicating the JPEG is not HDR, then the application may
 * repeat the call with a non-NULL third argument to load
 * the standard LDR JPEG.  Input will be read from stdin
 * if the first argument is NULL.
************************************************************/

/* Internal routines that may be useful to applications as well */

/* Desaturate an RGB floating point value */
EXTERN(void) jpeghdr_desaturate JPP((jh_compress_ptr jhinf, JHSAMPLE pv[3]));

/* Convert a floating point RGB value to a 24-bit YCC pixel */
EXTERN(void) jpeghdr_rgb2ycc JPP((JHSAMPLE rgb[3], JSAMPLE ycc[3]));

/* Resample an 8-bit planar image */
EXTERN(void) jpeghdr_resample8 JPP((const UINT8 *simg, int swidth, int sheight,
				UINT8 *dimg, int dwidth, int dheight));

/* Original (flawed) upsampling algorithm */
EXTERN(void) jpeghdr_resample8orig JPP((const UINT8 *simg, int swidth, int sheight,
				UINT8 *dimg, int dwidth, int dheight));

/* Get tone-mapped sRGB pixel from last decoded scanline */
EXTERN(void) jpeghdr_scan_rgb24 JPP((jh_decompress_ptr jhinf,
				JSAMPLE rgb[3], int x));

/* Get tone-mapped YCC pixel from last decoded scanline */
EXTERN(void) jpeghdr_scan_ycc JPP((jh_decompress_ptr jhinf,
				JSAMPLE ycc[3], int x));

/* Convert a 24-bit YCC pixel to a floating point RGB value */
EXTERN(void) jpeghdr_ycc2rgb JPP((const JSAMPLE ycc[3], JHSAMPLE rgb[3]));

/* Resaturate an RGB floating point value */
EXTERN(void) jpeghdr_resaturate JPP((jh_decompress_ptr jhinf, JHSAMPLE pv[3]));

#ifdef __cplusplus
}
#endif

#endif  /* ! JPEGHDR_H */
