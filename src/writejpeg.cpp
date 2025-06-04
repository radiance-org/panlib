/*
 *  writejpeg.cpp
 *  panlib
 *
 *  JPEG image writer
 *
 *  Created by gward on Fri Jun 08 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <iostream>
#include <fstream>
using namespace std;
#include "system.h"
#include "imgio.h"
#include "imgwriter.h"
#include "color.h"
#include "jstreamsrc.h"
#include "jpeghdr.h"

// Callback to convert 32-bit floating-point tristimulus values
static void
cvt_floatColor(JHSAMPLE *rv, const hdr_image_struct *im, int h, int v)
{
	const float 	*iv = (const float *)hdr_iptr(im, h, v);
	const float	(*mat)[3] = (const float (*)[3])im->c_data;
	
        rv[0] = mat[0][0]*iv[0] + mat[0][1]*iv[1] + mat[0][2]*iv[2];
        rv[1] = mat[1][0]*iv[0] + mat[1][1]*iv[1] + mat[1][2]*iv[2];
        rv[2] = mat[2][0]*iv[0] + mat[2][1]*iv[1] + mat[2][2]*iv[2];
}

// Write a JPEG image from floating-point RGB color space
static long
JPEGwriteImageHDR(const char *fname, const ImgWriteBuf *wb)
{
	COLORMAT	myMat;
	char		pval[256];
	int		n;
						// open file for output
	ofstream	ostr(fname, ios::out|ios::trunc|ios::binary);
	if (!ostr.is_open())
		return 0;
						// set up input & error handling
	jpeghdr_compress_struct		jhinf;
	struct jpeg_error_mgr		jerr;
	jumper_struct			myjumper;
	myjumper.w = JMSG_NOMESSAGE;
	jhinf.cinfo.err = jpeg_std_error(&jerr);
	jerr.output_message = jpeg_error_output;
	jerr.emit_message = jpeg_emit_message;
	jerr.error_exit = jpeg_error_jump;
	if (setjmp(myjumper.b)) {		// JPEG error jumps to here
		jpeghdr_destroy_compress(&jhinf);
		return 0;			// all clean -- return failure
	}
	jhinf.cinfo.client_data = (void *)&myjumper;
	jpeghdr_create_compress(&jhinf);
	jpeg_stream_dest(&jhinf.cinfo, &ostr);
						// assign HDR source image
	if (!PmatchColorSpace(wb->csp, &ICS_RGB709, PICMptype|PICMprims)) {
		jhinf.cinfo.image_width = wb->xres;
		jhinf.cinfo.image_height = wb->yres;
		jhinf.hdr.im_base = (const void *)wb->img;
		jhinf.hdr.h_stride = sizeof(float) * 3;
		jhinf.hdr.v_stride = wb->rowsize;
		jhinf.hdr.getval = &cvt_floatColor;
		if (wb->csp->format == IPFrgb)
			comprgb2rgbWBmat(myMat, (RGBPRIMP)wb->csp->chroma,
					(RGBPRIMP)ICS_RGB709.chroma);
		else /* if (wb->csp->format == IPFxyz) */
			compxyz2rgbWBmat(myMat, (RGBPRIMP)ICS_RGB709.chroma);
		jhinf.hdr.c_data = (void *)myMat;
	} else {
		jpeghdr_src_floatRGB(&jhinf, (const float *)wb->img,
						wb->xres, wb->yres);
		jhinf.hdr.v_stride = wb->rowsize;
	}
						// assign sample-to-nits factor
	if (wb->info.flags & IIFstonits)
		jhinf.samp2nits = wb->info.stonits;
						// compute tone-mapping
	if ((n = GetImgInfoParam(&wb->info, "JHtmo", pval)) > 0 &&
			!strncasecmp(pval, "multiscale", n))
		jpeghdr_tonemap_multiscale(&jhinf);
	else
		jpeghdr_tonemap_default(&jhinf);
						// set compression params.
	jhinf.correction = JHfullsamp;
	if (wb->info.flags & IIFhdensity) {
		jhinf.cinfo.density_unit = 1;
		jhinf.cinfo.X_density = wb->info.hdensity + .5f;
		jhinf.cinfo.Y_density = wb->info.hdensity*wb->pixAspect + .5f;
	}
	if (wb->info.flags & IIFquality)
		jhinf.quality = wb->info.quality;
	else
		jhinf.quality = DEF_IQUALITY;
	if (jhinf.quality < 80)
		jhinf.correction = JHprecorr;
	if (wb->csp->format == IPFxyz) {	// attempt to preserve gamut
		jhinf.alpha = 0.67f;
		jhinf.beta = 0.75f;
	}
						// check for overrides
	if ((n = GetImgInfoParam(&wb->info, "JHcorr", pval)) >= 3) {
		if (!strncasecmp(pval, "fullsamp", n))
			jhinf.correction = JHfullsamp;
		else if (!strncasecmp(pval, "precorrect", n))
			jhinf.correction = JHprecorr;
		else if (!strncasecmp(pval, "postcorrect", n))
			jhinf.correction = JHpostcorr;
	}
	if (GetImgInfoParam(&wb->info, "JHalpha", pval))
		jhinf.alpha = atof(pval);
	if (GetImgInfoParam(&wb->info, "JHbeta", pval))
		jhinf.beta = atof(pval);
						// compress HDR image
	jpeghdr_do_compress(&jhinf);
						// clean up
	jpeghdr_destroy_compress(&jhinf);

	return ostr.bad() ? 0L : (long)ostr.tellp();
}

// Write a JPEG image from 24-bit RGB color space
static long
JPEGwriteImage24(const char *fname, const ImgWriteBuf *wb)
{
	char		pval[256];
	int		n;
						// open file for output
	ofstream	ostr(fname, ios::out|ios::trunc|ios::binary);
	if (!ostr.is_open())
		return 0;

	struct jpeg_compress_struct	cinfo;
	struct jpeg_error_mgr		jerr;
	jumper_struct			myjumper;
	myjumper.w = JMSG_NOMESSAGE;
	cinfo.err = jpeg_std_error(&jerr);
	jerr.output_message = jpeg_error_output;
	jerr.emit_message = jpeg_emit_message;
	jerr.error_exit = jpeg_error_jump;
	if (setjmp(myjumper.b)) {		// JPEG error jumps to here
		jpeg_destroy_compress(&cinfo);
		return 0;			// all clean -- return failure
	}
	cinfo.client_data = (void *)&myjumper;
	jpeg_create_compress(&cinfo);
	jpeg_stream_dest(&cinfo, &ostr);
	cinfo.image_width = wb->xres;
	cinfo.image_height = wb->yres;
	if (wb->csp->format == IPFrgb) {
		cinfo.input_components = 3;
		cinfo.in_color_space = JCS_RGB;
	} else /* wb->csp->format == IPFy */ {
		cinfo.input_components = 1;
		cinfo.in_color_space = JCS_GRAYSCALE;
	}
	jpeg_set_defaults(&cinfo);
	// Set optional compression parameters here...
	if ((n = GetImgInfoParam(&wb->info, "JCent", pval)) > 0 &&
			!strncasecmp(pval, "optimize", n))
		cinfo.optimize_coding = TRUE;
	if ((n = GetImgInfoParam(&wb->info, "JCscan", pval)) > 0 &&
			!strncasecmp(pval, "progressive", n))
		jpeg_simple_progression(&cinfo);
	if (wb->info.flags & IIFhdensity) {
		cinfo.density_unit = 1;
		cinfo.X_density = wb->info.hdensity + .5f;
		cinfo.Y_density = wb->info.hdensity*wb->pixAspect + .5f;
	}
	if (wb->info.flags & IIFquality)
		jpeg_set_quality(&cinfo, wb->info.quality, TRUE);
	else
		jpeg_set_quality(&cinfo, DEF_IQUALITY, TRUE);
		
	// Start the compressed file output
	jpeg_start_compress(&cinfo, TRUE);
	JSAMPLE *	row_pointer;		// output scanlines
	while (cinfo.next_scanline < cinfo.image_height) {
		row_pointer = (JSAMPLE *)(wb->img +
					(ssize_t)cinfo.next_scanline*wb->rowsize);
		jpeg_write_scanlines(&cinfo, &row_pointer, 1);
	}
	jpeg_finish_compress(&cinfo);		// flushes and checks stream
	jpeg_destroy_compress(&cinfo);
	return ostr.bad() ? 0L : (long)ostr.tellp();
}

// Write a JPEG image from an RGB color space
static long
JPEGwriteImage(const char *fname, const ImgWriteBuf *wb)
{
	if ((fname == NULL) | (wb == NULL) || !*fname)
		return 0;
						// check resolution
	if ((wb->xres <= 0) | (wb->yres <= 0) | (wb->img == NULL))
		return 0;
						// HDR output?
	if ((wb->csp->dtype == IDTfloat) &
			((wb->csp->format == IPFrgb) | (wb->csp->format == IPFxyz)))
		return JPEGwriteImageHDR(fname, wb);
						// LDR output?
	if ((wb->csp->dtype == IDTubyte) &
			((wb->csp->format == IPFrgb) | (wb->csp->format == IPFy)))
		return JPEGwriteImage24(fname, wb);
						// invalid color space
	return 0;
}

// Check if the given color space is supported by our writer
static const char *
JPEGsupportedCS(const ImgColorSpace *csp, int)
{
	if (csp->logorig > 0)
		return NULL;
	if ((csp->dtype == IDTfloat) &
			((csp->format == IPFrgb) | (csp->format == IPFxyz)))
		return "JPEG HDR YCbCr";
	if (csp->dtype != IDTubyte)
		return NULL;
	if (csp->format == IPFrgb)
		return "JPEG DCT YCbCr";
	if (csp->format == IPFy)
		return "JPEG DCT grayscale";
	return NULL;
}

// JPEG writer interface
extern "C" {
extern const ImgWriterInterface	IWInterfaceJPEG;
const ImgWriterInterface	IWInterfaceJPEG = {
	"jpg", &JPEGsupportedCS, &JPEGwriteImage
};
}
