/*
 *  writemtx.c
 *  panlib
 *
 *  Write floating-point RGB, XYZ, or Y map as Radiance-style float matrix.
 *
 *  Created by Greg Ward on 3/29/20.
 *  Copyright 2020 Anyhere Software. All rights reserved.
 *
 */

#include <math.h>
#include "imgio.h"
#include "imgwriter.h"
#include "radheader.h"
#include "color.h"

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

/* Write out Radiance information header */
static void
MatrixWriteHeader(const ImgWriteBuf *wb, FILE *fout)
{
	RHstartHeader(&wb->info, fout);		/* create info header */
	RHoptViewAngle(wb, fout);
						/* color space */
	if (wb->csp->format == IPFrgb)
		fputprims((RGBPRIMP)wb->csp->chroma, fout);
	else if (wb->csp->format == IPFxyz)
		fputprims((RGBPRIMP)ICS_XYZ.chroma, fout);
						/* non-square pixels? */
	if ((wb->pixAspect < 0.99) | (wb->pixAspect > 1.01))
		fputaspect(1./wb->pixAspect, fout);
						/* dimensions, matrix-style */
	fprintf(fout, "NROWS=%d\nNCOLS=%d\nNCOMP=%d\n", wb->yres, wb->xres,
			ImgPixelLen[wb->csp->format]);

	fputendian(fout);			/* output byte order+format */
	fputformat((char *)"float", fout);

	fputc('\n', fout);			/* end of info. header */
}

/* Do actual write to file and close, returning file length, or 0 on error */
static long
fputmatrix(const ImgWriteBuf *wb, const int xstride, const int ystride,
			const float sf, FILE *fp)
{
	const int	psiz = ImgPixelSize(wb->csp);
	const int	doscale = ((sf < 0.999f) | (sf > 1.001f)) *
					ImgPixelLen[wb->csp->format];
	const uby8	*img = wb->img;
	int		y;
						/* adjust image origin */
	if ((xstride == -wb->rowsize) | (ystride == -wb->rowsize))
		img += (wb->yres - 1)*wb->rowsize;
	if ((xstride == -psiz) | (ystride == -psiz))
		img += (wb->xres - 1)*psiz;
	if (!doscale & (xstride == psiz)) {	/* optimize common case */
	    for (y = 0; y < wb->yres; y++)
		if (putbinary(img + y*ystride, psiz, wb->xres, fp) != wb->xres)
		    goto writerr;
	} else {				/* general case */
	    int		x, p;
	    for (y = 0; y < wb->yres; y++)
		for (x = 0; x < wb->xres; x++) {
		    float	val[MaxPixelLen];
		    memcpy(val, img + x*xstride + y*ystride, psiz);
		    for (p = doscale; p--; )
			val[p] *= sf;
		    if (putbinary(val, psiz, 1, fp) != 1)
			goto writerr;
		}
	}
	if (fflush(fp) == 0) {			/* finish output */
		long	flen = ftell(fp);
		fclose(fp);
		return flen;
	}
writerr:
	fclose(fp);
	return 0;
}

/* Write out a 32-bit float map */
static long
MatrixWriteImage(const char *fname, const ImgWriteBuf *wb)
{
	ImgOrientation	orient = IOtopleft;
	double		calib = 1;
	int		psiz;
	FILE *		fout;

	if ((fname == NULL) | (wb == NULL) ||
			!*fname | (wb->img == NULL) | (wb->csp == NULL))
		return 0;
						/* check resolution */
	if ((wb->xres <= 0) | (wb->yres <= 0))
		return 0;
						/* verify color space */
	if (wb->csp->dtype != IDTfloat || (wb->csp->format != IPFy) &
			(wb->csp->format != IPFxyz) & (wb->csp->format != IPFrgb))
		return 0;
						/* open output file */
	if ((fout = fopen(fname, "wb")) == NULL)
		return 0;
						/* write header, resolution */
	MatrixWriteHeader(wb, fout);
	
	if (wb->info.flags & IIFstonits) {
		calib = wb->info.stonits;
		if (wb->csp->format == IPFrgb)	/* hack for watts/sr/m2 */
			calib *= 1./WHTEFFICACY;
	}
	if (wb->info.flags & IIForientation)
		orient = wb->info.orientation;

	psiz = ImgPixelSize(wb->csp);

	switch (orient) {			/* write standard orientation */
	case IOtopleft:
		return fputmatrix(wb, psiz, wb->rowsize, calib, fout);
	case IOrighttop:
		return fputmatrix(wb, -wb->rowsize, psiz, calib, fout);
	case IObotright:
		return fputmatrix(wb, -psiz, -wb->rowsize, calib, fout);
	case IOleftbot:
		return fputmatrix(wb, wb->rowsize, -psiz, calib, fout);
	case IOtopright:
		return fputmatrix(wb, -psiz, wb->rowsize, calib, fout);
	case IObotleft:
		return fputmatrix(wb, psiz, -wb->rowsize, calib, fout);
	case IOlefttop:
		return fputmatrix(wb, wb->rowsize, psiz, calib, fout);
	case IOrightbot:
		return fputmatrix(wb, -wb->rowsize, -psiz, calib, fout);
	}
	return 0;				/* should never be reached */
}

/* Check if the given color space is supported by our writer */
static const char *
MatrixSupportedCS(const ImgColorSpace *csp, int qual)
{
	if ((csp->gamma != 1) | (csp->logorig > 0))
		return NULL;
	if (csp->dtype != IDTfloat)
		return NULL;

	switch (csp->format) {
	case IPFy:
		return "32-bit-float luminance matrix";
	case IPFxyz:
		return "32-bit-float CIE XYZ matrix";
	case IPFrgb:
		return "32-bit-float RGB radiance matrix";
	default:
		break;
	}
	return NULL;
}

/* Depth writer interface */
const ImgWriterInterface	IWInterfaceMTX = {
	"mtx", &MatrixSupportedCS, &MatrixWriteImage
};
