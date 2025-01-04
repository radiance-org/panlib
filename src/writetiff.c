/*
 *  writetiff.c
 *  panlib
 *
 *  Routines for writing out TIFF image files.
 *
 *  Created by gward on Thu Oct 04 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "system.h"
#include "tiffio.h"
#include "tiffmsg.h"
#include "color.h"
#include "imgio.h"
#include "imgwriter.h"

typedef int	(*TiffScanWriter)(TIFF *, const char *, int);

typedef struct scanConversion {
	int		fno;		/* open file descriptor (ID) */
	float *		scanbuf;	/* conversion/output buffer */
	COLORMAT	mat;		/* conversion matrix */
	struct scanConversion *
			next;		/* next conversion structure */
} TScanConversion;

static TScanConversion *	tscList = NULL;

/* Set up RGB->XYZ conversion matrix and allocate buffer */
static int
TiffSetConversion(TIFF *tfo, const ImgWriteBuf *wb)
{
	const int		idmatch = TIFFFileno(tfo);
	TScanConversion *	tsc;
					/* find/allocate structure */
	for (tsc = tscList; tsc != NULL; tsc = tsc->next)
		if (tsc->fno == idmatch)
			break;
	if (tsc == NULL) {
		tsc = (TScanConversion *)malloc(sizeof(TScanConversion));
		if (tsc == NULL)
			return 0;
		tsc->fno = idmatch;
		tsc->next = tscList;
		tscList = tsc;
	} else
		free(tsc->scanbuf);
	tsc->scanbuf = (float *)malloc(3*sizeof(float)*wb->xres);
	if (tsc->scanbuf == NULL)
		return 0;
					/* compute conversion matrix */
	comprgb2xyzWBmat(tsc->mat, (const RGBPRIMP)wb->csp->chroma);
	return 1;
}

/* Free conversion buffer and close TIFF */
static void
TiffDone(TIFF *tfo)
{
	const int		idmatch = TIFFFileno(tfo);
	TScanConversion *	tsclast = NULL;
	TScanConversion *	tsc;
					/* find conversion structure */
	for (tsc = tscList; tsc != NULL; tsc = tsc->next) {
		if (tsc->fno == idmatch)
			break;
		tsclast = tsc;
	}
	if (tsc != NULL) {		/* free resources & unlink */
		free(tsc->scanbuf);
		if (tsclast != NULL)
			tsclast->next = tsc->next;
		else
			tscList = tsc->next;
		free(tsc);
	}
	TIFFClose(tfo);			/* close TIFF */
}

/* Convert scanline buffer to output color space and write it out */
static int
TiffConvertScan(TIFF *tfo, const char *buf, int y)
{
	const int		idmatch = TIFFFileno(tfo);
	TScanConversion *	tsc;
	int32			xpos;
					/* find conversion structure */
	for (tsc = tscList; tsc != NULL; tsc = tsc->next)
		if (tsc->fno == idmatch)
			break;
	if (tsc == NULL)
		return 0;
					/* convert scanline */
	if (!TIFFGetField(tfo, TIFFTAG_IMAGEWIDTH, &xpos))
		return 0;
	for (xpos *= 3; (xpos -= 3) >= 0; )
		colortrans(tsc->scanbuf + xpos, tsc->mat, (float *)buf + xpos);
					/* write converted scanline */
	return (TIFFWriteScanline(tfo, (tdata_t)tsc->scanbuf, (uint32)y, 0) >= 0);
}

/* Write opaque scanline buffer, assuming everything's set up already */
static int
TiffWriteScanBuf(TIFF *tfo, const char *buf, int y)
{
	return (TIFFWriteScanline(tfo, (tdata_t)buf, (uint32)y, 0) >= 0);
}

/* Software version hack for 16-bit/sample linear gamma */
static const char *
gammaHack(const ImgWriteBuf *wb)
{
	const char *	defVers = NULL;
	if (wb->info.flags & IIFsource) {
		defVers = wb->info.source.vers;
		if (strstr(defVers, "dcraw") != NULL)
			return defVers;
	}
	if (wb->csp->dtype != IDTushort)
		return defVers;
	if ((wb->csp->logorig > 0) | (wb->csp->gamma > 1.1f))
		return defVers;
	return "dcraw";
}

/* Set up TIFF writer by setting appropriate fields */
static TiffScanWriter
TiffInitWrite(TIFF *tfo, const ImgWriteBuf *wb)
{
	int		qual = DEF_IQUALITY;
	const char *	vers = gammaHack(wb);
	uint32		rows_strip = 8192/(wb->xres*ImgPixelSize(wb->csp));
	float		hdens = 72.f;
					/* check image dimensions */
	if ((wb->xres <= 0) | (wb->yres <= 0))
		return NULL;
					/* get quality */
	if (wb->info.flags & IIFquality)
		qual = wb->info.quality;
					/* basic fields */
	TIFFSetField(tfo, TIFFTAG_IMAGEWIDTH, (uint32)wb->xres);
	TIFFSetField(tfo, TIFFTAG_IMAGELENGTH, (uint32)wb->yres);
	TIFFSetField(tfo, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	if (!rows_strip) rows_strip = 1;
	else if (rows_strip > wb->yres) rows_strip = wb->yres;
	TIFFSetField(tfo, TIFFTAG_ROWSPERSTRIP, rows_strip);
	TIFFSetField(tfo, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
	TIFFSetField(tfo, TIFFTAG_RESOLUTIONUNIT, RESUNIT_INCH);
	if (wb->info.flags & IIFhdensity)
		hdens = wb->info.hdensity;
	TIFFSetField(tfo, TIFFTAG_XRESOLUTION, hdens);
	TIFFSetField(tfo, TIFFTAG_YRESOLUTION, (float)(hdens * wb->pixAspect));
	if (wb->info.flags & IIForientation)
		TIFFSetField(tfo, TIFFTAG_ORIENTATION, (short)wb->info.orientation);
					/* info fields */
	if (wb->info.flags & IIFcapdate)
		TIFFSetField(tfo, TIFFTAG_DATETIME, wb->info.capdate);
	if (wb->info.flags & IIFsource) {
		TIFFSetField(tfo, TIFFTAG_MAKE, wb->info.source.make);
		TIFFSetField(tfo, TIFFTAG_MODEL, wb->info.source.model);
	}
	if (vers != NULL)
		TIFFSetField(tfo, TIFFTAG_SOFTWARE, vers);
	if (wb->info.flags & IIFstonits)
		TIFFSetField(tfo, TIFFTAG_STONITS, wb->info.stonits);
	if (wb->info.flags & IIFowner)
		TIFFSetField(tfo, TIFFTAG_ARTIST, wb->info.owner);
	if (wb->info.flags & (IIFparams|IIFcomments)) {
		char	mydescr[4096];
		char *	cp;
		mydescr[0] = '\0';
		if (wb->info.flags & IIFparams)
			strcpy(mydescr, wb->info.params);
		if (wb->info.flags & IIFcomments)
			strcat(mydescr, wb->info.comments);
		for (cp = mydescr; *cp; cp++)
			if (*cp == '\n')
				*cp = cp[1] ? ';' : '\0';
		TIFFSetField(tfo, TIFFTAG_IMAGEDESCRIPTION, mydescr);
	}
					/* 24-bit RGB color space */
	if (wb->csp->format == IPFrgb && wb->csp->dtype == IDTubyte) {
		TIFFSetField(tfo, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
		TIFFSetField(tfo, TIFFTAG_SAMPLESPERPIXEL, 3);
		TIFFSetField(tfo, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
		TIFFSetField(tfo, TIFFTAG_BITSPERSAMPLE, 8);
		TIFFSetField(tfo, TIFFTAG_PRIMARYCHROMATICITIES, wb->csp->chroma);
		TIFFSetField(tfo, TIFFTAG_WHITEPOINT, wb->csp->chroma+3);
		if (qual < 100)
			TIFFSetField(tfo, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
		return &TiffWriteScanBuf;
	}
					/* 8-bit grayscale */
	if (wb->csp->format == IPFy && wb->csp->dtype == IDTubyte) {
		TIFFSetField(tfo, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
		TIFFSetField(tfo, TIFFTAG_SAMPLESPERPIXEL, 1);
		TIFFSetField(tfo, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
		TIFFSetField(tfo, TIFFTAG_BITSPERSAMPLE, 8);
		if (qual < 100)
			TIFFSetField(tfo, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
		return &TiffWriteScanBuf;
	}
					/* 48-bit RGB color space */
	if (wb->csp->format == IPFrgb && wb->csp->dtype == IDTushort) {
		TIFFSetField(tfo, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
		TIFFSetField(tfo, TIFFTAG_SAMPLESPERPIXEL, 3);
		TIFFSetField(tfo, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
		TIFFSetField(tfo, TIFFTAG_BITSPERSAMPLE, 16);
		TIFFSetField(tfo, TIFFTAG_PRIMARYCHROMATICITIES, wb->csp->chroma);
		TIFFSetField(tfo, TIFFTAG_WHITEPOINT, wb->csp->chroma+3);
		return &TiffWriteScanBuf;
	}
					/* 16-bit grayscale */
	if (wb->csp->format == IPFy && wb->csp->dtype == IDTushort) {
		TIFFSetField(tfo, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
		TIFFSetField(tfo, TIFFTAG_SAMPLESPERPIXEL, 1);
		TIFFSetField(tfo, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
		TIFFSetField(tfo, TIFFTAG_BITSPERSAMPLE, 16);
		return &TiffWriteScanBuf;
	}
					/* floating-point XYZ -> 24-bit LogLuv */
	if (qual <= 50 && wb->csp->format == IPFxyz && wb->csp->dtype == IDTfloat) {
		rows_strip = 8192/(wb->xres*3);
		if (!rows_strip) rows_strip = 1;
		else if ((int)rows_strip > wb->yres) rows_strip = wb->yres;
		TIFFSetField(tfo, TIFFTAG_ROWSPERSTRIP, rows_strip);
		TIFFSetField(tfo, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_LOGLUV);
		TIFFSetField(tfo, TIFFTAG_COMPRESSION, COMPRESSION_SGILOG24);
		TIFFSetField(tfo, TIFFTAG_SGILOGDATAFMT, SGILOGDATAFMT_FLOAT);
		TIFFSetField(tfo, TIFFTAG_SAMPLESPERPIXEL, 3);
		TIFFSetField(tfo, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
		TIFFSetField(tfo, TIFFTAG_BITSPERSAMPLE, 8*sizeof(float));
		return &TiffWriteScanBuf;
	}
					/* floating-point XYZ -> LogLuv */
	if (wb->csp->format == IPFxyz && wb->csp->dtype == IDTfloat) {
		rows_strip = 8192/(wb->xres*4);
		if (!rows_strip) rows_strip = 1;
		else if ((int)rows_strip > wb->yres) rows_strip = wb->yres;
		TIFFSetField(tfo, TIFFTAG_ROWSPERSTRIP, rows_strip);
		TIFFSetField(tfo, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_LOGLUV);
		TIFFSetField(tfo, TIFFTAG_COMPRESSION, COMPRESSION_SGILOG);
		TIFFSetField(tfo, TIFFTAG_SGILOGDATAFMT, SGILOGDATAFMT_FLOAT);
		TIFFSetField(tfo, TIFFTAG_SAMPLESPERPIXEL, 3);
		TIFFSetField(tfo, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
		TIFFSetField(tfo, TIFFTAG_BITSPERSAMPLE, 8*sizeof(float));
		return &TiffWriteScanBuf;
	}
					/* floating point grayscale */
	if (qual >= 100 && wb->csp->format == IPFy && wb->csp->dtype == IDTfloat) {
		TIFFSetField(tfo, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
		TIFFSetField(tfo, TIFFTAG_SAMPLESPERPIXEL, 1);
		TIFFSetField(tfo, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
		TIFFSetField(tfo, TIFFTAG_BITSPERSAMPLE, 8*sizeof(float));
		return &TiffWriteScanBuf;
	}
					/* floating point RGB */
	if (qual >= 100 && wb->csp->format == IPFrgb && wb->csp->dtype == IDTfloat) {
		TIFFSetField(tfo, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
		TIFFSetField(tfo, TIFFTAG_SAMPLESPERPIXEL, 3);
		TIFFSetField(tfo, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
		TIFFSetField(tfo, TIFFTAG_BITSPERSAMPLE, 8*sizeof(float));
		TIFFSetField(tfo, TIFFTAG_PRIMARYCHROMATICITIES, wb->csp->chroma);
		TIFFSetField(tfo, TIFFTAG_WHITEPOINT, wb->csp->chroma+3);
		return &TiffWriteScanBuf;
	}
					/* floating-point grayscale -> LogL */
	if (wb->csp->format == IPFy && wb->csp->dtype == IDTfloat) {
		rows_strip = 8192/(wb->xres*2);
		if (!rows_strip) rows_strip = 1;
		else if ((int)rows_strip > wb->yres) rows_strip = wb->yres;
		TIFFSetField(tfo, TIFFTAG_ROWSPERSTRIP, rows_strip);
		TIFFSetField(tfo, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_LOGL);
		TIFFSetField(tfo, TIFFTAG_COMPRESSION, COMPRESSION_SGILOG);
		TIFFSetField(tfo, TIFFTAG_SGILOGDATAFMT, SGILOGDATAFMT_FLOAT);
		TIFFSetField(tfo, TIFFTAG_SAMPLESPERPIXEL, 1);
		TIFFSetField(tfo, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
		TIFFSetField(tfo, TIFFTAG_BITSPERSAMPLE, 8*sizeof(float));
		return &TiffWriteScanBuf;
	}
					/* floating-point RGB -> 24-bit LogLuv */
	if (qual <= 50 && wb->csp->format == IPFrgb && wb->csp->dtype == IDTfloat) {
		rows_strip = 8192/(wb->xres*3);
		if (!rows_strip) rows_strip = 1;
		else if ((int)rows_strip > wb->yres) rows_strip = wb->yres;
		TIFFSetField(tfo, TIFFTAG_ROWSPERSTRIP, rows_strip);
		TIFFSetField(tfo, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_LOGLUV);
		TIFFSetField(tfo, TIFFTAG_COMPRESSION, COMPRESSION_SGILOG24);
		TIFFSetField(tfo, TIFFTAG_SGILOGDATAFMT, SGILOGDATAFMT_FLOAT);
		TIFFSetField(tfo, TIFFTAG_SAMPLESPERPIXEL, 3);
		TIFFSetField(tfo, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
		TIFFSetField(tfo, TIFFTAG_BITSPERSAMPLE, 8*sizeof(float));
		if (TiffSetConversion(tfo, wb))
			return &TiffConvertScan;
	}
					/* floating-point RGB -> LogLuv */
	if (wb->csp->format == IPFrgb && wb->csp->dtype == IDTfloat) {
		rows_strip = 8192/(wb->xres*4);
		if (!rows_strip) rows_strip = 1;
		else if ((int)rows_strip > wb->yres) rows_strip = wb->yres;
		TIFFSetField(tfo, TIFFTAG_ROWSPERSTRIP, rows_strip);
		TIFFSetField(tfo, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_LOGLUV);
		TIFFSetField(tfo, TIFFTAG_COMPRESSION, COMPRESSION_SGILOG);
		TIFFSetField(tfo, TIFFTAG_SGILOGDATAFMT, SGILOGDATAFMT_FLOAT);
		TIFFSetField(tfo, TIFFTAG_SAMPLESPERPIXEL, 3);
		TIFFSetField(tfo, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
		TIFFSetField(tfo, TIFFTAG_BITSPERSAMPLE, 8*sizeof(float));
		if (TiffSetConversion(tfo, wb))
			return &TiffConvertScan;
	}
	return NULL;			/* unsupported color space */
}

/* Write out a TIFF image */
static long
TiffWriteImage(const char *fname, const ImgWriteBuf *wb)
{
	off_t		flen;
	TIFF *		tfo;
	TiffScanWriter	scanwrite;
	int		y;
	const char *	bpos;
					/* check arguments */
	if ((fname == NULL) | (wb == NULL) ||
			!*fname | (wb->img == NULL) | (wb->csp == NULL))
		return 0;
	/* replace TIFF stderr reporting with buffer writers */
 	TIFFSetErrorHandler(BufferTiffError);
	TIFFSetWarningHandler(BufferTiffWarning);
	TiffMessageBuffer[0] = '\0';
	tfo = TIFFOpen(fname, "w");
	if (tfo == NULL)
		return 0;
	flen = 0;
	scanwrite = TiffInitWrite(tfo, wb);
	if (scanwrite == NULL)
		goto cleanup;
					/* write each scanline */
	for (y = 0, bpos = (char *)wb->img; y < wb->yres;
					y++, bpos += wb->rowsize)
		if (!(*scanwrite)(tfo, bpos, y))
			goto cleanup;
	if (!TIFFFlush(tfo))		/* flush buffers */
		goto cleanup;
	flen = lseek(TIFFFileno(tfo), 0L, SEEK_END);
cleanup:
	TiffDone(tfo);
	return (long)flen;
}

/* Check if the given color space is supported by our writer */
static const char *
TiffSupportedCS(const ImgColorSpace *csp, int qual)
{
	if (qual < 0)
		qual = DEF_IQUALITY;
	if (csp->logorig > 0)
		return NULL;
	if (qual >= 100 && (csp->dtype == IDTubyte) & (csp->format == IPFrgb))
		return "TIFF 24-bit RGB";
	if ((csp->dtype == IDTubyte) & (csp->format == IPFrgb))
		return "TIFF LZW RGB";
	if (qual >= 100 && (csp->dtype == IDTubyte) & (csp->format == IPFy))
		return "TIFF 8-bit grayscale";
	if ((csp->dtype == IDTubyte) & (csp->format == IPFy))
		return "TIFF LZW grayscale";
	if ((csp->dtype == IDTushort) & (csp->format == IPFrgb))
		return "TIFF 48-bit RGB";
	if ((csp->dtype == IDTushort) & (csp->format == IPFy))
		return "TIFF 16-bit grayscale";
	if (qual >= 100 && (csp->dtype == IDTfloat) & (csp->format == IPFrgb))
		return "TIFF 96-bit RGB";
	if (qual >= 100 && (csp->dtype == IDTfloat) & (csp->format == IPFy))
		return "TIFF 32-bit grayscale";
	if (qual <= 50 && csp->dtype == IDTfloat &&
			(csp->format == IPFxyz) | (csp->format == IPFrgb))
		return "TIFF 24-bit-Log Luv";
	if (csp->dtype == IDTfloat &&
			(csp->format == IPFxyz) | (csp->format == IPFrgb))
		return "TIFF RLE-Log Luv";
	if ((csp->dtype == IDTfloat) & (csp->format == IPFy))
		return "TIFF RLE-Log luminance";
	return NULL;
}	

const ImgWriterInterface	IWInterfaceTIFF = {
	"tif", &TiffSupportedCS, &TiffWriteImage
};
