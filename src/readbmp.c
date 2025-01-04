/*
 *  readbmp.c
 *  panlib
 *
 *  Read Windows BMP files.
 *
 *  Created by Greg Ward on Mon Apr 12 2004.
 *  Copyright (c) 2004 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "imgreader.h"
#include "bmpfile.h"

/* BMP image reader struct */
typedef struct {
	BASE_IMGREADER;			/* base struct members (first!) */
	BMPReader       *bmp;		/* open BMP reader */
} BMReader;

/* Open a BMP file for reading */
static ImgReader *
BMopen(const char *fname)
{
	extern const ImgReaderInterface IRInterfaceBMP;
	BMReader       *br;
	
	if (!fname || !*fname)		/* see if we can open it, first */
		return NULL;
	br = (BMReader *)calloc(1, sizeof(BMReader));
	if (br == NULL)
		return NULL;
	br->bmp = BMPopenInputFile(fname);
	if (br->bmp == NULL) {
		free(br);
		return NULL;
	}
	br->ri = &IRInterfaceBMP;       /* set reader fields */
	strcpy(br->file, fname);
	if (br->bmp->hdr->hRes > 0 && br->bmp->hdr->vRes > 0)
		br->pixAspect = (double)br->bmp->hdr->vRes/br->bmp->hdr->hRes;
	else
		br->pixAspect = 1.;
	br->xres = br->bmp->hdr->width;
	br->yres = br->bmp->hdr->height;
	if (BMPisGrayscale(br->bmp->hdr))
		PcopyCS(&br->cs, &ICS_Y8);
	else
		PcopyCS(&br->cs, &ICS_sRGB);
	switch (br->bmp->hdr->bpp) {
	case 32:
		br->encoding = "32-bit RGB";
		break;
	case 24:
		br->encoding = "24-bit RGB";
		break;
	case 16:
		br->encoding = "16-bit RGB";
		break;
	case 8:
		br->encoding = br->cs.format==IPFy ?
				"8-bit grayscale" : "8-bit color";
		break;
	case 4:
		br->encoding = br->cs.format==IPFy ?
				"4-bit grayscale" : "4-bit color";
		break;
	case 1:
		br->encoding = "bilevel bitmap";
		break;
	default:
		br->errCode = IREunsupported;
		sprintf(br->errMsg, "unsupported bit depth (%d)",
				br->bmp->hdr->bpp);
		return (ImgReader *)br;
	}
	br->frame = 0; br->nframes = 1;
	br->frameType = IRFnone; br->frameRate = 0;
	br->fr.xleft = 0; br->fr.xright = br->xres;
	if (br->bmp->hdr->yIsDown) {
		br->fr.ytop = 0; br->fr.ybottom = 1;
	} else {
		br->fr.ytop = br->yres-1; br->fr.ybottom = br->yres;
	}
	br->nr = br->fr;		/* we're ready to roll */
	return (ImgReader *)br;
}

/* Close our BMP reader */
static void
BMclose(ImgReader *ir)
{
	BMReader	*br = (BMReader *)ir;
	if (br == NULL)
		return;
	if (br->bmp != NULL)
		BMPcloseInput(br->bmp);
	free(br);
}

/* Convert our scanline once it's been read */
static void
BMconvertScanline(const BMReader *br, ImgReadBuf *rb)
{
	int     psiz = ImgPixelSize(&rb->cs);
	uby8    *pp = rb->buf;
	int     y = br->bmp->yscan;
	int     x;
	if (!br->bmp->hdr->yIsDown)
		y = br->yres-1 - y;
	if (y > rb->r.ytop)
		pp += psiz * ((rb->r.xright - rb->r.xleft)/rb->subsample) *
				((y - rb->r.ytop)/rb->subsample);
	for (x = rb->r.xleft; x < rb->r.xright; x += rb->subsample) {
		RGBquad rgbq = BMPdecodePixel(x, br->bmp);
		if (rb->cs.format == IPFrgb) {
			*pp++ = rgbq.r;
			*pp++ = rgbq.g;
			*pp++ = rgbq.b;
		} else /* rb->cs.format == IPFy */
			*pp++ = rgbq.g;
	}
}

/* Interpret BMP library error */
static ImgReadErr
BMinterpErr(BMReader *br, int ec)
{
	switch (ec) {
	case BIR_OK:
		return IREnone;
	case BIR_EOF:
		strcpy(br->errMsg, "unexpected read past end");
		return br->errCode = IREunknown;
	case BIR_TRUNCATED:
		br->errCode = IREtruncated;
		break;
	case BIR_UNSUPPORTED:
		br->errCode = IREunsupported;
		break;
	case BIR_RLERROR:
		br->errCode = IREformat;
		break;
	case BIR_SEEKERR:
		br->errCode = IREread;
		break;
	default:
		br->errCode = IREunknown;
		break;
	}
	strcpy(br->errMsg, BMPerrorMessage(ec));
	return br->errCode;
}

/* Read the indicated record from open BMP file */
static ImgReadErr
BMreadRec(ImgReader *ir, ImgReadBuf *rb)
{
	BMReader	*br = (BMReader *)ir;
	int		y, ystep, yend;
	if (br == NULL || br->bmp == NULL || rb == NULL)
		return IREunknown;

	if (!PmatchColorSpace(&rb->cs, &br->cs, PICMptype)) {
		rb->cs = br->cs;		/* make caller work */
		rb->buf = NULL;
	}
	if (!PlegalRect(&rb->r, br->xres, br->yres)) {
		strcpy(br->errMsg, "illegal read rectangle");
		return br->errCode = IREunknown;
	}
	ImgFixSampling(rb);			/* fix rectangle bounds */
	if (!rb->buf) {				/* need new buffer */
		int     buflen;
		rb->r.xleft = 0;		/* use scanline-sized buffer */
		rb->r.xright = br->xres;
		rb->r.ybottom = rb->r.ytop + rb->subsample;
		ImgFixSampling(rb);
		buflen = ImgReadBufLen(rb);
		if (buflen <= 0) {
			strcpy(br->errMsg, "internal buffer size error");
			return br->errCode = IREunknown;
		}
						/* do not free client memory */
		rb->buf = (uby8 *)malloc(buflen);
		if (rb->buf == NULL) {
			strcpy(br->errMsg, "cannot allocate new read buffer");
			return br->errCode = IREmemory;
		}
	}
						/* establish read direction */
	if (br->bmp->hdr->yIsDown) {
		y = rb->r.ytop;
		ystep = rb->subsample;
		yend = rb->r.ybottom;
	} else {
		y = rb->r.ybottom - 1;
		ystep = -rb->subsample;
		yend = rb->r.ytop - 1;
	}
						/* read our rectangle */
	for ( ; y != yend; y += ystep) {
		int     ec = BMPseekScanline(br->bmp->hdr->yIsDown ?
						y : br->yres-1-y, br->bmp);
		if (ec == BIR_OK)
			BMconvertScanline(br, rb);
		else
			return BMinterpErr(br, ec);

		if (br->bmp->hdr->yIsDown)
			br->nr.ytop = br->nr.ybottom++;
		else
			br->nr.ybottom = br->nr.ytop--;
	}
	return IREnone;				/* all is well */
}

/* Interface for BMP image reader */
const ImgReaderInterface IRInterfaceBMP = {
	"BMP.rle",
	&BMopen, NULL, NULL, NULL,
	NULL, &BMreadRec, &BMclose
};
