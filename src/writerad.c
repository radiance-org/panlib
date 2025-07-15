/*
 *  writerad.c
 *  panlib
 *
 *  Write Radiance (RGBE or XYZE) image.
 *
 *  Created by gward on Fri Sep 07 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#include <math.h>
#include "imgio.h"
#include "imgwriter.h"
#include "radheader.h"
#include "resolu.h"
#include "color.h"

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

/* Write out a Radiance information header */
static void
RadWriteHeader(const ImgWriteBuf *wb, FILE *pout)
{
	RHstartHeader(&wb->info, pout);		/* create info header */
	RHoptViewAngle(wb, pout);

	if (wb->info.flags & IIFstonits)	/* add EXPOSURE */
		fputexpos((wb->csp->format==IPFxyz ? 1. : WHTEFFICACY) /
				wb->info.stonits, pout);
						/* primary chromaticities */
	if (wb->csp->format == IPFrgb)
		fputprims((RGBPRIMP)wb->csp->chroma, pout);
						/* pixel aspect ratio */
	if ((wb->pixAspect < .995) | (wb->pixAspect > 1.005))
		fputaspect(1./wb->pixAspect, pout);
						/* pixel format */
	fputformat((char *)(wb->csp->format==IPFxyz ? CIEFMT : COLRFMT), pout);
						/* end of info. header */
	fputc('\n', pout);
}

/* Write out a Radiance RGBE or XYZE picture */
static long
RadWriteImage(const char *fname, const ImgWriteBuf *wb)
{
	COLR *		clrs = NULL;
	int		orient = PIXSTANDARD;
	FILE *		pout;
	const uby8 *	slpos;
	int		y;
	long		flen;

	if ((fname == NULL) | (wb == NULL) ||
			!*fname | (wb->img == NULL) | (wb->csp == NULL))
		return 0;
						/* check resolution */
	if ((wb->xres <= 0) | (wb->yres <= 0))
		return 0;
						/* verify color space */
	if (wb->csp->dtype != IDTfloat || (wb->csp->format != IPFrgb) &
						(wb->csp->format != IPFxyz))
		return 0;
						/* open output file */
	if ((pout = fopen(fname, "wb")) == NULL)
		return 0;
						/* write header + resolution */
	RadWriteHeader(wb, pout);
	if (wb->info.flags & IIForientation)
		orient = RHortab[wb->info.orientation-1];
	fputresolu(orient, wb->xres, wb->yres, pout);
						/* allocate buffer if flat */
	if (wb->info.flags & IIFquality && wb->info.quality >= 100)
		clrs = (COLR *)tempbuffer(sizeof(COLR)*wb->xres);
						/* write scanlines */
	for (slpos = wb->img, y = 0; y < wb->yres; slpos += wb->rowsize, y++) {
		if (clrs != NULL) {
			const float *	clr = (const float *)slpos;
			COLOR		gclr;
			int		x;
			if (wb->csp->format == IPFxyz)
				for (x = 0; x < wb->xres; x++, clr += 3)
					setcolr(clrs[x], clr[CIEX],
							clr[CIEY], clr[CIEZ]);
			else
				for (x = 0; x < wb->xres; x++, clr += 3) {
					copycolor(gclr, clr);
					clipgamut(gclr, bright(gclr),
							CGAMUT_LOWER,
							cblack, cwhite);
					setcolr(clrs[x], gclr[RED],
							gclr[GRN], gclr[BLU]);
				}
			if ((int)fwrite(clrs, sizeof(COLR), wb->xres, pout) != wb->xres)
				goto writerr;
		} else if (fwritescan((COLOR *)slpos, wb->xres, pout) < 0)
			goto writerr;
	}
	if (fflush(pout) == EOF)		/* clean up */
		goto writerr;
	flen = ftell(pout);
	fclose(pout);
	return flen;				/* return file length */
writerr:
	fclose(pout);
	return 0;
}

/* Check if the given color space is supported by our writer */
static const char *
RadSupportedCS(const ImgColorSpace *csp, int qual)
{
	if ((csp->gamma != 1) | (csp->logorig > 0))
		return NULL;
	if (csp->dtype != IDTfloat)
		return NULL;
	if (qual >= 100) {
		if (csp->format == IPFrgb)
			return "Radiance 32-bit RGBE";
		if (csp->format == IPFxyz)
			return "Radiance 32-bit XYZE";
		return NULL;
	}
	if (csp->format == IPFrgb)
		return "Radiance RLE RGBE";
	if (csp->format == IPFxyz)
		return "Radiance RLE XYZE";
	return NULL;
}

/* Radiance writer interface */
const ImgWriterInterface	IWInterfaceRad = {
	"hdr", &RadSupportedCS, &RadWriteImage
};
