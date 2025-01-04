/*
 *  writedpt.c
 *  panlib
 *
 *  Write 16-bit encoded depth image.
 *
 *  Created by gward on July 27 2019.
 *  Copyright (c) 2019 Anyhere Software. All rights reserved.
 *
 */

#include <math.h>
#include <stdlib.h>
#include "imgio.h"
#include "imgwriter.h"
#include "radheader.h"
#include "fvect.h"
#include "depthcodec.h"

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

/* Write out Radiance information header */
static void
DepthWriteHeader(const ImgWriteBuf *wb, FILE *pout)
{
	RHstartHeader(&wb->info, pout);		/* create info header */
						/* depth format */
	fputformat((char *)DEPTH16FMT, pout);
						/* end of info. header */
	fputc('\n', pout);
}

/* Write out a 16-bit encoded depth map */
static long
DepthWriteImage(const char *fname, const ImgWriteBuf *wb)
{
	int		orient = PIXSTANDARD;
	double		refDepth = 1;
	const char *	ref_depth_unit;
	FILE *		pout;
	int		x, y;
	long		flen;

	if ((fname == NULL) | (wb == NULL) ||
			!*fname | (wb->img == NULL) | (wb->csp == NULL))
		return 0;
						/* check resolution */
	if ((wb->xres <= 0) | (wb->yres <= 0))
		return 0;
						/* verify color space */
	if ((wb->csp->dtype != IDTfloat) | (ImgPixelLen[wb->csp->format] != 1))
		return 0;
						/* get reference depth */
	ref_depth_unit = FindImgInfoParam(&wb->info, DEPTHSTR);
	if (ref_depth_unit) {
		ref_depth_unit += LDEPTHSTR;
		refDepth = atof(ref_depth_unit);
		if (refDepth <= 0) refDepth = 1;
	}
						/* open output file */
	if ((pout = fopen(fname, "wb")) == NULL)
		return 0;
						/* write header + resolution */
	DepthWriteHeader(wb, pout);
	if (wb->info.flags & IIForientation)
		orient = RHortab[wb->info.orientation-1];
	fputresolu(orient, wb->xres, wb->yres, pout);
						/* write each scanline */
	for (y = 0; y < wb->yres; y++) {
		const float *	dp = (const float *)(wb->img + y*wb->rowsize);
		for (x = wb->xres; x--; dp++)
			if (putint(depth2code(*dp, refDepth), 2, pout) == EOF) {
				fclose(pout);
				return 0;
			}
	}
	flen = ftell(pout);			/* return file length */
	return flen * (fclose(pout) != EOF);	/* or 0 on write error */
}

/* Check if the given color space is supported by our writer */
static const char *
DepthSupportedCS(const ImgColorSpace *csp, int qual)
{
	if ((csp->gamma != 1) | (csp->logorig > 0))
		return NULL;
	if (csp->dtype != IDTfloat)
		return NULL;
	if ((csp->format != IPFd) & (csp->format != IPFy))
		return NULL;

	return "16-bit encoded depth";
}

/* Depth writer interface */
const ImgWriterInterface	IWInterfaceDPT = {
	"dpt", &DepthSupportedCS, &DepthWriteImage
};
