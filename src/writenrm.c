/*
 *  writenrm.c
 *  panlib
 *
 *  Write 32-bit encoded direction image.
 *
 *  Created by gward on July 27 2019.
 *  Copyright (c) 2019 Anyhere Software. All rights reserved.
 *
 */

#define _USE_MATH_DEFINES
#include <math.h>
#include "system.h"
#include "imgio.h"
#include "imgwriter.h"
#include "radheader.h"
#include "rtmath.h"
#include "normcodec.h"

/* Check if the given color space is supported by our writer */
static const char *
NormSupportedCS(const ImgColorSpace *csp, int qual)
{
	if (!PmatchColorSpace(csp, &ICS_ENORM, PICMptype) &&
			!PmatchColorSpace(csp, &ICS_VEC3, PICMptype))
		return NULL;

	return "32-bit encoded direction";
}

/* Write out a 32-bit encoded normal map */
static long
NormWriteImage(const char *fname, const ImgWriteBuf *wb)
{
	int		orient = PIXSTANDARD;
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
	if (!NormSupportedCS(wb->csp, -1))
		return 0;
						/* open output file */
	if ((pout = fopen(fname, "wb")) == NULL)
		return 0;
						/* write header + resolution */
	RHstartHeader(&wb->info, pout);
	RHoptViewAngle(wb, pout);
	fputformat((char *)NORMAL32FMT, pout);
	fputc('\n', pout);			/* end of info. header */
	if (wb->info.flags & IIForientation)
		orient = RHortab[wb->info.orientation-1];
	fputresolu(orient, wb->xres, wb->yres, pout);
						/* write each scanline */
	if (wb->csp->dtype == IDTint)
	    for (y = 0; y < wb->yres; y++) {
		const uint32 *	ip = (const uint32 *)(wb->img + y*wb->rowsize);
		for (x = wb->xres; x--; ip++)
			putint(*ip, 4, pout);
	    }
	else /* wb->csp->dtype == IDTfloat */
	    for (y = 0; y < wb->yres; y++) {
		const float *	vp = (const float *)(wb->img + y*wb->rowsize);
		for (x = wb->xres; x--; vp += 3) {
			FVECT	v;
			VCOPY(v, vp);
			normalize(v);
			if (putint(encodedir(v), 4, pout) == EOF) {
				fclose(pout);
				return 0;
			}
		}
	    }
	flen = ftell(pout);			/* return file length */
	return flen * (fclose(pout) != EOF);	/* or 0 on write error */
}

/* Norm writer interface */
const ImgWriterInterface	IWInterfaceNRM = {
	"nrm", &NormSupportedCS, &NormWriteImage
};
