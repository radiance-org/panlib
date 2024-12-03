/*
 *  radheader.h
 *  pan
 *
 *  Definitions and declarations for routines handling Radiance header i/o
 *  Include after "imgio.h"
 *
 *  Created by Greg Ward on 7/27/19.
 *  Copyright 2019 Anyhere Software. All rights reserved.
 *
 */

#ifndef _RADHEADER_H_
#define _RADHEADER_H_

#include "rtio.h"

#ifdef __cplusplus
extern "C" {
#endif

extern const short	RHortab[8];	/* orientation conversion table */

/* Extract horizontal view angle from view parameter string */
extern int	RHgetAngle(ImgInfo *info);
/* Check if header line fits "VAR= value" format */
extern int	RHisParameter(const char *hl);
/* Scan header line, modifying ImgInfo struct appropriately */
extern int	RHscanHeadline(char *hl, void *p);

/* Get header information and bundle for calling application */
extern int	RHgetInfo(FILE *fp, ImgInfo *info);

/* Create Radiance information header */
extern void	RHstartHeader(const ImgInfo *info, FILE *pout);

/* Put out view angles if no view */
#define RHoptViewAngle(wb,pout) \
	if (((wb)->info.flags & (IIFview|IIFhvangle)) == IIFhvangle && \
			(.1 < (wb)->info.hvangle) & ((wb)->info.hvangle < 125.)) \
		fprintf(pout, "VIEW= -vtv -vh %f -vv %f\n", (wb)->info.hvangle, \
			360./M_PI*atan((wb)->yres/((wb)->pixAspect*(wb)->xres) \
				* tan(M_PI/360.*(wb)->info.hvangle) )); else

#ifdef __cplusplus
}
#endif

#endif	/* _RADHEADER_H_ */
