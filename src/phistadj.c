/*
 *  phistadj.c
 *  panlib
 *
 *  Pancine histogram adjustment tone-mapping operator.
 *
 *  Created by Greg Ward on 1/5/22.
 *  Copyright 2021 Anyhere Software. All rights reserved.
 */

#include <stdio.h>
#include "dmessage.h"
#include "pimage.h"
#include "color.h"
#include "tonemap.h"

#ifndef true
#define true	1
#define false	0
#endif

#ifdef _WIN32
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif


/* check for and report a tone-mapping error */
static int
checkTMOerr(int erc, int lineno)
{
	DMsgClass	eclass;		/* translate to error class */

	switch (erc) {
	case TM_E_OK:			/* no error, so return false */
		return false;
	case TM_E_NOMEM:		/* out of memory */
		eclass = DMCmemory;
		break;
	case TM_E_BADFILE:		/* cannot open or understand file */
		eclass = DMCresource;
		break;
	case TM_E_TMFAIL:		/* cannot compute tone mapping */
		eclass = DMCdata;
		break;
	case TM_E_ILLEGAL:		/* illegal argument value */
	case TM_E_TMINVAL:		/* no valid tone mapping */
	case TM_E_CODERR1:		/* code consistency error 1 */
	case TM_E_CODERR2:		/* code consistency error 2 */
		eclass = DMCparameter;
		break;
	default:
		dmessage(DMCparameter, "Unknown TMO error!", __FILE__, lineno);
		return true;
	}
					/* report error and return true */
	dmessage(eclass, tmErrorMessage[erc], __FILE__, lineno);
	return true;
}

#define TMOerr(call)	checkTMOerr(call, __LINE__)

/* Compute histogram adjustment tone-mapping operator (camera) */
int
PhistAdjTMO(ImgStruct *ib, const ImgStruct *ia, double LdDyn)
{
	int		tflags = TM_F_CAMERA|TM_F_NOSTDERR;
	int		ok = false;
	TMbright	*limg = NULL;
	TMstruct	*tms = NULL;
	int		y;
						/* check arguments */
	if (!ib || !ib->csp)
		return false;
	if (!ia || !ia->img)
		return false;
	if ((ib->csp->dtype != IDTubyte) | ((ia->csp->dtype != IDTfloat) &
						(ia->csp->dtype != IDTushort))) {
		DMESG(DMCparameter, "Improper input/output color space in PhistAdjTMO");
		return false;
	}
	if (!ib->img) {				/* allocate destination? */
		ib->xres = ia->xres; ib->yres = ia->yres;
		if (!PnewImage(ib, .0))
			return false;
	} else if ((ib->xres != ia->xres) | (ib->yres != ia->yres)) {
		DMESG(DMCparameter, "Mismatched image dimensions");
		return false;
	}
	tflags |= TM_F_BW*((ib->csp->format == IPFy) | (ia->csp->format == IPFy));
	tms = tmInit(tflags, tflags&TM_F_BW ? stdprims : (RGBPRIMP)ib->csp->chroma,
			ib->csp->gamma);
	if (!tms) {
		DMESG(DMCparameter, "Cannot initialize tone-mapping");
		goto cleanup;
	}
	if (TMOerr(tmSetSpace(tms, tflags&TM_F_BW ? stdprims :
			ia->csp->format==IPFxyz ? TM_XYZPRIM :
					(RGBPRIMP)ia->csp->chroma, 1.f)) != TM_E_OK)
		goto cleanup;
						/* convert image data */
	limg = (TMbright *)malloc(sizeof(TMbright)*ia->xres*ia->yres);
	if (!limg) {
		DMESG(DMCmemory, "Cannot allocate temp image in PhistAdjTMO");
		goto cleanup;
	}
	switch ((ia->csp->dtype == IDTushort) << 2 |
			(ia->csp->format == IPFy) << 1 |
			(ib->csp->format == IPFy)) {
	case 0:					/* float RGB input, RGB output */
		for (y = ia->yres; y--; )
			if (TMOerr(tmCvColors(tms, limg + y*ia->xres, ProwPtr(ib,y),
						(COLOR *)ProwPtr(ia,y), ia->xres)))
				goto cleanup;
		break;
	case 04:				/* short RGB input, RGB output */
		for (y = ia->yres; y--; )
			if (TMOerr(tmCvRGB48(tms, limg + y*ia->xres, ProwPtr(ib,y),
						(uint16 (*)[3])ProwPtr(ia,y), ia->xres,
						ia->csp->gamma)))
				goto cleanup;
		break;
	case 01:				/* float RGB input, Y output */
		for (y = ia->yres; y--; )
			if (TMOerr(tmCvColors(tms, limg + y*ia->xres, TM_NOCHROM,
						(COLOR *)ProwPtr(ia,y), ia->xres)))
				goto cleanup;
		break;
	case 05:				/* short RGB input, Y output */
		for (y = ia->yres; y--; )
			if (TMOerr(tmCvRGB48(tms, limg + y*ia->xres, TM_NOCHROM,
						(uint16 (*)[3])ProwPtr(ia,y), ia->xres,
						ia->csp->gamma)))
				goto cleanup;
		break;
	case 02:				/* float Y input, RGB output */
	case 03:				/* float Y input, Y output */
		for (y = ia->yres; y--; )
			if (TMOerr(tmCvGrays(tms, limg + y*ia->xres,
						(float *)ProwPtr(ia,y), ia->xres)))
				goto cleanup;
		break;
	case 06:				/* short Y input, RGB output */
	case 07:				/* short Y input, Y output */
		for (y = ia->yres; y--; )
			if (TMOerr(tmCvGray16(tms, limg + y*ia->xres,
						(uint16 *)ProwPtr(ia,y), ia->xres,
						ia->csp->gamma)))
				goto cleanup;
		break;
	}
						/* get histogram */
	if (TMOerr(tmAddHisto(tms, limg, ia->xres*ia->yres, 1)))
		goto cleanup;
	if (!tms->histo) {			/* black image special case? */
		PclearImage(ib, NULL);
		ok = true;
		goto cleanup;
	}					/* else compute tone mapping */
	if (TMOerr(tmComputeMapping(tms, ib->csp->gamma, LdDyn, 0)))
		goto cleanup;

	for (y = ib->yres; y--; ) {		/* map pixels using global curve */
		uby8	*dp = ProwPtr(ib,y);
		if (TMOerr(tmMapPixels(tms, dp, limg + y*ib->xres,
					tflags&TM_F_BW ? TM_NOCHROM : dp, ib->xres)))
			goto cleanup;
		if (tflags&TM_F_BW && ib->csp->format == IPFrgb)
			PgetRGB24fromY8(dp, dp, ib->xres);
	}
	ok = true;				/* successful operation */
cleanup:
	if (!ok) PfreeImage(ib);		/* free dest. on failure */
	if (limg) free(limg);
	tmDone(tms);
	return ok;
}
