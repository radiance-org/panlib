/*
 *  piresamp.c
 *  panlib
 *
 *  Pancine image resampling routines.
 *
 *  Created by Greg Ward on 1/4/05.
 *  Modified by plonghurst on Thursday Feb 15 2007.
 *  Copyright 2005 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "system.h"
#include "dmessage.h"
#include "pimage.h"
#include "tiff.h"		/* needed for uint16, etc. */

/* pl: Fix for Windows */
#if defined(_WIN32) || defined(_WIN64)
#undef expf
#define expf(x)		(float)exp((double)(x))
#endif

#ifndef true
#define true	1
#define false	0
#endif

#define ifloor(x)	(((x) < 0) ? (int)((x)-.999999f) : (int)(x))

#define MAXINTERP	1.62		/* max. downsampling ratio for interp. */

#define OPTRAD		0.42f		/* optimal Gaussian radius */
#define	MAXRADM		3.2f		/* multiplier for search radius */
#define MINWT		0.002f		/* minimum useful weight */
#define MINIWT		20		/* integer version of MINWT */

#define comp_interp(hv,t,b)	((b)==PSlinear ? linear_interp(hv,t) : \
					cubic_interp(hv,t))
#define comp_interp8(hv,t,b)	((b)==PSlinear ? linear_interp8(hv,t) : \
					cubic_interp8(hv,t))

static const char	basisNam[][10] = {"Fast", "Best", "Nearest", "Box",
						"Gaussian", "Linear", "Cubic"};

/* Compute Catmull-Rom interpolant basis */
static void
cubic_interp(float hv[4], const float t)
{
	const float	t2 = t*t;
	const float	t3 = t2*t;

	hv[0] = .5f*(-t3 + 2.f*t2 - t);
	hv[1] = .5f*(3.f*t3 - 5.f*t2 + 2.f);
	hv[2] = .5f*(-3.f*t3 + 4.f*t2 + t);
	hv[3] = .5f*(t3 - t2);
}

/* Compute Catmull-Rom interpolant basis (0-1 mapped to 0-256) */
static void
cubic_interp8(short hv[4], const float t)
{
	const float	t2 = t*t;
	const float	t3 = t2*t;

	hv[0] = 128.f*(-t3 + 2.f*t2 - t) - .5f;
	hv[1] = 128.f*(3.f*t3 - 5.f*t2 + 2.f) + .5f;
	hv[2] = 128.f*(-3.f*t3 + 4.f*t2 + t) + .5f;
	hv[3] = 128.f*(t3 - t2) - .5f;
}

/* Linear interpolant basis */
static void
linear_interp(float hv[4], const float t)
{
	hv[0] = hv[3] = 0;
	hv[1] = 1.f - t;
	hv[2] = t;
}

/* Linear interpolant basis (0-1 mapped to 0-256) */
static void
linear_interp8(short hv[4], const float t)
{
	hv[0] = hv[3] = 0;
	hv[2] = 256.f*t + .5f;
	hv[1] = 256 - hv[2];
}

/* Nearest-neighbor resampling */
static int
nearestN(ImgStruct *dim, const ImgStruct *sim)
{
	const int	psiz = ImgPixelSize(sim->csp);
	int		xd, yd;
					/* resample image */
	sprintf(dmessage_buf, "Nearest neighbor resample of %dx%d to %dx%d image",
			sim->xres, sim->yres, dim->xres, dim->yres);
	DMESG(DMCtrace, dmessage_buf);
	for (yd = 0; yd < dim->yres; yd++) {
		const int	ys = (int)((yd+.5f)*sim->yres/dim->yres);
		const uby8	*srow = ProwPtr(sim,ys);
		uby8		*drow = ProwPtr(dim,yd);
					/* copy pixels along scanline */
		for (xd = 0; xd < dim->xres; xd++) {
			int	xs = (int)((xd+.5f)*sim->xres/dim->xres);
			int	i;
			
			xs *= psiz;
			for (i = psiz; i--; )
				*drow++ = srow[xs++];
		}
	}
	return true;
}

/* Box filter downsampling (only exact for even image size multiples) */
static int
boxDownsample(ImgStruct *dim, const ImgStruct *sim)
{
	const int		hsiz = sim->xres/dim->xres;
	const int		vsiz = sim->yres/dim->yres;
	const float		bwt = 1.f/(float)(hsiz*vsiz);
	const ImgColorSpace	*dst_csp = dim->csp;
	PweightMeasure		*wml = (PweightMeasure *)malloc((hsiz*vsiz+1) *
							sizeof(PweightMeasure));
	int			i, j;

	if (wml == NULL) {
		DMESG(DMCmemory, "malloc() failed in boxDownsample");
		return false;
	}
	sprintf(dmessage_buf, "Box downsample of %dx%d to %dx%d image",
			sim->xres, sim->yres, dim->xres, dim->yres);
	DMESG(DMCtrace, dmessage_buf);
	wml += hsiz*vsiz;
	wml->wt = 0;
	for (j = 0; j < vsiz; j++)	/* assign uniform weights */
		for (i = 0; i < hsiz; i++) {
			--wml;
			wml->wt = bwt;
			wml->mx = i - (hsiz>>1);
			wml->my = j - (vsiz>>1);
		}
					/* perform actual box filter */
	dim->csp = sim->csp;		/* avoid color conversion! */
	i = PsampleImage(dim, sim, &wml, 1, 1,
				.0, .0, NULL, NULL, PHCrgb, Pblack);
	dim->csp = dst_csp;
	free(wml);
	return i;
}

/* Resample a floating-point grayscale image using interpolation */
static int
interpY(ImgStruct *dim, const ImgStruct *sim, int basis)
{
	float	*rowv;
	float	*colwta, *cwp;
	int	xd, yd;
					/* allocate row vector holder */
	rowv = (float *)malloc(sizeof(float)*(sim->xres+4));
	colwta = (float *)malloc(sizeof(float)*4*dim->xres);
	if ((rowv == NULL) | (colwta == NULL)) {
		DMESG(DMCmemory, "malloc() failed in interpY");
		return false;
	}
	rowv += 2;
	cwp = colwta;			/* assign column weight array */
	for (xd = 0; xd < dim->xres; xd++) {
		const float	xsc = (xd+.5f)*sim->xres/dim->xres - .5f;
		comp_interp(cwp, xsc-ifloor(xsc), basis);
		cwp += 4;
	}
					/* interpolate image */
	sprintf(dmessage_buf, "%s interpolation of %dx%d to %dx%d float Y image",
			basisNam[basis], sim->xres, sim->yres, dim->xres, dim->yres);
	DMESG(DMCtrace, dmessage_buf);
	for (yd = 0; yd < dim->yres; yd++) {
		const float	ysc = (yd+.5f)*sim->yres/dim->yres - .5f;
		const int	ys1 = ifloor(ysc);
		const float	*srow1 = (ys1 <= 0) ? (const float *)sim->img :
					(const float *)ProwPtr(sim,ys1);
		const float	*srow0 = (ys1 <= 0) ? srow1 :
					(const float *)ProwPtr(sim,ys1-1);
		const float	*srow2 = (ys1 >= sim->yres-1) ? srow1 :
					(const float *)ProwPtr(sim,ys1+1);
		const float	*srow3 = (ys1 >= sim->yres-2) ? srow2 :
					(const float *)ProwPtr(sim,ys1+2);
		float		*drow = (float *)ProwPtr(dim,yd);
		float		hvec[4];
		int		xs;
					/* premultiply scanline at this yd */
		comp_interp(hvec, ysc-ys1, basis);

		for (xs = 0; xs < sim->xres; xs++) {
			rowv[xs] =	*srow0++ * hvec[0] +
					*srow1++ * hvec[1] +
					*srow2++ * hvec[2] +
					*srow3++ * hvec[3] ;
		}
		rowv[-2] = rowv[-1] = rowv[0];
		rowv[sim->xres+1] = rowv[sim->xres] = rowv[sim->xres-1];
					/* compute pixels along scanline */
		cwp = colwta;
		for (xd = 0; xd < dim->xres; xd++) {
			const float	xsc = (xd+.5f)*sim->xres/dim->xres - .5f;
			const int	xs1 = ifloor(xsc);
			
			*drow++ =	rowv[xs1-1] * cwp[0] +
					rowv[xs1  ] * cwp[1] +
					rowv[xs1+1] * cwp[2] + 
					rowv[xs1+2] * cwp[3] ;
			cwp += 4;
		}
	}
					/* clean up */
	free(rowv-2);
	free(colwta);
	return true;
}

/* Resample a floating-point RGB image using interpolation */
static int
interpRGB(ImgStruct *dim, const ImgStruct *sim, int basis)
{
	float	(*rowv)[3];
	float	*colwta, *cwp;
	int	xd, yd;
					/* allocate row vector holder */
	rowv = (float (*)[3])malloc(sizeof(float)*3*(sim->xres+4));
	colwta = (float *)malloc(sizeof(float)*4*dim->xres);
	if ((rowv == NULL) | (colwta == NULL)) {
		DMESG(DMCmemory, "malloc() failed in interpRGB");
		return false;
	}
	rowv += 2;
	cwp = colwta;			/* assign column weight array */
	for (xd = 0; xd < dim->xres; xd++) {
		const float	xsc = (xd+.5f)*sim->xres/dim->xres - .5f;
		comp_interp(cwp, xsc-ifloor(xsc), basis);
		cwp += 4;
	}
					/* interpolate image */
	sprintf(dmessage_buf, "%s interpolation of %dx%d to %dx%d float 3-channel image",
			basisNam[basis], sim->xres, sim->yres, dim->xres, dim->yres);
	DMESG(DMCtrace, dmessage_buf);
	for (yd = 0; yd < dim->yres; yd++) {
		const float	ysc = (yd+.5f)*sim->yres/dim->yres - .5f;
		const int	ys1 = ifloor(ysc);
		const float	*srow1 = (ys1 <= 0) ? (const float *)sim->img :
					(const float *)ProwPtr(sim,ys1);
		const float	*srow0 = (ys1 <= 0) ? srow1 :
					(const float *)ProwPtr(sim,ys1-1);
		const float	*srow2 = (ys1 >= sim->yres-1) ? srow1 :
					(const float *)ProwPtr(sim,ys1+1);
		const float	*srow3 = (ys1 >= sim->yres-2) ? srow2 :
					(const float *)ProwPtr(sim,ys1+2);
		float		*drow = (float *)ProwPtr(dim,yd);
		float		hvec[4];
		int		xs;
					/* premultiply scanline at this yd */
		comp_interp(hvec, ysc-ys1, basis);

		for (xs = 0; xs < sim->xres; xs++) {
			rowv[xs][0] =	*srow0++ * hvec[0] +
					*srow1++ * hvec[1] +
					*srow2++ * hvec[2] +
					*srow3++ * hvec[3] ;
			rowv[xs][1] =	*srow0++ * hvec[0] +
					*srow1++ * hvec[1] +
					*srow2++ * hvec[2] +
					*srow3++ * hvec[3] ;
			rowv[xs][2] =	*srow0++ * hvec[0] +
					*srow1++ * hvec[1] +
					*srow2++ * hvec[2] +
					*srow3++ * hvec[3] ;
		}
		rowv[-2][0] = rowv[-1][0] = rowv[0][0];
		rowv[-2][1] = rowv[-1][1] = rowv[0][1];
		rowv[-2][2] = rowv[-1][2] = rowv[0][2];
		rowv[sim->xres+1][0] = rowv[sim->xres][0] = rowv[sim->xres-1][0];
		rowv[sim->xres+1][1] = rowv[sim->xres][1] = rowv[sim->xres-1][1];
		rowv[sim->xres+1][2] = rowv[sim->xres][2] = rowv[sim->xres-1][2];
					/* compute pixels along scanline */
		cwp = colwta;
		for (xd = 0; xd < dim->xres; xd++) {
			const float	xsc = (xd+.5f)*sim->xres/dim->xres - .5f;
			const int	xs1 = ifloor(xsc);
			
			*drow++ =	rowv[xs1-1][0]*cwp[0] +
					rowv[xs1  ][0]*cwp[1] +
					rowv[xs1+1][0]*cwp[2] +
					rowv[xs1+2][0]*cwp[3] ;
			*drow++ =	rowv[xs1-1][1]*cwp[0] +
					rowv[xs1  ][1]*cwp[1] +
					rowv[xs1+1][1]*cwp[2] +
					rowv[xs1+2][1]*cwp[3] ;
			*drow++ =	rowv[xs1-1][2]*cwp[0] +
					rowv[xs1  ][2]*cwp[1] +
					rowv[xs1+1][2]*cwp[2] +
					rowv[xs1+2][2]*cwp[3] ;
			cwp += 4;
		}
	}
					/* clean up */
	free(rowv-2);
	free(colwta);
	return true;
}

/* Resample an 8-bit grayscale image using interpolation */
static int
interpY8(ImgStruct *dim, const ImgStruct *sim, int basis)
{
	int	*rowv;
	short	*colwta, *cwp;
	int	xd, yd;
					/* allocate row vector holder */
	rowv = (int *)malloc(sizeof(int)*(sim->xres+4));
	colwta = (short *)malloc(sizeof(short)*4*dim->xres);
	if ((rowv == NULL) | (colwta == NULL)) {
		DMESG(DMCmemory, "malloc() failed in interpY8");
		return false;
	}
	rowv += 2;
	cwp = colwta;			/* assign column weight array */
	for (xd = 0; xd < dim->xres; xd++) {
		const float	xsc = (xd+.5f)*sim->xres/dim->xres - .5f;
		comp_interp8(cwp, xsc-ifloor(xsc), basis);
		cwp += 4;
	}
					/* interpolate image */
	sprintf(dmessage_buf, "%s interpolation of %dx%d to %dx%d byte Y image",
			basisNam[basis], sim->xres, sim->yres, dim->xres, dim->yres);
	DMESG(DMCtrace, dmessage_buf);
	for (yd = 0; yd < dim->yres; yd++) {
		const float	ysc = (yd+.5f)*sim->yres/dim->yres - .5f;
		const int	ys1 = ifloor(ysc);
		const uby8	*srow1 = (ys1 <= 0) ? sim->img :
						ProwPtr(sim,ys1);
		const uby8	*srow0 = (ys1 <= 0) ? srow1 :
						ProwPtr(sim,ys1-1);
		const uby8	*srow2 = (ys1 >= sim->yres-1) ? srow1 :
						ProwPtr(sim,ys1+1);
		const uby8	*srow3 = (ys1 >= sim->yres-2) ? srow2 :
						ProwPtr(sim,ys1+2);
		uby8		*drow = ProwPtr(dim,yd);
		short		hvec[4];
		int		xs;
					/* premultiply scanline at this yd */
		comp_interp8(hvec, ysc-ys1, basis);

		for (xs = 0; xs < sim->xres; xs++) {
			rowv[xs] =	(int)*srow0++ * hvec[0] +
					(int)*srow1++ * hvec[1] +
					(int)*srow2++ * hvec[2] +
					(int)*srow3++ * hvec[3] ;
		}
		rowv[-2] = rowv[-1] = rowv[0];
		rowv[sim->xres+1] = rowv[sim->xres] = rowv[sim->xres-1];
					/* compute pixels along scanline */
		cwp = colwta;
		for (xd = 0; xd < dim->xres; xd++) {
			const float	xsc = (xd+.5f)*sim->xres/dim->xres - .5f;
			const int	xs1 = ifloor(xsc);
			int		res;
			
			res =	rowv[xs1-1]*cwp[0] + rowv[xs1  ]*cwp[1] +
				rowv[xs1+1]*cwp[2] + rowv[xs1+2]*cwp[3] ;

			*drow++ = (res <= 0) ? 0 :
					(res >= 0xff0000) ? 0xff :
					(res + (1<<15)) >> 16;
			cwp += 4;
		}
	}
					/* clean up */
	free(rowv-2);
	free(colwta);
	return true;
}

/* Resample a 24-bit RGB image using interpolation */
static int
interpRGB24(ImgStruct *dim, const ImgStruct *sim, int basis)
{
	int	(*rowv)[3];
	short	*colwta, *cwp;
	int	xd, yd;
					/* allocate row vector holder */
	rowv = (int (*)[3])malloc(sizeof(int)*3*(sim->xres+4));
	colwta = (short *)malloc(sizeof(short)*4*dim->xres);
	if ((rowv == NULL) | (colwta == NULL)) {
		DMESG(DMCmemory, "malloc() failed in interpRGB24");
		return false;
	}
	rowv += 2;
	cwp = colwta;			/* assign column weight array */
	for (xd = 0; xd < dim->xres; xd++) {
		const float	xsc = (xd+.5f)*sim->xres/dim->xres - .5f;
		comp_interp8(cwp, xsc-ifloor(xsc), basis);
		cwp += 4;
	}
					/* interpolate image */
	sprintf(dmessage_buf, "%s interpolation of %dx%d to %dx%d byte 3-channel image",
			basisNam[basis], sim->xres, sim->yres, dim->xres, dim->yres);
	DMESG(DMCtrace, dmessage_buf);
	for (yd = 0; yd < dim->yres; yd++) {
		const float	ysc = (yd+.5f)*sim->yres/dim->yres - .5f;
		const int	ys1 = ifloor(ysc);
		const uby8	*srow1 = (ys1 <= 0) ? sim->img :
						ProwPtr(sim,ys1);
		const uby8	*srow0 = (ys1 <= 0) ? srow1 :
						ProwPtr(sim,ys1-1);
		const uby8	*srow2 = (ys1 >= sim->yres-1) ? srow1 :
						ProwPtr(sim,ys1+1);
		const uby8	*srow3 = (ys1 >= sim->yres-2) ? srow2 :
						ProwPtr(sim,ys1+2);
		uby8		*drow = ProwPtr(dim,yd);
		short		hvec[4];
		int		xs;
					/* premultiply scanline at this yd */
		comp_interp8(hvec, ysc-ys1, basis);

		for (xs = 0; xs < sim->xres; xs++) {
			rowv[xs][0] =	(int)*srow0++ * hvec[0] +
					(int)*srow1++ * hvec[1] +
					(int)*srow2++ * hvec[2] +
					(int)*srow3++ * hvec[3] ;
			rowv[xs][1] =	(int)*srow0++ * hvec[0] +
					(int)*srow1++ * hvec[1] +
					(int)*srow2++ * hvec[2] +
					(int)*srow3++ * hvec[3] ;
			rowv[xs][2] =	(int)*srow0++ * hvec[0] +
					(int)*srow1++ * hvec[1] +
					(int)*srow2++ * hvec[2] +
					(int)*srow3++ * hvec[3] ;
		}
		rowv[-2][0] = rowv[-1][0] = rowv[0][0];
		rowv[-2][1] = rowv[-1][1] = rowv[0][1];
		rowv[-2][2] = rowv[-1][2] = rowv[0][2];
		rowv[sim->xres+1][0] = rowv[sim->xres][0] = rowv[sim->xres-1][0];
		rowv[sim->xres+1][1] = rowv[sim->xres][1] = rowv[sim->xres-1][1];
		rowv[sim->xres+1][2] = rowv[sim->xres][2] = rowv[sim->xres-1][2];
					/* compute pixels along scanline */
		cwp = colwta;
		for (xd = 0; xd < dim->xres; xd++) {
			const float	xsc = (xd+.5f)*sim->xres/dim->xres - .5f;
			const int	xs1 = ifloor(xsc);
			int		res;

			res =	rowv[xs1-1][0]*cwp[0] + rowv[xs1  ][0]*cwp[1] +
				rowv[xs1+1][0]*cwp[2] + rowv[xs1+2][0]*cwp[3] ;
			*drow++ = (res <= 0) ? 0 :
					(res >= 0xff0000) ? 0xff :
					(res + (1<<15)) >> 16;
			res =	rowv[xs1-1][1]*cwp[0] + rowv[xs1  ][1]*cwp[1] +
				rowv[xs1+1][1]*cwp[2] + rowv[xs1+2][1]*cwp[3] ;
			*drow++ = (res <= 0) ? 0 :
					(res >= 0xff0000) ? 0xff :
					(res + (1<<15)) >> 16;
			res =	rowv[xs1-1][2]*cwp[0] + rowv[xs1  ][2]*cwp[1] +
				rowv[xs1+1][2]*cwp[2] + rowv[xs1+2][2]*cwp[3] ;
			*drow++ = (res <= 0) ? 0 :
					(res >= 0xff0000) ? 0xff :
					(res + (1<<15)) >> 16;
			cwp += 4;
		}
	}
					/* clean up */
	free(rowv-2);
	free(colwta);
	return true;
}

/* Resample a 16-bit grayscale image using interpolation */
static int
interpY16(ImgStruct *dim, const ImgStruct *sim, int basis)
{
	int	*rowv;
	short	*colwta, *cwp;
	int	xd, yd;
					/* allocate row vector holder */
	rowv = (int *)malloc(sizeof(int)*(sim->xres+4));
	colwta = (short *)malloc(sizeof(short)*4*dim->xres);
	if ((rowv == NULL) | (colwta == NULL)) {
		DMESG(DMCmemory, "malloc() failed in interpY16");
		return false;
	}
	rowv += 2;
	cwp = colwta;			/* assign column weight array */
	for (xd = 0; xd < dim->xres; xd++) {
		const float	xsc = (xd+.5f)*sim->xres/dim->xres - .5f;
		comp_interp8(cwp, xsc-ifloor(xsc), basis);
		cwp += 4;
	}
					/* interpolate image */
	sprintf(dmessage_buf, "%s interpolation of %dx%d to %dx%d 16-bit Y image",
			basisNam[basis], sim->xres, sim->yres, dim->xres, dim->yres);
	DMESG(DMCtrace, dmessage_buf);
	for (yd = 0; yd < dim->yres; yd++) {
		const float	ysc = (yd+.5f)*sim->yres/dim->yres - .5f;
		const int	ys1 = ifloor(ysc);
		const uint16	*srow1 = (ys1 <= 0) ? (const uint16 *)sim->img :
					(const uint16 *)ProwPtr(sim,ys1);
		const uint16	*srow0 = (ys1 <= 0) ? srow1 :
					(const uint16 *)ProwPtr(sim,ys1-1);
		const uint16	*srow2 = (ys1 >= sim->yres-1) ? srow1 :
					(const uint16 *)ProwPtr(sim,ys1+1);
		const uint16	*srow3 = (ys1 >= sim->yres-2) ? srow2 :
					(const uint16 *)ProwPtr(sim,ys1+2);
		uint16		*drow = (uint16 *)ProwPtr(dim,yd);
		short		hvec[4];
		int		xs;
					/* premultiply scanline at this yd */
		comp_interp8(hvec, ysc-ys1, basis);

		for (xs = 0; xs < sim->xres; xs++) {
			rowv[xs] = (	(int)*srow0++ * hvec[0] +
					(int)*srow1++ * hvec[1] +
					(int)*srow2++ * hvec[2] +
					(int)*srow3++ * hvec[3] + (1<<7) ) >> 8;
		}
		rowv[-2] = rowv[-1] = rowv[0];
		rowv[sim->xres+1] = rowv[sim->xres] = rowv[sim->xres-1];
					/* compute pixels along scanline */
		cwp = colwta;
		for (xd = 0; xd < dim->xres; xd++) {
			const float	xsc = (xd+.5f)*sim->xres/dim->xres - .5f;
			const int	xs1 = ifloor(xsc);
			int		res;
			
			res =	rowv[xs1-1]*cwp[0] + rowv[xs1  ]*cwp[1] +
				rowv[xs1+1]*cwp[2] + rowv[xs1+2]*cwp[3] ;

			*drow++ = (res <= 0) ? 0 :
					(res >= 0xffff00) ? 0xffff :
					(res + (1<<7)) >> 8;
			cwp += 4;
		}
	}
					/* clean up */
	free(rowv-2);
	free(colwta);
	return true;
}

/* Resample a 48-bit RGB image using interpolation */
static int
interpRGB48(ImgStruct *dim, const ImgStruct *sim, int basis)
{
	int	(*rowv)[3];
	short	*colwta, *cwp;
	int	xd, yd;
					/* allocate row vector holder */
	rowv = (int (*)[3])malloc(sizeof(int)*3*(sim->xres+4));
	colwta = (short *)malloc(sizeof(short)*4*dim->xres);
	if ((rowv == NULL) | (colwta == NULL)) {
		DMESG(DMCmemory, "malloc() failed in interpRGB48");
		return false;
	}
	rowv += 2;
	cwp = colwta;			/* assign column weight array */
	for (xd = 0; xd < dim->xres; xd++) {
		const float	xsc = (xd+.5f)*sim->xres/dim->xres - .5f;
		comp_interp8(cwp, xsc-ifloor(xsc), basis);
		cwp += 4;
	}
					/* interpolate image */
	sprintf(dmessage_buf, "%s interpolation of %dx%d to %dx%d 16-bit by 3-channel image",
			basisNam[basis], sim->xres, sim->yres, dim->xres, dim->yres);
	DMESG(DMCtrace, dmessage_buf);
	for (yd = 0; yd < dim->yres; yd++) {
		const float	ysc = (yd+.5f)*sim->yres/dim->yres - .5f;
		const int	ys1 = ifloor(ysc);
		const uint16	*srow1 = (ys1 <= 0) ? (const uint16 *)sim->img :
					(const uint16 *)ProwPtr(sim,ys1);
		const uint16	*srow0 = (ys1 <= 0) ? srow1 :
					(const uint16 *)ProwPtr(sim,ys1-1);
		const uint16	*srow2 = (ys1 >= sim->yres-1) ? srow1 :
					(const uint16 *)ProwPtr(sim,ys1+1);
		const uint16	*srow3 = (ys1 >= sim->yres-2) ? srow2 :
					(const uint16 *)ProwPtr(sim,ys1+2);
		uint16		*drow = (uint16 *)ProwPtr(dim,yd);
		short		hvec[4];
		int		xs;
					/* premultiply scanline at this yd */
		comp_interp8(hvec, ysc-ys1, basis);

		for (xs = 0; xs < sim->xres; xs++) {
			rowv[xs][0] = (	(int)*srow0++ * hvec[0] +
					(int)*srow1++ * hvec[1] +
					(int)*srow2++ * hvec[2] +
					(int)*srow3++ * hvec[3] + (1<<7) ) >> 8;
			rowv[xs][1] = (	(int)*srow0++ * hvec[0] +
					(int)*srow1++ * hvec[1] +
					(int)*srow2++ * hvec[2] +
					(int)*srow3++ * hvec[3] + (1<<7) ) >> 8;
			rowv[xs][2] = (	(int)*srow0++ * hvec[0] +
					(int)*srow1++ * hvec[1] +
					(int)*srow2++ * hvec[2] +
					(int)*srow3++ * hvec[3] + (1<<7) ) >> 8;
		}
		rowv[-2][0] = rowv[-1][0] = rowv[0][0];
		rowv[-2][1] = rowv[-1][1] = rowv[0][1];
		rowv[-2][2] = rowv[-1][2] = rowv[0][2];
		rowv[sim->xres+1][0] = rowv[sim->xres][0] = rowv[sim->xres-1][0];
		rowv[sim->xres+1][1] = rowv[sim->xres][1] = rowv[sim->xres-1][1];
		rowv[sim->xres+1][2] = rowv[sim->xres][2] = rowv[sim->xres-1][2];
					/* compute pixels along scanline */
		cwp = colwta;
		for (xd = 0; xd < dim->xres; xd++) {
			const float	xsc = (xd+.5f)*sim->xres/dim->xres - .5f;
			const int	xs1 = ifloor(xsc);
			int		res;
			
			res =	rowv[xs1-1][0]*cwp[0] + rowv[xs1  ][0]*cwp[1] +
				rowv[xs1+1][0]*cwp[2] + rowv[xs1+2][0]*cwp[3] ;
			*drow++ = (res <= 0) ? 0 :
					(res >= 0xffff00) ? 0xffff :
					(res + (1<<7)) >> 8;
			res =	rowv[xs1-1][1]*cwp[0] + rowv[xs1  ][1]*cwp[1] +
				rowv[xs1+1][1]*cwp[2] + rowv[xs1+2][1]*cwp[3] ;
			*drow++ = (res <= 0) ? 0 :
					(res >= 0xffff00) ? 0xffff :
					(res + (1<<7)) >> 8;
			res =	rowv[xs1-1][2]*cwp[0] + rowv[xs1  ][2]*cwp[1] +
				rowv[xs1+1][2]*cwp[2] + rowv[xs1+2][2]*cwp[3] ;
			*drow++ = (res <= 0) ? 0 :
					(res >= 0xffff00) ? 0xffff :
					(res + (1<<7)) >> 8;
			cwp += 4;
		}
	}
					/* clean up */
	free(rowv-2);
	free(colwta);
	return true;
}

/* Resample a floating point grayscale image using Gaussian average */
static int
gaussianY(ImgStruct *dim, const ImgStruct *sim, float rad)
{
	const float	xdens = (float)dim->xres/sim->xres;
	const float	ydens = (float)dim->yres/sim->yres;
	const int	hrad = (int)(MAXRADM*rad/xdens + .5f);
	const int	vrad = (int)(MAXRADM*rad/ydens + .5f);
	float		*rowv;
	float		*colwta, *cwp;
	int		xd, yd;
	int		j;
	/* pl: The following needed to be done for legal compilation in Windows */
	float *slwt = (float *)malloc(sizeof(float) * (2*vrad+1));

					/* allocate row vector holder */
	rowv = (float *)malloc(sizeof(float)*(sim->xres + 2*hrad));
	colwta = (float *)malloc(sizeof(float)*(2*hrad+1)*dim->xres);
	if ((slwt == NULL) | (rowv == NULL) | (colwta == NULL)) {
		DMESG(DMCmemory, "malloc() failed in gaussianY");
		return false;
	}
	rowv += hrad;
	cwp = colwta;			/* assign column weight array */
	for (xd = 0; xd < dim->xres; xd++) {
		const float	xsc = (xd+.5f)/xdens;
		const int	xs0 = (int)xsc;
		for (j = -hrad; j <= hrad; j++) {
			const float	rx = (xs0+j+.5f - xsc)*xdens;
			*cwp++ = expf(-rx*rx/(2.f*rad*rad));
		}
	}
					/* filter our image */
	sprintf(dmessage_buf, "Gaussian resample of %dx%d to %dx%d float Y image",
			sim->xres, sim->yres, dim->xres, dim->yres);
	DMESG(DMCtrace, dmessage_buf);
	for (yd = 0; yd < dim->yres; yd++) {
		const float	ysc = (yd+.5f)/ydens;
		const int	ys0 = (int)ysc;
		float		*drow = (float *)ProwPtr(dim,yd);
		float		wsum;
		int		xs;
					/* compute row vector at this yd */
		wsum = 0;
		for (j = -vrad; j <= vrad; j++) {
			const float	ry = (ys0+j+.5f - ysc)*ydens;
			slwt[j+vrad] = expf(-ry*ry/(2.f*rad*rad));
			if (slwt[j+vrad] < MINWT)
				slwt[j+vrad] = 0;
			else
				wsum += slwt[j+vrad];
		}
		memset(rowv, 0, sizeof(float)*sim->xres);
		for (j = -vrad; j <= vrad; j++) {
			const float	wt = slwt[j+vrad] / wsum;
			const float	*srow;
			if (wt <= 1e-8f)
				continue;
			srow = (const float *)ProwPtr(sim, (ys0+j < 0) ? 0 :
					(ys0+j >= sim->yres) ? sim->yres-1 : ys0+j);
			for (xs = 0; xs < sim->xres; xs++)
				rowv[xs] += *srow++ * wt;
		}
		for (xs = -hrad; xs < 0; xs++)
			rowv[xs] = rowv[0];
		for (xs = sim->xres; xs < sim->xres+hrad; xs++)
			rowv[xs] = rowv[sim->xres-1];
					/* compute pixels along scanline */
		cwp = colwta;
		for (xd = 0; xd < dim->xres; xd++) {
			const int	xs0 = (int)((xd+.5f)/xdens);
			*drow = 0;
			wsum = 0;
			for (j = -hrad; j <= hrad; j++) {
				*drow += rowv[xs0+j] * *cwp;
				wsum += *cwp++;
			}
			*drow++ /= wsum;
		}
	}
					/* clean up */
	free(slwt);
	free(rowv-hrad);
	free(colwta);
	return true;
}

/* Resample a floating point RGB image using Gaussian average */
static int
gaussianRGB(ImgStruct *dim, const ImgStruct *sim, float rad)
{
	const float	xdens = (float)dim->xres/sim->xres;
	const float	ydens = (float)dim->yres/sim->yres;
	const int	hrad = (int)(MAXRADM*rad/xdens + .5f);
	const int	vrad = (int)(MAXRADM*rad/ydens + .5f);
	float		(*rowv)[3];
	float		*colwta, *cwp;
	int		xd, yd;
	int		j;
	/* pl: The following needed to be done for legal compilation in Windows */
	float *slwt = (float *)malloc(sizeof(float) * (2*vrad+1));
	
					/* allocate row vector holder */
	rowv = (float (*)[3])malloc(sizeof(float)*3*(sim->xres + 2*hrad));
	colwta = (float *)malloc(sizeof(float)*(2*hrad+1)*dim->xres);
	if ((slwt == NULL) | (rowv == NULL) | (colwta == NULL)) {
		DMESG(DMCmemory, "malloc() failed in gaussianRGB");
		return false;
	}
	rowv += hrad;
	cwp = colwta;			/* assign column weight array */
	for (xd = 0; xd < dim->xres; xd++) {
		const float	xsc = (xd+.5f)/xdens;
		const int	xs0 = (int)xsc;
		for (j = -hrad; j <= hrad; j++) {
			const float	rx = (xs0+j+.5f - xsc)*xdens;
			*cwp++ = expf(-rx*rx/(2.f*rad*rad));
		}
	}
					/* filter our image */
	sprintf(dmessage_buf, "Gaussian resample of %dx%d to %dx%d float 3-channel image",
			sim->xres, sim->yres, dim->xres, dim->yres);
	DMESG(DMCtrace, dmessage_buf);
	for (yd = 0; yd < dim->yres; yd++) {
		const float	ysc = (yd+.5f)/ydens;
		const int	ys0 = (int)ysc;
		float		*drow = (float *)ProwPtr(dim,yd);
		float		wsum;
		int		xs;
					/* compute row vector at this yd */
		wsum = 0;
		for (j = -vrad; j <= vrad; j++) {
			const float	ry = (ys0+j+.5f - ysc)*ydens;
			slwt[j+vrad] = expf(-ry*ry/(2.f*rad*rad));
			if (slwt[j+vrad] < MINWT)
				slwt[j+vrad] = 0;
			else
				wsum += slwt[j+vrad];
		}
		memset(rowv, 0, sizeof(float)*3*sim->xres);
		for (j = -vrad; j <= vrad; j++) {
			const float	wt = slwt[j+vrad] / wsum;
			const float	*srow;
			if (wt <= 1e-8f)
				continue;
			srow = (const float *)ProwPtr(sim, (ys0+j < 0) ? 0 :
					(ys0+j >= sim->yres) ? sim->yres-1 : ys0+j);
			for (xs = 0; xs < sim->xres; xs++) {
				rowv[xs][0] += *srow++ * wt;
				rowv[xs][1] += *srow++ * wt;
				rowv[xs][2] += *srow++ * wt;
			}
		}
		for (xs = -hrad; xs < 0; xs++) {
			rowv[xs][0] = rowv[0][0];
			rowv[xs][1] = rowv[0][1];
			rowv[xs][2] = rowv[0][2];
		}
		for (xs = sim->xres; xs < sim->xres+hrad; xs++) {
			rowv[xs][0] = rowv[sim->xres-1][0];
			rowv[xs][1] = rowv[sim->xres-1][1];
			rowv[xs][2] = rowv[sim->xres-1][2];
		}
					/* compute pixels along scanline */
		cwp = colwta;
		for (xd = 0; xd < dim->xres; xd++) {
			const int	xs0 = (int)((xd+.5f)/xdens);
			wsum = 0;
			drow[0] = drow[1] = drow[2] = 0;
			for (j = -hrad; j <= hrad; j++) {
				drow[0] += rowv[xs0+j][0] * *cwp;
				drow[1] += rowv[xs0+j][1] * *cwp;
				drow[2] += rowv[xs0+j][2] * *cwp;
				wsum += *cwp++;
			}
			wsum = 1.f/wsum;
			*drow++ *= wsum;
			*drow++ *= wsum;
			*drow++ *= wsum;
		}
	}
					/* clean up */
	free(slwt);
	free(rowv-hrad);
	free(colwta);
	return true;
}

/* Resample an 8-bit grayscale image using Gaussian average */
static int
gaussianY8(ImgStruct *dim, const ImgStruct *sim, float rad)
{
	const float	xdens = (float)dim->xres/sim->xres;
	const float	ydens = (float)dim->yres/sim->yres;
	const int	hrad = (int)(MAXRADM*rad/xdens + .5f);
	const int	vrad = (int)(MAXRADM*rad/ydens + .5f);
	const float	hwsca = 0x7fff / (float)hrad;
	const float	vwsca = 0x7fff / (float)vrad;
	uint32		*rowv;
	short		*colwta, *cwp;
	int		xd, yd;
	int		j;
					/* allocate row vector holder */
	rowv = (uint32 *)malloc(sizeof(uint32)*(sim->xres + 2*hrad));
	colwta = (short *)malloc(sizeof(short)*(2*hrad+1)*dim->xres);
	if ((rowv == NULL) | (colwta == NULL)) {
		DMESG(DMCmemory, "malloc() failed in gaussianY8");
		return false;
	}
	rowv += hrad;
	cwp = colwta;			/* assign column weight array */
	for (xd = 0; xd < dim->xres; xd++) {
		const float	xsc = (xd+.5f)/xdens;
		const int	xs0 = (int)xsc;
		for (j = -hrad; j <= hrad; j++) {
			const float	rx = (xs0+j+.5f - xsc)*xdens;
			*cwp++ = (int)(hwsca*expf(-rx*rx/(2.f*rad*rad)) + .5f);
		}
	}
					/* filter our image */
	sprintf(dmessage_buf, "Gaussian resample of %dx%d to %dx%d byte Y image",
			sim->xres, sim->yres, dim->xres, dim->yres);
	DMESG(DMCtrace, dmessage_buf);
	for (yd = 0; yd < dim->yres; yd++) {
		const float	ysc = (yd+.5f)/ydens;
		const int	ys0 = (int)ysc;
		uby8		*drow = ProwPtr(dim,yd);
		uint32		wsum;
		int		xs, j;
					/* compute row vector at this yd */
		memset(rowv, 0, sizeof(uint32)*sim->xres);
		wsum = 0;
		for (j = -vrad; j <= vrad; j++) {
			const float	ry = (ys0+j+.5f - ysc)*ydens;
			const uint32	wt = (int)(vwsca*expf(-ry*ry/(2.f*rad*rad)) + .5f);
			const uby8	*srow;
			if (wt < MINIWT)
				continue;
			srow = ProwPtr(sim, (ys0+j < 0) ? 0 :
					(ys0+j >= sim->yres) ? sim->yres-1 : ys0+j);
			for (xs = 0; xs < sim->xres; xs++)
				rowv[xs] += *srow++ * wt;
			wsum += wt;
		}
		for (xs = sim->xres; xs--; )
			rowv[xs] = (rowv[xs] << 8) / wsum;
		for (xs = -hrad; xs < 0; xs++)
			rowv[xs] = rowv[0];
		for (xs = sim->xres; xs < sim->xres+hrad; xs++)
			rowv[xs] = rowv[sim->xres-1];
					/* compute pixels along scanline */
		cwp = colwta;
		for (xd = 0; xd < dim->xres; xd++) {
			const int	xs0 = (int)((xd+.5f)/xdens);
			uint32	res = 0;
			wsum = 0;
			for (j = -hrad; j <= hrad; j++) {
				res += rowv[xs0+j] * *cwp;
				wsum += *cwp++;
			}
			res /= wsum;
			*drow++ = (res >= 0xff00) ? 0xff : res >> 8;
		}
	}
					/* clean up */
	free(rowv-hrad);
	free(colwta);
	return true;
}

/* Resample a 24-bit RGB image using Gaussian average */
static int
gaussianRGB24(ImgStruct *dim, const ImgStruct *sim, float rad)
{
	const float	xdens = (float)dim->xres/sim->xres;
	const float	ydens = (float)dim->yres/sim->yres;
	const int	hrad = (int)(MAXRADM*rad/xdens + .5f);
	const int	vrad = (int)(MAXRADM*rad/ydens + .5f);
	const float	hwsca = 0x7fff / (float)hrad;
	const float	vwsca = 0x7fff / (float)vrad;
	uint32		(*rowv)[3];
	short		*colwta, *cwp;
	int		xd, yd;
	int		j;
					/* allocate row vector holder */
	rowv = (uint32 (*)[3])malloc(sizeof(uint32)*3*(sim->xres + 2*hrad));
	colwta = (short *)malloc(sizeof(short)*(2*hrad+1)*dim->xres);
	if ((rowv == NULL) | (colwta == NULL)) {
		DMESG(DMCmemory, "malloc() failed in gaussianRGB24");
		return false;
	}
	rowv += hrad;
	cwp = colwta;			/* assign column weight array */
	for (xd = 0; xd < dim->xres; xd++) {
		const float	xsc = (xd+.5f)/xdens;
		const int	xs0 = (int)xsc;
		for (j = -hrad; j <= hrad; j++) {
			const float	rx = (xs0+j+.5f - xsc)*xdens;
			*cwp++ = (int)(hwsca*expf(-rx*rx/(2.f*rad*rad)) + .5f);
		}
	}
					/* filter our image */
	sprintf(dmessage_buf, "Gaussian resample of %dx%d to %dx%d byte 3-channel image",
			sim->xres, sim->yres, dim->xres, dim->yres);
	DMESG(DMCtrace, dmessage_buf);
	for (yd = 0; yd < dim->yres; yd++) {
		const float	ysc = (yd+.5f)/ydens;
		const int	ys0 = (int)ysc;
		uby8		*drow = ProwPtr(dim,yd);
		uint32		wsum;
		int		xs;
					/* compute row vector at this yd */
		memset(rowv, 0, sizeof(uint32)*3*sim->xres);
		wsum = 0;
		for (j = -vrad; j <= vrad; j++) {
			const float	ry = (ys0+j+.5f - ysc)*ydens;
			const uint32	wt = (int)(vwsca*expf(-ry*ry/(2.f*rad*rad)) + .5f);
			const uby8	*srow;
			if (wt < MINIWT)
				continue;
			srow = ProwPtr(sim, (ys0+j < 0) ? 0 :
					(ys0+j >= sim->yres) ? sim->yres-1 : ys0+j);
			for (xs = 0; xs < sim->xres; xs++) {
				rowv[xs][0] += *srow++ * wt;
				rowv[xs][1] += *srow++ * wt;
				rowv[xs][2] += *srow++ * wt;
			}
			wsum += wt;
		}
		for (xs = sim->xres; xs--; ) {
			rowv[xs][0] = (rowv[xs][0] << 8) / wsum;
			rowv[xs][1] = (rowv[xs][1] << 8) / wsum;
			rowv[xs][2] = (rowv[xs][2] << 8) / wsum;
		}
		for (xs = -hrad; xs < 0; xs++) {
			rowv[xs][0] = rowv[0][0];
			rowv[xs][1] = rowv[0][1];
			rowv[xs][2] = rowv[0][2];
		}
		for (xs = sim->xres; xs < sim->xres+hrad; xs++) {
			rowv[xs][0] = rowv[sim->xres-1][0];
			rowv[xs][1] = rowv[sim->xres-1][1];
			rowv[xs][2] = rowv[sim->xres-1][2];
		}
					/* compute pixels along scanline */
		cwp = colwta;
		for (xd = 0; xd < dim->xres; xd++) {
			const int	xs0 = (int)((xd+.5f)/xdens);
			uint32		res[3];
			wsum = 0;
			res[0] = res[1] = res[2] = 0;
			for (j = -hrad; j <= hrad; j++) {
				res[0] += rowv[xs0+j][0] * *cwp;
				res[1] += rowv[xs0+j][1] * *cwp;
				res[2] += rowv[xs0+j][2] * *cwp;
				wsum += *cwp++;
			}
			res[0] /= wsum;
			res[1] /= wsum;
			res[2] /= wsum;
			*drow++ = (res[0] >= 0xff00) ? 0xff : res[0] >> 8;
			*drow++ = (res[1] >= 0xff00) ? 0xff : res[1] >> 8;
			*drow++ = (res[2] >= 0xff00) ? 0xff : res[2] >> 8;
		}
	}
					/* clean up */
	free(rowv-hrad);
	free(colwta);
	return true;
}

/* Resample a 16-bit grayscale image using Gaussian average */
static int
gaussianY16(ImgStruct *dim, const ImgStruct *sim, float rad)
{
	const float	xdens = (float)dim->xres/sim->xres;
	const float	ydens = (float)dim->yres/sim->yres;
	const int	hrad = (int)(MAXRADM*rad/xdens + .5f);
	const int	vrad = (int)(MAXRADM*rad/ydens + .5f);
	const float	hwsca = 0x7fff / (float)hrad;
	const float	vwsca = 0x7fff / (float)vrad;
	uint32		*rowv;
	uint16		*colwta, *cwp;
	int		xd, yd;
	int		j;
					/* allocate row vector holder */
	rowv = (uint32 *)malloc(sizeof(uint32)*(sim->xres + 2*hrad));
	colwta = (uint16 *)malloc(sizeof(uint16)*(2*hrad+1)*dim->xres);
	if ((rowv == NULL) | (colwta == NULL)) {
		DMESG(DMCmemory, "malloc() failed in gaussianY16");
		return false;
	}
	rowv += hrad;
	cwp = colwta;			/* assign column weight array */
	for (xd = 0; xd < dim->xres; xd++) {
		const float	xsc = (xd+.5f)/xdens;
		const int	xs0 = (int)xsc;
		for (j = -hrad; j <= hrad; j++) {
			const float	rx = (xs0+j+.5f - xsc)*xdens;
			*cwp++ = (int)(hwsca*expf(-rx*rx/(2.f*rad*rad)) + .5f);
		}
	}
					/* filter our image */
	sprintf(dmessage_buf, "Gaussian resample of %dx%d to %dx%d 16-bit Y image",
			sim->xres, sim->yres, dim->xres, dim->yres);
	DMESG(DMCtrace, dmessage_buf);
	for (yd = 0; yd < dim->yres; yd++) {
		const float	ysc = (yd+.5f)/ydens;
		const int	ys0 = (int)ysc;
		uint16		*drow = (uint16 *)ProwPtr(dim,yd);
		uint32		wsum;
		int		xs;
					/* compute row vector at this yd */
		memset(rowv, 0, sizeof(uint32)*sim->xres);
		wsum = 0;
		for (j = -vrad; j <= vrad; j++) {
			const float	ry = (ys0+j+.5f - ysc)*ydens;
			const uint32	wt = (int)(vwsca*expf(-ry*ry/(2.f*rad*rad)) + .5f);
			const uint16	*srow;
			if (wt < MINIWT)
				continue;
			srow = (const uint16 *)ProwPtr(sim, (ys0+j < 0) ? 0 :
					(ys0+j >= sim->yres) ? sim->yres-1 : ys0+j);
			for (xs = 0; xs < sim->xres; xs++)
				rowv[xs] += *srow++ * wt;
			wsum += wt;
		}
		for (xs = sim->xres; xs--; )
			rowv[xs] /= wsum;
		for (xs = -hrad; xs < 0; xs++)
			rowv[xs] = rowv[0];
		for (xs = sim->xres; xs < sim->xres+hrad; xs++)
			rowv[xs] = rowv[sim->xres-1];
					/* compute pixels along scanline */
		cwp = colwta;
		for (xd = 0; xd < dim->xres; xd++) {
			const int	xs0 = (int)((xd+.5f)/xdens);
			uint32		res = 0;
			wsum = 0;
			for (j = -hrad; j <= hrad; j++) {
				res += rowv[xs0+j] * *cwp;
				wsum += *cwp++;
			}
			res /= wsum;
			*drow++ = (res >= 0xffff) ? 0xffff : res;
		}
	}
					/* clean up */
	free(rowv-hrad);
	free(colwta);
	return true;
}

/* Resample a 48-bit RGB image using Gaussian average */
static int
gaussianRGB48(ImgStruct *dim, const ImgStruct *sim, float rad)
{
	const float	xdens = (float)dim->xres/sim->xres;
	const float	ydens = (float)dim->yres/sim->yres;
	const int	hrad = (int)(MAXRADM*rad/xdens + .5f);
	const int	vrad = (int)(MAXRADM*rad/ydens + .5f);
	const float	hwsca = 0x7fff / (float)hrad;
	const float	vwsca = 0x7fff / (float)vrad;
	uint32		(*rowv)[3];
	uint16		*colwta, *cwp;
	int		xd, yd;
	int		j;
					/* allocate row vector holder */
	rowv = (uint32 (*)[3])malloc(sizeof(uint32)*3*(sim->xres + 2*hrad));
	colwta = (uint16 *)malloc(sizeof(uint16)*(2*hrad+1)*dim->xres);
	if ((rowv == NULL) | (colwta == NULL)) {
		DMESG(DMCmemory, "malloc() failed in gaussianRGB48");
		return false;
	}
	rowv += hrad;
	cwp = colwta;			/* assign column weight array */
	for (xd = 0; xd < dim->xres; xd++) {
		const float	xsc = (xd+.5f)/xdens;
		const int	xs0 = (int)xsc;
		for (j = -hrad; j <= hrad; j++) {
			const float	rx = (xs0+j+.5f - xsc)*xdens;
			*cwp++ = (int)(hwsca*expf(-rx*rx/(2.f*rad*rad)) + .5f);
		}
	}
					/* filter our image */
	sprintf(dmessage_buf, "Gaussian resample of %dx%d to %dx%d 16-bit by 3-channel image",
			sim->xres, sim->yres, dim->xres, dim->yres);
	DMESG(DMCtrace, dmessage_buf);
	for (yd = 0; yd < dim->yres; yd++) {
		const float	ysc = (yd+.5f)/ydens;
		const int	ys0 = (int)ysc;
		uint16		*drow = (uint16 *)ProwPtr(dim,yd);
		uint32		wsum;
		int		xs;
					/* compute row vector at this yd */
		memset(rowv, 0, sizeof(uint32)*3*sim->xres);
		wsum = 0;
		for (j = -vrad; j <= vrad; j++) {
			const float	ry = (ys0+j+.5f - ysc)*ydens;
			const uint32	wt = (int)(vwsca*expf(-ry*ry/(2.f*rad*rad)) + .5f);
			const uint16	*srow;
			if (wt < MINIWT)
				continue;
			srow = (const uint16 *)ProwPtr(sim, (ys0+j < 0) ? 0 :
					(ys0+j >= sim->yres) ? sim->yres-1 : ys0+j);
			for (xs = 0; xs < sim->xres; xs++) {
				rowv[xs][0] += *srow++ * wt;
				rowv[xs][1] += *srow++ * wt;
				rowv[xs][2] += *srow++ * wt;
			}
			wsum += wt;
		}
		for (xs = sim->xres; xs--; ) {
		    rowv[xs][0] /= wsum;
		    rowv[xs][1] /= wsum;
		    rowv[xs][2] /= wsum;
		}
		for (xs = -hrad; xs < 0; xs++) {
			rowv[xs][0] = rowv[0][0];
			rowv[xs][1] = rowv[0][1];
			rowv[xs][2] = rowv[0][2];
		}
		for (xs = sim->xres; xs < sim->xres+hrad; xs++) {
			rowv[xs][0] = rowv[sim->xres-1][0];
			rowv[xs][1] = rowv[sim->xres-1][1];
			rowv[xs][2] = rowv[sim->xres-1][2];
		}
					/* compute pixels along scanline */
		cwp = colwta;
		for (xd = 0; xd < dim->xres; xd++) {
			const int	xs0 = (int)((xd+.5f)/xdens);
			uint32	res[3];
			wsum = 0;
			res[0] = res[1] = res[2] = 0;
			for (j = -hrad; j <= hrad; j++) {
				res[0] += rowv[xs0+j][0] * *cwp;
				res[1] += rowv[xs0+j][1] * *cwp;
				res[2] += rowv[xs0+j][2] * *cwp;
				wsum += *cwp++;
			}
			res[0] /= wsum;
			res[1] /= wsum;
			res[2] /= wsum;
			*drow++ = (res[0] >= 0xffff) ? 0xffff : res[0];
			*drow++ = (res[1] >= 0xffff) ? 0xffff : res[1];
			*drow++ = (res[2] >= 0xffff) ? 0xffff : res[2];
		}
	}
					/* clean up */
	free(rowv-hrad);
	free(colwta);
	return true;
}

/* Call appropriate image resampler */
static int
resampleImage(ImgStruct *ib, const ImgStruct *ia, int basis, float rad)
{
	int			ok = false;
	const ImgColorSpace	*target_csp = NULL;
	ImgColorSpace		input_cs;
	ImgStruct		intermed;
	/*
	 * If we are downsampling from a floating-point space
	 * into an integer space, the results are better if we
	 * do the color conversion last.  If the output image
	 * has already been allocated, then we call ourselves
	 * to resample the image before remapping the CS.
	 */
	if (ia->xres > ib->xres && (ia->csp->dtype == IDTfloat) >
					(ib->csp->dtype == IDTfloat)) {
		if (ib->img != NULL) {
			intermed = *ib;		/* resample first */
			intermed.img = NULL;
			intermed.csp = ia->csp;
			ok = resampleImage(&intermed, ia, basis, rad) &&
				PmapImage(ib, &intermed, 1.f);
			PfreeImage(&intermed);
			return ok;
		}
		target_csp = ib->csp;		/* no recursive call needed */
		ib->csp = ia->csp;
	}
	if (!PnewImage(ib, .0))			/* make sure we're allocated */
		return false;
	if (target_csp == NULL)
		target_csp = ib->csp;
	/*
	 * Perform color conversion first if we cannot
	 * convert in situ without reallocating at the end.
	 * This also ensures we resample in a floating-point
	 * space when it can improve the results.
	 */
	if (ImgPixelSize(ib->csp) != ImgPixelSize(ia->csp)) {
		intermed = *ia;
		intermed.csp = ib->csp;
		intermed.img = NULL;
		if (!PmapImage(&intermed, ia, 1.f)) {
			PfreeImage(ib);
			return false;
		}
		ia = &intermed;
	}
	PcopyCS(&input_cs, ia->csp);		/* remember starting CS */
						/* fix upsampling radius */
	if (basis == PSgaussian && (ib->xres > ia->xres) & (ib->yres > ia->yres))
		rad *= sqrt((double)(ib->xres*ib->yres)/(ia->xres*ia->yres));

	if (basis == PSnearest)			/* nearest neighbor? */
		ok = nearestN(ib, ia);
	else if (basis == PSbox)		/* box filter? */
		ok = boxDownsample(ib, ia);
	else					/* else switch on pixel type */
		switch (ImgPixelLen[ia->csp->format] << 3 | ia->csp->dtype) {
		case 030 | IDTubyte:			/* 24-bit RGB */
			ok = (basis == PSgaussian)
				? gaussianRGB24(ib, ia, rad)
				: interpRGB24(ib, ia, basis) ;
			break;
		case 010 | IDTubyte:			/* 8-bit grayscale */
			ok = (basis == PSgaussian)
				? gaussianY8(ib, ia, rad)
				: interpY8(ib, ia, basis) ;
			break;
		case 030 | IDTfloat:			/* float RGB */
			ok = (basis == PSgaussian)
				? gaussianRGB(ib, ia, rad)
				: interpRGB(ib, ia, basis) ;
			break;
		case 010 | IDTfloat:			/* float grayscale */
			ok = (basis == PSgaussian)
				? gaussianY(ib, ia, rad)
				: interpY(ib, ia, basis) ;
			break;
		case 030 | IDTushort:			/* 48-bit RGB */
			ok = (basis == PSgaussian)
				? gaussianRGB48(ib, ia, rad)
				: interpRGB48(ib, ia, basis) ;
			break;
		case 010 | IDTushort:			/* 16-bit grayscale */
			ok = (basis == PSgaussian)
				? gaussianY16(ib, ia, rad)
				: interpY16(ib, ia, basis) ;
			break;
		default:
			DMESG(DMCparameter, "Unsupported pixel type");
			break;
		}
	if (ia == &intermed)			/* done with intermediate? */
		PfreeImage(&intermed);
	if (!ok) {				/* check status */
		PfreeImage(ib);
		return false;
	}
	if (PmatchColorSpace(target_csp, &input_cs, PICMall))
		return true;			/* color spaces match */
	ib->csp = &input_cs;			/* else need to convert */
	return PconvertColorSpace(ib, target_csp, 1.f);
}

/* Peform resampling on source image and copy to destination */
int
PsizeImage(ImgStruct *ib, const ImgStruct *ia, int basis)
{
	int	interp;
						/* check format */
	if (ia == NULL || (ia->img == NULL) | (ia->csp == NULL))
		return false;
	if (ib == ia)
		return true;
	if (ib == NULL || ib->csp == NULL || (ib->xres <= 0) | (ib->yres <= 0))
		return false;
						/* check for 1-1 mapping */
	if ((ib->xres == ia->xres) & (ib->yres == ia->yres))
		return PmapImage(ib, ia, 1.f);
						/* check for self-filter */
	if (PimagesOverlap(ia, ib)) {
		DMESG(DMCparameter, "Cannot resample image onto itself");
		return false;
	}
	if (basis == PSgaussian)
		interp = ((ib->xres >= ia->xres) & (ib->yres >= ia->yres));
	else
		interp = (((int)(ib->xres*MAXINTERP) >= ia->xres) &
				((int)(ib->yres*MAXINTERP) >= ia->yres));
	if (interp)
		switch (basis) {		/* reassign basis */
		case PSbest:
			basis = PScubic;
			break;
		case PSfast:
		case PSbox:
			basis = PSnearest;
			break;
		case PSgaussian:
			basis = PSlinear;	/* same speed as bicubic */
			break;
		case PSnearest:
		case PSlinear:
		case PScubic:
			break;
		default:
			DMESG(DMCparameter, "Unknown resampling basis");
			return false;
		}
	else
		switch (basis) {		/* reassign basis */
		case PSbest:
		case PSlinear:
		case PScubic:
			basis = PSgaussian;
			break;
		case PSfast:
			basis = PSnearest;
			break;
		case PSgaussian:
		case PSnearest:
		case PSbox:
			break;
		default:
			DMESG(DMCparameter, "Unknown resampling basis");
			return false;
		}
						/* call switch routine */
	return resampleImage(ib, ia, basis, OPTRAD);
}

/* Blur image using Gaussian kernel of given radius */
int
PblurImage(ImgStruct *ib, const ImgStruct *ia, float rad)
{
	int		basis = PSgaussian;
	int		ok;
	ImgStruct	tmpImg;
						/* check format */
	if (ia == NULL || (ia->img == NULL) | (ia->csp == NULL))
		return false;
	if (ib == NULL || ib->csp == NULL)
		return false;
	if (rad <= 0)
		return false;
	tmpImg.img = NULL;
	if (ia->csp->dtype == ib->csp->dtype)
		tmpImg.csp = (ImgPixelLen[ia->csp->format] <
				ImgPixelLen[ib->csp->format]) ?
					ia->csp : ib->csp;
	else
		tmpImg.csp = (ia->csp->dtype > ib->csp->dtype) ?
					ia->csp : ib->csp;
						/* optimize for large kernel */
	if (rad >= 6.f && (ia->xres > 3*rad) & (ia->yres > 3*rad)) {
		float		sf = sqrt(rad);

		tmpImg.xres = ia->xres/sf + 1;
		tmpImg.yres = ia->yres/sf + 1;
						/* downsample & blur */
		if (!resampleImage(&tmpImg, ia, basis, sf))
			return false;
		basis = PSlinear;		/* interpolate next */
		rad = 1.f;			/* ignored, anyway */
		ia = &tmpImg;
	} else if (PimagesOverlap(ia, ib)) {	/* self-filter? */
		tmpImg.xres = ia->xres;
		tmpImg.yres = ia->yres;
						/* copy/convert */
		if (!PmapImage(&tmpImg, ia, 1.f))
			return false;
		ia = &tmpImg;
	}
	if (rad <= OPTRAD)
		DMESGF(DMCwarning, "Non-blurring radius (%f)", rad);
						/* call switch routine */
	ok = resampleImage(ib, ia, basis, rad);
	PfreeImage(&tmpImg);
	return ok;
}
