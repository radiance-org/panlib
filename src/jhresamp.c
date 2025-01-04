/*
 *  jhresamp.c
 *
 *  Created by Greg Ward on 1/2/05.
 *  Copyright 2004 Sunnybrook Technologies, <www.sunnybrooktech.com>. 
 *  All rights reserved.
 *
 *  Routines for resampling 8-bit planar image data.
 */

#define JPEG_INTERNALS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "jpeghdr.h"

#if defined(_WIN32) || defined(_WIN64)
#undef powf
#define powf(x,e)	(float)pow((double)(x),(double)(e))
#undef logf
#define logf(x)		(float)log((double)(x))
#undef expf
#define expf(x)		(float)exp((double)(x))
#endif

#define ifloor(x)	(((x) < 0) ? (int)((x)-.999999f) : (int)(x))

/* Compute Catmull-Rom interpolant vector (0-1 mapped to 0-1024) */
LOCAL(void)
comp_interp (short hv[4], const float t)
{
	const float	t2 = t*t;
	const float	t3 = t2*t;

	hv[0] = (int)(512.f*(-t3 + 2.f*t2 - t) - .5f);
	hv[1] = (int)(512.f*(3.f*t3 - 5.f*t2 + 2.f) + .5f);
	hv[2] = (int)(512.f*(-3.f*t3 + 4.f*t2 + t) + .5f);
	hv[3] = (int)(512.f*(t3 - t2) - .5f);
}

/* Compute Catmull-Rom interpolant vector (0-1 mapped to 0-255) */
LOCAL(void)
comp_interp_orig (short hv[4], const float t)
{
	const float	t2 = t*t;
	const float	t3 = t2*t;

	hv[0] = (int)(128.f*(-t3 + 2.f*t2 - t) - .5f);
	hv[1] = (int)(128.f*(3.f*t3 - 5.f*t2 + 2.f) + .5f);
	hv[2] = (int)(128.f*(-3.f*t3 + 4.f*t2 + t) + .5f);
	hv[3] = (int)(128.f*(t3 - t2) - .5f);
}

/* Resample an 8-bit planar image using bicubic interpolation */
LOCAL(void)
bicubic8 (const UINT8 *simg, const int swidth, const int sheight,
		UINT8 *dimg, const int dwidth, const int dheight)
{
	int	*rowv;
	short	*colwta, *cwp;
	int	xd, yd;
					/* allocate row vector holder */
	rowv = (int *)malloc(sizeof(int)*(swidth+4));
	colwta = (short *)malloc(sizeof(short)*4*dwidth);
	if ((rowv == NULL) | (colwta == NULL))
		return;
	rowv += 2;
	cwp = colwta;			/* assign column weight array */
	for (xd = 0; xd < dwidth; xd++) {
		float	xsc = (xd+.5f)*swidth/dwidth - .5f;
		comp_interp(cwp, xsc-ifloor(xsc));
		cwp += 4;
	}
					/* interpolate image */
	srand(1);
	for (yd = 0; yd < dheight; yd++) {
		const float	ysc = (yd+.5f)*sheight/dheight - .5f;
		const int	ys1 = ifloor(ysc);
		const UINT8	*srow1 = (ys1 <= 0) ?
						simg : simg + ys1*swidth;
		const UINT8	*srow0 = (ys1 <= 0) ?
						srow1 : srow1-swidth;
		const UINT8	*srow2 = (ys1 >= sheight-1) ?
						srow1 : srow1+swidth;
		const UINT8	*srow3 = (ys1 >= sheight-2) ?
						srow2 : srow2+swidth;
		UINT8		*drow = dimg + yd*dwidth;
		short		hvec[4];
		int		xs;
					/* premultiply scanline at this yd */
		comp_interp(hvec, ysc-ys1);

		for (xs = 0; xs < swidth; xs++) {
			rowv[xs] =	(int)*srow0++ * hvec[0] +
					(int)*srow1++ * hvec[1] +
					(int)*srow2++ * hvec[2] +
					(int)*srow3++ * hvec[3] ;
		}
		rowv[-2] = rowv[-1] = rowv[0];
		rowv[swidth+1] = rowv[swidth] = rowv[swidth-1];
					/* compute pixels along scanline */
		cwp = colwta;
		for (xd = 0; xd < dwidth; xd++) {
			float		xsc = (xd+.5f)*swidth/dwidth - .5f;
			const int	xs1 = ifloor(xsc);
			int		res;
			
			res =	rowv[xs1-1]*cwp[0] + rowv[xs1]*cwp[1] +
				rowv[xs1+1]*cwp[2] + rowv[xs1+2]*cwp[3] ;

			*drow++ = (res <= 0) ? 0 :
					(res >= 0xff00000) ? 0xff :
#if (RAND_MAX>>20)
					(res + rand()/(RAND_MAX>>20))>>20;
#else
					(res + (1L<<19))>>20;
#endif
			cwp += 4;
		}
	}
					/* clean up */
	free((void *)(rowv-2));
	free((void *)colwta);
}

/* Original (flawed) 8-bit bicubic resampling */
LOCAL(void)
bicubic8orig (const UINT8 *simg, const int swidth, const int sheight,
		UINT8 *dimg, const int dwidth, const int dheight)
{
	int	*rowv;
	short	*colwta, *cwp;
	int	xd, yd;
					/* allocate row vector holder */
	rowv = (int *)malloc(sizeof(int)*(swidth+3));
	colwta = (short *)malloc(sizeof(short)*4*dwidth);
	if ((rowv == NULL) | (colwta == NULL))
		return;
	++rowv;
	cwp = colwta;			/* assign column weight array */
	for (xd = 0; xd < dwidth; xd++) {
		float	xsc = (float)xd*swidth/dwidth;
		comp_interp_orig(cwp, xsc-(int)xsc);
		cwp += 4;
	}
					/* interpolate image */
	for (yd = 0; yd < dheight; yd++) {
		const float	ysc = (float)yd*sheight/dheight;
		const int	ys1 = (int)ysc;
		const UINT8	*srow1 = simg + ys1*swidth;
		const UINT8	*srow0 = (ys1 <= 0) ?
						srow1 : srow1-swidth;
		const UINT8	*srow2 = (ys1 >= sheight-1) ?
						srow1 : srow1+swidth;
		const UINT8	*srow3 = (ys1 >= sheight-2) ?
						srow2 : srow2+swidth;
		UINT8		*drow = dimg + yd*dwidth;
		short		hvec[4];
		int		xs;
					/* premultiply scanline at this yd */
		comp_interp_orig(hvec, ysc-ys1);

		for (xs = 0; xs < swidth; xs++) {
			rowv[xs] =	(int)*srow0++ * hvec[0] +
					(int)*srow1++ * hvec[1] +
					(int)*srow2++ * hvec[2] +
					(int)*srow3++ * hvec[3] ;
		}
		rowv[-1] = rowv[0];
		rowv[swidth+1] = rowv[swidth] = rowv[swidth-1];
					/* compute pixels along scanline */
		cwp = colwta;
		for (xd = 0; xd < dwidth; xd++) {
			const int	xs1 = xd*swidth/dwidth;
			int		res;
			
			res =	rowv[xs1-1]*cwp[0] + rowv[xs1]*cwp[1] +
				rowv[xs1+1]*cwp[2] + rowv[xs1+2]*cwp[3] ;

			*drow++ = (res <= 0) ? 0 :
					(res >= 1<<24) ? 255 : res>>16;
			cwp += 4;
		}
	}
					/* clean up */
	free((void *)(rowv-1));
	free((void *)colwta);
}

/* Resample an 8-bit planar image using Gaussian average */
LOCAL(void)
gaussian8 (const UINT8 *simg, const int swidth, const int sheight,
		UINT8 *dimg, const int dwidth, const int dheight)
{
	const float	gaussrf = -1.f/(.6f*.6f);
	const float	maxradr = 1.5f;
	const float	xdens = (float)dwidth/swidth;
	const float	ydens = (float)dheight/sheight;
	const int	hrad = (int)(maxradr/xdens + .5f);
	const int	vrad = (int)(maxradr/ydens + .5f);
	const float	hwsca = 0x7fff / (float)hrad;
	const float	vwsca = 0x7fff / (float)vrad;
	short		*colwta, *cwp;
	INT32		*rowv;
	int		xd, yd;
	int		j;
					/* allocate row vector holder */
	rowv = (INT32 *)malloc(sizeof(INT32)*(swidth + 2*hrad));
	colwta = (short *)malloc(sizeof(short)*(2*hrad+1)*dwidth);
	if ((rowv == NULL) | (colwta == NULL))
		return;
	rowv += hrad;
					/* assign column weight array */
	cwp = colwta;
	for (xd = 0; xd < dwidth; xd++) {
		const float	xsc = (xd+.5f)/xdens;
		const int	xs0 = (int)xsc;
		for (j = -hrad; j <= hrad; j++) {
			const float	rx = (xs0+j+.5f - xsc)*xdens;
			*cwp++ = (int)(hwsca*expf(gaussrf*rx*rx) + .5f);
		}
	}
					/* filter our image */
	for (yd = 0; yd < dheight; yd++) {
		const float	ysc = (yd+.5f)/ydens;
		const int	ys0 = (int)ysc;
		UINT8		*drow = dimg + yd*dwidth;
		int		wsum;
		int		xs;
					/* compute row vector at this yd */
		memset((void *)rowv, 0, sizeof(INT32)*swidth);
		wsum = 0;
		for (j = -vrad; j <= vrad; j++) {
			const float	ry = (ys0+j+.5f - ysc)*ydens;
			const INT32	wt = (int)(vwsca*expf(gaussrf*ry*ry) + .5f);
			const UINT8	*srow;
			if (wt <= 0)
				continue;
			srow = simg + swidth*((ys0+j < 0) ? 0 :
					(ys0+j >= sheight) ? sheight-1 : ys0+j);
			for (xs = 0; xs < swidth; xs++)
				rowv[xs] += *srow++ * wt;
			wsum += wt;
		}
		for (xs = swidth; xs--; )
			rowv[xs] = (rowv[xs] << 8) / wsum;
		for (xs = -hrad; xs < 0; xs++)
			rowv[xs] = rowv[0];
		for (xs = swidth; xs < swidth+hrad; xs++)
			rowv[xs] = rowv[swidth-1];
					/* compute pixels along scanline */
		cwp = colwta;
		for (xd = 0; xd < dwidth; xd++) {
			const int	xs0 = (int)((xd+.5f)/xdens);
			INT32		res = 0;
			wsum = 0;
			for (j = -hrad; j <= hrad; j++) {
				res += rowv[xs0+j] * *cwp;
				wsum += *cwp++;
			}
			res /= wsum;
			*drow++ = (res >= 0xff00) ? 0xff : res>>8;
		}
	}
					/* clean up */
	free((void *)(rowv-hrad));
	free((void *)colwta);
}

/* Resample an 8-bit planar image */
GLOBAL(void)
jpeghdr_resample8 (const UINT8 *simg, int swidth, int sheight,
		UINT8 *dimg, int dwidth, int dheight)
{
	if ((dwidth >= swidth) & (dheight >= sheight))
		bicubic8(simg, swidth, sheight, dimg, dwidth, dheight);
	else
		gaussian8(simg, swidth, sheight, dimg, dwidth, dheight);
}

/* Original (flawed) upsampling algorithm */
GLOBAL(void)
jpeghdr_resample8orig (const UINT8 *simg, int swidth, int sheight,
		UINT8 *dimg, int dwidth, int dheight)
{
	if ((dwidth >= swidth) & (dheight >= sheight))
		bicubic8orig(simg, swidth, sheight, dimg, dwidth, dheight);
	else
		gaussian8(simg, swidth, sheight, dimg, dwidth, dheight);
}
