/*
 *  piwarp.c
 *  panlib
 *
 *  Pancine image warping and rotation routines.
 *
 *  Created by Greg Ward on 11/16/06.
 *  Copyright 2006 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "system.h"
#include "dmessage.h"
#include "pimage.h"
#include "tiff.h"		/* needed for uint16, etc. */

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif
#ifndef true
#define true	1
#define false	0
#endif

/* Subimage warping structure */
typedef struct {
	double	a, b, c, d, e, f, g, h;
	int	ures, vres;
	double	base_scale;
} WarpXform;

#if 0	/* XXX use bilinear rather than projective mapping to avoid tearing */
/* Compute projective transform [Heckbert 1989] */
static void
warpCompute(WarpXform *pxf, int ur, int vr, float corn[4][2])
{
	const double	x0 = corn[0][0], y0 = corn[0][1];
	const double	x1 = corn[1][0], y1 = corn[1][1];
	const double	x2 = corn[3][0], y2 = corn[3][1];
	const double	x3 = corn[2][0], y3 = corn[2][1];
	const double	xsum = x0 - x1 + x2 - x3;
	const double	ysum = y0 - y1 + y2 - y3;
	
	pxf->ures = ur;
	pxf->vres = vr;
	pxf->c = x0;
	pxf->f = y0;
					/* check for affine */
	if (fabs(xsum) <= 1e-4 && fabs(ysum) <= 1e-4) {
		pxf->a = x1 - x0;
		pxf->b = x2 - x1;
		pxf->d = y1 - y0;
		pxf->e = y2 - y1;
		pxf->g = pxf->h = .0;
	} else {			/* projective */
		double	dx1 = x1 - x2, dy1 = y1 - y2;
		double	dx2 = x3 - x2, dy2 = y3 - y2;
		double	sca = 1./(dx1*dy2 - dy1*dx2);

		pxf->g = sca*(xsum*dy2 - ysum*dx2);
		pxf->h = sca*(dx1*ysum - dy1*xsum);
		pxf->a = x1 - x0 + pxf->g*x1;
		pxf->b = x3 - x0 + pxf->h*x3;
		pxf->d = y1 - y0 + pxf->g*y1;
		pxf->e = y3 - y0 + pxf->h*y3;
	}
					/* compute base scale factor */
	pxf->base_scale = sqrt(
			.5*( pxf->a*pxf->a + pxf->d*pxf->d +
			     pxf->b*pxf->b + pxf->e*pxf->e ) /
			     (ur*ur + vr*vr) );
}

/* Apply projective transform */
static double
warpMap(double xyp[2], int ui, int vi, const WarpXform *pxf)
{
	const double	u = (ui + .5)/pxf->ures;
	const double	v = (vi + .5)/pxf->vres;
	double		sca = pxf->g*u + pxf->h*v + 1.;

	if (sca == .0)
		return .0;
	sca = 1./sca;
	xyp[0] = sca*(pxf->a*u + pxf->b*v + pxf->c) - .5;
	xyp[1] = sca*(pxf->d*u + pxf->e*v + pxf->f) - .5;

	return sca * pxf->base_scale;
}
#else
/* Compute bilinear transform [Heckbert 1989] */
static void
warpCompute(WarpXform *pxf, int ur, int vr, float corn[4][2])
{
	pxf->a = corn[0][0] - corn[1][0] - corn[2][0] + corn[3][0];
	pxf->b = -corn[0][0] + corn[1][0];
	pxf->c = -corn[0][0] + corn[2][0];
	pxf->d = corn[0][0];
	pxf->e = corn[0][1] - corn[1][1] - corn[2][1] + corn[3][1];
	pxf->f = -corn[0][1] + corn[1][1];
	pxf->g = -corn[0][1] + corn[2][1];
	pxf->h = corn[0][1];
	pxf->ures = ur;
	pxf->vres = vr;
					/* compute base scale factor */
	pxf->base_scale = sqrt(
			.5*( pxf->b*pxf->b + pxf->c*pxf->c +
			     pxf->f*pxf->f + pxf->g*pxf->g ) /
			     (ur*ur + vr*vr) );
}

/* Apply bilinear transform */
static double
warpMap(double xyp[2], int ui, int vi, const WarpXform *pxf)
{
	const double	u = (ui + .5)/pxf->ures;
	const double	v = (vi + .5)/pxf->vres;

	xyp[0] = u*v*pxf->a + u*pxf->b + v*pxf->c + pxf->d - .5;
	xyp[1] = u*v*pxf->e + u*pxf->f + v*pxf->g + pxf->h - .5;

	return pxf->base_scale;
}
#endif

/* Linear interpolant basis */
static void
linear_interp(float hv[4], const float t)
{
	hv[0] = hv[3] = 0;
	hv[1] = 1.f - t;
	hv[2] = t;
}

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

/* Compute interpolation weight array after warp, box filter in limit */
static const uby8 *
interpMap(float wta[2][4], const uby8 *pia[4][4],
		int x, int y, const WarpXform *pxf,
		const ImgStruct *ia, int basis, const uby8 *fill)
{
	const int	psiz = ImgPixelSize(ia->csp);
	int		step = 1;
	double		fpos[2];
	double		psca;
	int		x0, y0;
	int		i, j;
						/* map coordinates */
	psca = warpMap(fpos, x, y, pxf);
	if (psca <= .0)				/* map to infinity */
		return fill;
						/* nearest neighbor? */
	if ((basis == PSnearest) | (basis == PSfast)) {
		if ((fpos[0] < -.5) | (fpos[1] < -.5))
			return fill;
		x0 = (int)(fpos[0] + .5);
		y0 = (int)(fpos[1] + .5);
		if ((x0 >= ia->xres) | (y0 >= ia->yres))
			return fill;
		return PpixP(ia,x0,y0,psiz);
	}
	x0 = fpos[0] < .0 ? (int)(fpos[0]-.999999999) : (int)fpos[0];
	y0 = fpos[1] < .0 ? (int)(fpos[1]-.999999999) : (int)fpos[1];
						/* compute weights */
	if (psca >= 4.) {
		for (j = 4; j--; )
			wta[0][j] = wta[1][j] = .25f;
		step = (int)(psca*.25 + .000001);
	} else {
		if (basis == PSlinear) {
			linear_interp(wta[0], (float)(fpos[0] - x0));
			linear_interp(wta[1], (float)(fpos[1] - y0));
		} else {
			cubic_interp(wta[0], (float)(fpos[0] - x0));
			cubic_interp(wta[1], (float)(fpos[1] - y0));
		}
		if (psca > 1.) {
			const float	sf = (float)((4. - psca)/3.);
			const float	cv = (1.f - sf)*.25f;
			for (j = 4; j--; ) {
				wta[0][j] = sf*wta[0][j] + cv;
				wta[1][j] = sf*wta[1][j] + cv;
			}
		}
	}
						/* assign data pointers */
	for (j = 0; j < 4; j++) {
		const int	ys = y0 + (j-1)*step;
		int		xs;
		if ((ys < 0) | (ys >= ia->yres)) {
			pia[j][0] = pia[j][1] = pia[j][2] = pia[j][3] = fill;
			continue;
		}
		for (i = 0, xs = x0 - step; i < 4; i++, xs += step)
			if ((xs < 0) | (xs >= ia->xres))
				pia[j][i] = fill;
			else
				pia[j][i] = PpixP(ia,xs,ys,psiz);
	}
	return NULL;
}

/* Warp 24-bit RGB image using warping map */
static int
mapRGB24(ImgStruct *ib, const WarpXform *pxf,
		const ImgStruct *ia, int basis, const void *pf)
{
	uby8	fill[3];
	int	x, y, i, j;

	if (pf != NULL)
		memcpy(fill, pf, 3*sizeof(uby8));
	else
		fill[0] = fill[1] = fill[2] = 0;

	for (y = 0; y < ib->yres; y++) {
	    uby8 *	pdest = ProwPtr(ib,y);
	    for (x = 0; x < ib->xres; x++) {
		float		wtarr[2][4];
		const uby8	*inp[4][4], *sval;
		float		res[3];

		sval = interpMap(wtarr, inp, x, y, pxf, ia, basis, fill);
		if (sval != NULL) {
			*pdest++ = *sval++;
			*pdest++ = *sval++;
			*pdest++ = *sval;
			continue;
		}
		res[0] = res[1] = res[2] = 0;
		for (j = 4; j--; )
		    for (i = 4; i--; ) {
			const float	wt = wtarr[0][i] * wtarr[1][j];
			const uby8 *	inpp = inp[j][i];

			res[0] += wt * inpp[0];
			res[1] += wt * inpp[1];
			res[2] += wt * inpp[2];
		    }
		for (i = 0; i < 3; i++)
			*pdest++ = res[i] <= 0 ? 0 :
				   res[i] >= 255.f ? 255 : (int)(res[i] + .5f);
	    }
	}
	return true;
}

/* Warp 8-bit grayscale image using warping map */
static int
mapY8(ImgStruct *ib, const WarpXform *pxf,
		const ImgStruct *ia, int basis, const void *pf)
{
	uby8	fill;
	int	x, y, i, j;

	if (pf != NULL)
		fill = *(const uby8 *)pf;
	else
		fill = 0;

	for (y = 0; y < ib->yres; y++) {
	    uby8 *	pdest = ProwPtr(ib,y);
	    for (x = 0; x < ib->xres; x++) {
		float		wtarr[2][4];
		const uby8	*inp[4][4], *sval;
		float		res;

		sval = interpMap(wtarr, inp, x, y, pxf, ia, basis, &fill);
		if (sval != NULL) {
			*pdest++ = *sval;
			continue;
		}
		res = 0;
		for (j = 4; j--; )
		    for (i = 4; i--; )
			res += wtarr[0][i] * wtarr[1][j] * *inp[j][i];

		*pdest++ = res <= 0 ? 0 :
			   res >= 255.f ? 255 : (int)(res + .5f);
	    }
	}
	return true;
}

/* Warp RGB float image using warping map */
static int
mapRGB(ImgStruct *ib, const WarpXform *pxf,
		const ImgStruct *ia, int basis, const void *pf)
{
	float	fill[3];
	int	x, y, i, j;

	if (pf != NULL)
		memcpy(fill, pf, 3*sizeof(float));
	else
		fill[0] = fill[1] = fill[2] = 0;

	for (y = 0; y < ib->yres; y++) {
	    float *	pdest = (float *)ProwPtr(ib,y);
	    for (x = 0; x < ib->xres; x++) {
		float		wtarr[2][4];
		const uby8	*inp[4][4];
		const float	*sval;

		sval = (const float *)interpMap(wtarr, inp, x, y, pxf, ia,
						basis, (const uby8 *)fill);
		if (sval != NULL) {
			*pdest++ = *sval++;
			*pdest++ = *sval++;
			*pdest++ = *sval;
			continue;
		}
		pdest[0] = pdest[1] = pdest[2] = 0;
		for (j = 4; j--; )
		    for (i = 4; i--; ) {
			const float	wt = wtarr[0][i] * wtarr[1][j];
			const float *	inpp = (const float *)inp[j][i];

			pdest[0] += wt * inpp[0];
			pdest[1] += wt * inpp[1];
			pdest[2] += wt * inpp[2];
		    }
		pdest += 3;
	    }
	}
	return true;
}

/* Warp grayscale float image using warping map */
static int
mapY(ImgStruct *ib, const WarpXform *pxf,
		const ImgStruct *ia, int basis, const void *pf)
{
	float	fill;
	int	x, y, i, j;

	if (pf != NULL)
		fill = *(const float *)pf;
	else
		fill = 0;

	for (y = 0; y < ib->yres; y++) {
	    float *	pdest = (float *)ProwPtr(ib,y);
	    for (x = 0; x < ib->xres; x++) {
		float		wtarr[2][4];
		const uby8	*inp[4][4];
		const float	*sval;

		sval = (const float *)interpMap(wtarr, inp, x, y, pxf, ia,
						basis, (const uby8 *)&fill);
		if (sval != NULL) {
			*pdest++ = *sval;
			continue;
		}
		*pdest = 0;
		for (j = 4; j--; )
		    for (i = 4; i--; )
			*pdest += wtarr[0][i] * wtarr[1][j] *
					*(const float *)inp[j][i];

		pdest++;
	    }
	}
	return true;
}

/* Warp 48-bit RGB image using warping map */
static int
mapRGB48(ImgStruct *ib, const WarpXform *pxf,
		const ImgStruct *ia, int basis, const void *pf)
{
	uint16	fill[3];
	int	x, y, i, j;

	if (pf != NULL)
		memcpy(fill, pf, 3*sizeof(uint16));
	else
		fill[0] = fill[1] = fill[2] = 0;

	for (y = 0; y < ib->yres; y++) {
	    uint16 *	pdest = (uint16 *)ProwPtr(ib,y);
	    for (x = 0; x < ib->xres; x++) {
		float		wtarr[2][4];
		const uby8	*inp[4][4];
		const uint16	*sval;
		float		res[3];

		sval = (const uint16 *)interpMap(wtarr, inp, x, y, pxf, ia,
						basis, (const uby8 *)fill);
		if (sval != NULL) {
			*pdest++ = *sval++;
			*pdest++ = *sval++;
			*pdest++ = *sval;
			continue;
		}
		res[0] = res[1] = res[2] = 0;
		for (j = 4; j--; )
		    for (i = 4; i--; ) {
			const float	wt = wtarr[0][i] * wtarr[1][j];
			const uint16 *	inpp = (const uint16 *)inp[j][i];

			res[0] += wt * inpp[0];
			res[1] += wt * inpp[1];
			res[2] += wt * inpp[2];
		    }
		for (i = 0; i < 3; i++)
			*pdest++ = res[i] <= 0 ? 0 :
				   res[i] >= 65535.f ? 65535 : (int)(res[i] + .5f);
	    }
	}
	return true;
}

/* Warp 16-bit grayscale image using warping map */
static int
mapY16(ImgStruct *ib, const WarpXform *pxf,
		const ImgStruct *ia, int basis, const void *pf)
{
	uint16	fill;
	int	x, y, i, j;

	if (pf != NULL)
		fill = *(const uint16 *)pf;
	else
		fill = 0;

	for (y = 0; y < ib->yres; y++) {
	    uint16 *	pdest = (uint16 *)ProwPtr(ib,y);
	    for (x = 0; x < ib->xres; x++) {
		float		wtarr[2][4];
		const uby8	*inp[4][4];
		const uint16	*sval;
		float		res;

		sval = (const uint16 *)interpMap(wtarr, inp, x, y, pxf, ia,
						basis, (const uby8 *)&fill);
		if (sval != NULL) {
			*pdest++ = *sval;
			continue;
		}
		res = 0;
		for (j = 4; j--; )
		    for (i = 4; i--; )
			res += wtarr[0][i] * wtarr[1][j] *
					*(const uint16 *)inp[j][i];

		*pdest++ = res <= 0 ? 0 :
			   res >= 65535.f ? 65535 : (int)(res + .5f);
	    }
	}
	return true;
}

/* Resample using warp over image quad */
static int
doWarp(ImgStruct *ib, const ImgStruct *ia,
			float corn[4][2], int basis, const void *pf)
{
	int		ok = false;
	WarpXform	myXF;
#if 0
						/* check arguments */
	if (ia == NULL || (ia->img == NULL) | (ia->csp == NULL))
		return false;
	if (ib == NULL || ib->csp == NULL)
		return false;
	if (!PmatchColorSpace(ib->csp, ia->csp, PICMptype)) {
		DMESG(DMCparameter, "Incompatible pixel type");
		return false;
	}
	if (corn == NULL)
		return false;
#endif
	if (!PnewImage(ib, .0))			/* check allocation */
		return false;
						/* compute mapping */
	warpCompute(&myXF, ib->xres, ib->yres, corn);
						/* apply warp */
	switch (ImgPixelLen[ib->csp->format] << 3 | ib->csp->dtype) {
	case 030 | IDTubyte:			/* 24-bit RGB */
		ok = mapRGB24(ib, &myXF, ia, basis, pf);
		break;
	case 010 | IDTubyte:			/* 8-bit grayscale */
		ok = mapY8(ib, &myXF, ia, basis, pf);
		break;
	case 030 | IDTfloat:			/* float RGB */
		ok = mapRGB(ib, &myXF, ia, basis, pf);
		break;
	case 010 | IDTfloat:			/* float grayscale */
		ok = mapY(ib, &myXF, ia, basis, pf);
		break;
	case 030 | IDTushort:			/* 48-bit RGB */
		ok = mapRGB48(ib, &myXF, ia, basis, pf);
		break;
	case 010 | IDTushort:			/* 16-bit grayscale */
		ok = mapY16(ib, &myXF, ia, basis, pf);
		break;
	default:
		DMESG(DMCparameter, "Unsupported pixel type");
		break;
	}
	if (!ok)
		PfreeImage(ib);
	return ok;
}

/*
 * Warp an image based on the given grid of source positions in scanline order.
 * Fill border with given pixel pfv.
 */
int
PwarpImage(ImgStruct *ib, const ImgStruct *ia,
			float warpGrid[][2], const int whres, const int wvres,
			int basis, PixelVal pfv)
{
	const void *		pf = NULL;
	const ImgColorSpace *	target_space;
	float			downsamp;
	ImgStruct		srcImg;
	int			h, v;
						/* check arguments */
	if (ia == NULL || (ia->img == NULL) | (ia->csp == NULL))
		return false;
	if (ib == NULL || ib->csp == NULL)
		return false;
	if ((ib->xres <= 0) | (ib->yres <= 0))
		return false;
	if ((warpGrid == NULL) | (whres <= 1) | (wvres <= 1))
		return false;
	if (PimagesOverlap(ia, ib)) {
		DMESG(DMCparameter, "Cannot warp image onto itself");
		return false;
	}
						/* get downsampling rate */
	if ((basis == PSnearest) | (basis == PSfast)) {
		downsamp = 1.f;
	} else {
		downsamp = ia->xres*ia->yres;
		for (v = wvres; --v > 0; )
		    for (h = whres; --h > 0; ) {
			float	srate;
			srate = (warpGrid[v*whres + h][0] -
					warpGrid[v*whres + (h-1)][0]) *
				(warpGrid[v*whres + h][1] -
					warpGrid[(v-1)*whres + h][1]) ;
			if (srate < .0) srate = -srate;
			if (srate < downsamp)
				downsamp = srate;
		    }
		downsamp = (float)sqrt( downsamp *
				(double)((whres-1)*(wvres-1)) /
							(ib->xres*ib->yres) );
		if (downsamp < 1.f)
			downsamp = 1.f;
	}
						/* create intermediate image */
	if (downsamp > 1.5f ||
			!PmatchColorSpace(ib->csp, ia->csp, PICMptype)) {
		srcImg.xres = (int)((double)ia->xres/downsamp + .5);
		srcImg.yres = (int)((double)ia->yres/downsamp + .5);
		srcImg.csp = ib->csp;
		srcImg.img = NULL;
		if (!PsizeImage(&srcImg, ia, basis))
			return false;
	} else {
		PlinkImage(&srcImg, ia);
		downsamp = 1.f;
	}
	if (!PnewImage(ib, .0)) {		/* make sure we're allocated */
		PfreeImage(&srcImg);
		return false;
	}
	target_space = ib->csp;
						/* figure out fill pixel */
	if (pfv.csp != NULL) {
		pfv = PconvertPixel(pfv, srcImg.csp);
		pf = (const void *)pfv.v.b;
	}
	sprintf(dmessage_buf, "Warping %dx%d image to %dx%d using %dx%d sample grid",
			ia->xres, ia->yres, ib->xres, ib->yres, whres, wvres);
	DMESG(DMCtrace, dmessage_buf);
	for (v = wvres; --v > 0; )		/* run through quads */
	    for (h = whres; --h > 0; ) {
		float		corner[4][2];
		ImgRect		destRect;
		ImgStruct	destImg;
		corner[0][0] = warpGrid[(v-1)*whres + (h-1)][0] / downsamp;
		corner[0][1] = warpGrid[(v-1)*whres + (h-1)][1] / downsamp;
		corner[1][0] = warpGrid[(v-1)*whres + h][0] / downsamp;
		corner[1][1] = warpGrid[(v-1)*whres + h][1] / downsamp;
		corner[2][0] = warpGrid[v*whres + (h-1)][0] / downsamp;
		corner[2][1] = warpGrid[v*whres + (h-1)][1] / downsamp;
		corner[3][0] = warpGrid[v*whres + h][0] / downsamp;
		corner[3][1] = warpGrid[v*whres + h][1] / downsamp;
		destRect.xleft = (h-1) * ib->xres / (whres-1);
		destRect.xright = h * ib->xres / (whres-1);
		destRect.ytop = (v-1) * ib->yres / (wvres-1);
		destRect.ybottom = v * ib->yres / (wvres-1);
		destImg.img = NULL;
		if (!PlinkSubimage(&destImg, ib, &destRect) ||
				!doWarp(&destImg, &srcImg, corner, basis, pf)) {
			PfreeImage(&srcImg);
			PfreeImage(ib);
			return false;
		}
		PfreeImage(&destImg);		/* free reference */
	    }
	v = PmatchColorSpace(target_space, srcImg.csp, PICMall);
	PfreeImage(&srcImg);			/* unlink intermediate */
	if (v)
		return true;			/* color spaces match */
	ib->csp = ia->csp;			/* else need to convert */
	return PconvertColorSpace(ib, target_space, 1.f);
}

/* Rotate image clockwise about its center, filling border with pfv */
int
ProtateImage(ImgStruct *ib, const ImgStruct *ia,
				double degCW, int basis, PixelVal pfv)
{
	const double	radiansCW = M_PI/180. * degCW;
	const double	radDiag = .5*sqrt((double)(ia->xres*ia->xres +
						ia->yres*ia->yres));
	const double	ang3 = atan2((double)ia->yres, (double)ia->xres);
	double		ang;
	float		rotMap[4][2];
						/* compute rotation map */
	ang = -M_PI + ang3 - radiansCW;
	rotMap[0][0] = .5*ia->xres + radDiag*cos(ang);
	rotMap[0][1] = .5*ia->yres + radDiag*sin(ang);
	ang = -ang3 - radiansCW;
	rotMap[1][0] = .5*ia->xres + radDiag*cos(ang);
	rotMap[1][1] = .5*ia->yres + radDiag*sin(ang);
	ang = M_PI - ang3 - radiansCW;
	rotMap[2][0] = .5*ia->xres + radDiag*cos(ang);
	rotMap[2][1] = .5*ia->yres + radDiag*sin(ang);
	ang = ang3 - radiansCW;
	rotMap[3][0] = .5*ia->xres + radDiag*cos(ang);
	rotMap[3][1] = .5*ia->yres + radDiag*sin(ang);
						/* perform rotation */
	return PwarpImage(ib, ia, rotMap, 2, 2, basis, pfv);
}
