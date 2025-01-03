/*
 *  phisto.c
 *  pan
 *
 *  Image histogram routines.
 *
 *  Created by Greg Ward on 6/12/08.
 *  Copyright 2008 Anyhere Software. All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "system.h"
#include "dmessage.h"
#include "pimage.h"
#include "tiff.h"		/* needed for uint16, etc. */

#ifndef true
#define true	1
#define false	0
#endif

/* Look for extreme values in byte stream */
static void
extremeBytes(uby8 *minmax, const int nchan, const int step,
				int nsteps, const uby8 *slp)
{
	int	c;
	while (nsteps-- > 0) {
		for (c = nchan; c-- > 0; ) {
			if (minmax[c*2] > slp[c])
				minmax[c*2] = slp[c];
			if (minmax[c*2+1] < slp[c])
				minmax[c*2+1] = slp[c];
		}
		slp += step;
	}
}

/* Add byte values to histogram */
static unsigned long
addBytes(const uby8 *minmax, unsigned long *histo, const int hlen,
		const int nchan, const int step, int nsteps, const uby8 *slp)
{
	unsigned long	cnt;
	int		c;
					/* special case */
	if ((nchan == 1) & (minmax[0] == 0) && minmax[1] == 255) {
		for (cnt = 0; nsteps-- > 0; slp += step) {
			++histo[*slp];
			++cnt;
		}
		return cnt;
	}
	for (cnt = 0; nsteps-- > 0; slp += step)
		for (c = nchan; c-- > 0; ) {
			if (slp[c] < minmax[c*2])
				continue;
			if (slp[c] > minmax[c*2+1])
				continue;
			++histo[c*hlen +
				hlen*(slp[c]-minmax[c*2])/
					(1+minmax[c*2+1]-minmax[c*2])];
			++cnt;
		}
	return cnt;
}

/* Compute byte histogram */
static unsigned long
byteHisto(uby8 *minmax, unsigned long *hist, const int hlen,
			const ImgStruct *ia, PHistoChan chan)
{
	const int	nchan = (chan == PHCrgb3) ? 3 : 1;
	const int	newextr = (minmax[0] >= minmax[1]);
	int		offset = (chan <= PHCblue) ? chan : 0;
	int		step = ImgPixelLen[ia->csp->format];
	int		nsteps = ia->xres;
	unsigned long	tot = 0;
	ImgStruct	myImg;
	int		y, i;
	
	if (nchan > step) {
		DMESG(DMCparameter, "Bad color space for PHCrgb3");
		return 0;
	}
	if ((chan == PHCrgb) & (step > 1)) {
		nsteps *= step;
		step = 1;
	} else if (offset >= step)
		offset = step-1;
					/* need luminance image? */
	if ((chan == PHCluminance) & (ia->csp->format != IPFy)) {
		ImgColorSpace	newCS;
		PcopyCS(&newCS, ia->csp);
		newCS.format = IPFy;
		myImg.csp = &newCS;
		myImg.img = NULL;
		if (!PmapImage(&myImg, ia, 1.f))
			return 0;
		ia = &myImg;
		step = 1;
	}
	if (newextr | (hlen <= 0)) {	/* need to compute min & max? */
		for (i = newextr*nchan; i--; ) {
			minmax[i*2] = 0xff;
			minmax[i*2+1] = 0;
		}
		for (y = ia->yres; y--; )
			extremeBytes(minmax, nchan, step, nsteps,
					(const uby8 *)ProwPtr(ia,y) + offset);
		if (hlen <= 0) {
			tot = 1;
			goto done;
		}
					/* starting fresh */
		memset(hist, 0, sizeof(unsigned long)*nchan*hlen);
	}
					/* add a scanline at a time */
	for (y = ia->yres; y--; )
		tot += addBytes(minmax, hist, hlen, nchan, step, nsteps,
					(const uby8 *)ProwPtr(ia,y) + offset);
done:
	if (ia == &myImg)
		PfreeImage(&myImg);
	return tot;
}

/* Look for extreme values in short stream */
static void
extremeShorts(uint16 *minmax, const int nchan, const int step,
				int nsteps, const uint16 *slp)
{
	int	c;
	while (nsteps-- > 0) {
		for (c = nchan; c-- > 0; ) {
			if (minmax[c*2] > slp[c])
				minmax[c*2] = slp[c];
			if (minmax[c*2+1] < slp[c])
				minmax[c*2+1] = slp[c];
		}
		slp += step;
	}
}

/* Add short values to histogram */
static unsigned long
addShorts(const uint16 *minmax, unsigned long *histo, const int hlen,
		const int nchan, const int step, int nsteps, const uint16 *slp)
{
	unsigned long	cnt;
	int		c;
					/* special case */
	if ((nchan == 1) & (minmax[0] == 0) && minmax[1] == 0xffff) {
		for (cnt = 0; nsteps-- > 0; slp += step) {
			++histo[*slp];
			++cnt;
		}
		return cnt;
	}
	for (cnt = 0; nsteps-- > 0; slp += step)
		for (c = nchan; c-- > 0; ) {
			if (slp[c] < minmax[c*2])
				continue;
			if (slp[c] > minmax[c*2+1])
				continue;
			++histo[c*hlen +
				(unsigned long)hlen*(slp[c]-minmax[c*2])/
					(1+minmax[c*2+1]-minmax[c*2])];
			++cnt;
		}
	return cnt;
}

/* Compute short histogram */
static unsigned long
shortHisto(uint16 *minmax, unsigned long *hist, const int hlen,
			const ImgStruct *ia, PHistoChan chan)
{
	const int	nchan = (chan == PHCrgb3) ? 3 : 1;
	const int	newextr = (minmax[0] >= minmax[1]);
	int		offset = (chan <= PHCblue) ? chan : 0;
	int		step = ImgPixelLen[ia->csp->format];
	int		nsteps = ia->xres;
	unsigned long	tot = 0;
	ImgStruct	myImg;
	int		y, i;
	
	if (nchan > step) {
		DMESG(DMCparameter, "Bad color space for PHCrgb3");
		return 0;
	}
	if ((chan == PHCrgb) & (step > 1)) {
		nsteps *= step;
		step = 1;
	} else if (offset >= step)
		offset = step-1;
					/* need luminance image? */
	if ((chan == PHCluminance) & (ia->csp->format != IPFy)) {
		ImgColorSpace	newCS;
		PcopyCS(&newCS, ia->csp);
		newCS.format = IPFy;
		myImg.csp = &newCS;
		myImg.img = NULL;
		if (!PmapImage(&myImg, ia, 1.f))
			return 0;
		ia = &myImg;
		step = 1;
	}
	if (newextr | (hlen <= 0)) {	/* need to compute min & max? */
		for (i = newextr*nchan; i--; ) {
			minmax[i*2] = 0xffff;
			minmax[i*2+1] = 0;
		}
		for (y = ia->yres; y--; )
			extremeShorts(minmax, nchan, step, nsteps,
					(const uint16 *)ProwPtr(ia,y) + offset);
		if (hlen <= 0) {
			tot = 1;
			goto done;
		}
					/* starting fresh */
		memset(hist, 0, sizeof(unsigned long)*nchan*hlen);
	}
					/* add a scanline at a time */
	for (y = ia->yres; y--; )
		tot += addShorts(minmax, hist, hlen, nchan, step, nsteps,
					(const uint16 *)ProwPtr(ia,y) + offset);
done:
	if (ia == &myImg)
		PfreeImage(&myImg);
	return tot;
}

/* Look for extreme values in float stream */
static void
extremeFloats(float *minmax, const int nchan, const int step,
			int nsteps, const float *slp, const int posonly)
{
	int	c;

	while (nsteps-- > 0) {
		for (c = nchan; c-- > 0; ) {
			if (posonly && slp[c] <= 0)
				continue;
			if (minmax[c*2] > slp[c])
				minmax[c*2] = slp[c];
			if (minmax[c*2+1] < slp[c])
				minmax[c*2+1] = slp[c];
		}
		slp += step;
	}
}

/* Add float values to histogram */
static unsigned long
addFloats(const float *minmax, unsigned long *histo, const int hlen,
		const int nchan, const int step, int nsteps, const float *slp)
{
	unsigned long	cnt;
	int		c;

	for (cnt = 0; nsteps-- > 0; slp += step)
		for (c = nchan; c-- > 0; ) {
			if (slp[c] < minmax[c*2])
				continue;
			if (slp[c] > minmax[c*2+1])
				continue;
			++histo[c*hlen +
				(int)(hlen*(slp[c]-minmax[c*2])*0.999999f/
					(minmax[c*2+1]-minmax[c*2]))];
			++cnt;
		}
	return cnt;
}

/* Compute float histogram */
static unsigned long
floatHisto(float *minmax, unsigned long *hist, const int hlen,
			const ImgStruct *ia, PHistoChan chan)
{
	const int	nchan = (chan == PHCrgb3) ? 3 : 1;
	const int	newextr = (minmax[0] >= minmax[1]);
	int		offset = (chan <= PHCblue) ? chan : 0;
	int		step = ImgPixelLen[ia->csp->format];
	int		nsteps = ia->xres;
	float		myMinMax[2];
	unsigned long	tot = 0;
	ImgStruct	myImg;
	ImgColorSpace	newCS;
	int		y, i;
	
	if (nchan > step) {
		DMESG(DMCparameter, "Bad color space for PHCrgb3");
		return 0;
	}
	if ((chan == PHCrgb) & (step > 1)) {
		nsteps *= step;
		step = 1;
	} else if (chan == PHCluminance && ia->csp->format == IPFxyz &&
				(hlen <= 1) | (ia->csp->logorig > 0))
		offset = 1;		/* special case */
	else if (offset >= step)
		offset = step-1;
	myImg.img = NULL;		/* compute/update min & max */
	if (newextr | (hlen <= 0)) {
		for (i = newextr*nchan; i--; ) {
			minmax[i*2] = 1e10f;
			minmax[i*2+1] = -1e10f;
		}
					/* need to convert to Y channel? */
		if (!offset & (chan == PHCluminance) && ia->csp->format != IPFy) {
			PcopyCS(&newCS, ia->csp);
			newCS.format = IPFy;
			myImg.csp = &newCS;
			if (!PmapImage(&myImg, ia, 1.f))
				return 0;
			ia = &myImg;
			step = 1;
		}
		for (y = ia->yres; y--; )
			extremeFloats(minmax, nchan, step, nsteps,
					(const float *)ProwPtr(ia,y) + offset,
					chan==PHCluminance & ia->csp->logorig==0);
					/* make sure we got something */
		for (i = newextr*nchan; i--; )
			if (minmax[i*2] > minmax[i*2+1])
				return 0;
		if (hlen <= 0)		/* no histogram? */
			return 1;
					/* else starting fresh */
		memset(hist, 0, sizeof(unsigned long)*nchan*hlen);
	}
					/* need (log) luminance image? */
	if (!offset & (chan == PHCluminance) && (ia->csp->format != IPFy ||
				(hlen > 1) & (ia->csp->logorig == 0))) {
		PcopyCS(&newCS, ia->csp);
		newCS.format = IPFy;	/* need to convert to log space? */
		if ((hlen > 1) & (ia->csp->logorig == 0)) {
			if (minmax[0] <= 0)
				return 0;
			if (!PmapPixels((uby8 *)myMinMax, (uby8 *)minmax, 2,
					PgetColorConv(&ICS_YdB, &newCS)))
				return 0;
			minmax = myMinMax;
			PcopyCS(&newCS, &ICS_YdB);
		}
		if (ia != &myImg) {	/* first conversion? */
			myImg.csp = &newCS;
			if (!PmapImage(&myImg, ia, 1.f))
				return 0;
			ia = &myImg;
			step = 1;
		} else if (!PconvertColorSpace(&myImg, &newCS, 1.f))
			return 0;
	}
					/* add a scanline at a time */
	for (y = ia->yres; y--; )
		tot += addFloats(minmax, hist, hlen, nchan, step, nsteps,
					(const float *)ProwPtr(ia,y) + offset);
	PfreeImage(&myImg);
	return tot;
}

/* Compute histogram for image:
 *  If minmax[0] < minmax[1] on call, then these are taken as the
 *  historgram range and the histogram is only appended (not cleared).
 *  Otherwise, the minmax values are set based on the image and the
 *  histogram is cleared to all zeroes before the tally begins.
 *  If hlen <= 0, then the minmax calculation is all that happens, and
 *  and the extrema will be updated even if starting minmax[0] < minmax[1].
 *  The actual type of the minmax array is determined by ia->csp->dtype.
 *  If chan == PHCrgb3, the extrema array is dimensioned minmax[3][2],
 *  the histogram is hist[3][histlen].  If chan == PHCluminance and
 *  ia->csp->dtype == IDTfloat, then the histogram is partitioned using
 *  logarithmic steps, i.e.:
 *	hist_bin_low[i] = min*pow(max/min, (double)i/hlen)
 *  and the minimum is set to the smallest positive Y value.
 *  The total count of pixels added into the histogram during the
 *  call is returned, or zero if there was an error.  This total will be
 *  less than ia->xres*ia->yres if pixels are outside the minmax range.
 */
unsigned long
PcomputeHisto(void *minmax, unsigned long *hist, int hlen, const ImgStruct *ia,
					PHistoChan chan)
{
	if ((minmax == NULL) | (ia == NULL))
		return 0;
	if ((ia->csp == NULL) | (ia->img == NULL))
		return 0;
	if ((hlen > 0) & (hist == NULL))
		return 0;
	switch (ia->csp->dtype) {
	case IDTubyte:
		return byteHisto((uby8 *)minmax, hist, hlen, ia, chan);
	case IDTushort:
		return shortHisto((uint16 *)minmax, hist, hlen, ia, chan);
	case IDTfloat:
		return floatHisto((float *)minmax, hist, hlen, ia, chan);
	}
	DMESG(DMCparameter, "Unsupported image data type for histogram");
	return 0;
}

/* Compute pixel percentiles:
 *  The pctl[] array specifies the desired percentiles between 0 and 100.
 *  Results are returned in the res[] array, which is dimensioned res[3][n]
 *  if chan == PHCrgb3.
 */
int
PcomputePercentiles(void *res, const float *pctl, int n,
					const ImgStruct *ia, PHistoChan chan)
{
#define	HLEN		10240
	int		hlen = HLEN;
	unsigned long	histo[3*HLEN], htot = 0;
	union {
		uby8	b[3*2];
		uint16	s[3*2];
		float	f[3*2];
	}		mm;
	int		i, j, k;
						/* check arguments */
	if (n <= 0)
		return false;
	if ((res == NULL) | (pctl == NULL) | (ia == NULL) || ia->csp == NULL)
		return false;
	for (k = n; k--; )
		if ((pctl[k] < 0) | (pctl[k] > 100.f)) {
			DMESGF(DMCparameter, "Illegal percentile: %f", pctl[i]);
			return false;
		}
	switch (ia->csp->dtype) {		/* compute histogram */
	case IDTubyte:
		mm.b[0] = 0; mm.b[1] = 255;	/* special case */
		if (chan == PHCrgb3) {
			mm.b[4] = mm.b[2] = mm.b[0];
			mm.b[5] = mm.b[3] = mm.b[1];
		}
		hlen = 256;
		memset(histo, 0, 3*sizeof(histo[0])*hlen);
		htot = byteHisto(mm.b, histo, hlen, ia, chan);
		break;
	case IDTushort:
		mm.s[0] = 1; mm.s[1] = 0;
		htot = shortHisto(mm.s, histo, hlen, ia, chan);
		break;
	case IDTfloat:
		mm.f[0] = 1; mm.f[1] = 0;
		htot = floatHisto(mm.f, histo, hlen, ia, chan);
		break;
	default:
		DMESG(DMCparameter, "Unsupported image data type");
		return false;
	}
	if (!htot)
		return false;
	for (k = n; k--; )			/* derive percentiles */
		for (j = chan==PHCrgb3 ? 3 : 1; j--; ) {
			unsigned long	cnt = pctl[k]*.01*htot + .5;
			unsigned long *	hp = histo + j*hlen;
			float		pos;
			for (i = 0; cnt > *hp; i++)
				cnt -= *hp++;
			if (cnt)
				pos = (i + cnt/(float)*hp) / hlen;
			else
				pos = (float)i / hlen;
						/* convert position to value */
			switch (ia->csp->dtype) {
			case IDTubyte:
				((uby8 *)res)[j*n+k] = 255.999f*pos;
				break;
			case IDTushort:
				((uint16 *)res)[j*n+k] = mm.s[2*j]*(1.f-pos) +
							mm.s[2*j+1]*pos;
				break;
			case IDTfloat:
				if (chan == PHCluminance &&
						ia->csp->logorig == 0) {
					((float *)res)[j*n+k] = mm.f[2*j]*pow(
							mm.f[2*j+1]/mm.f[2*j], pos);
					break;
				}
				((float *)res)[j*n+k] = mm.f[2*j]*(1.f-pos) +
							mm.f[2*j+1]*pos;
				break;
			}
		}
	return true;
#undef HLEN
}
