/*
 *  phdalign.cpp
 *  panlib
 *
 *  Image alignment as required for HDR creation.
 *
 *  Created by gward on Wed Sep 19 2001.
 *  Copyright (c) 2006 Anyhere Software. All rights reserved.
 *
 */

#include <math.h>
#include "pancine.h"
#include "phdrimg.h"

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

#define P_CHAN		1		// channel to compare for alignment
#define P_MAXSHIFT	6		// move at most 1<<P_MAXSHIFT pixels
#define P_MAXDIVERGENCE	0.005		// maximum split shift divergence

// Compute fast rotation of image
static bool
fastRotate(ImgStruct *dim, const ImgStruct *sim, double degCW)
{
	if (dim->img != NULL)
		PfreeImage(dim);
	*dim = *sim;			// do not convert to Y
	dim->img = NULL;
					// XXX should do this myself
	return ProtateImage(dim, sim, degCW, PSfast, Pblack);
}

// Shrink image a factor of two in each dimension -- grayscale box filter
static bool
shrinkImage(ImgStruct *dim, const ImgStruct *sim)
{
	bool	srcIsRGB;
	if (dim->img != NULL)
		PfreeImage(dim);
	if (PmatchColorSpace(sim->csp, &ICS_sRGB, PICMptype))
		srcIsRGB = true;
	else if (PmatchColorSpace(sim->csp, &ICS_Y8, PICMptype))
		srcIsRGB = false;
	else
		return false;
	dim->xres = sim->xres / 2;
	dim->yres = sim->yres / 2;
	dim->csp = &ICS_Y8;
	if (!PnewImage(dim, .0))
		return false;
	const int	pstep = (srcIsRGB ? 3 : 1);
	uby8 *		dp = dim->img;
	int		x, y;
	for (y = 0; y < dim->yres; y++) {
		const uby8 *	sp = ProwPtr(sim, y<<1);
		if (srcIsRGB) sp += P_CHAN;
		for (x = 0; x < dim->xres; x++, sp += pstep<<1)
			*dp++ = (sp[0] + sp[pstep] +
					sp[sim->rowsize] + sp[sim->rowsize+pstep]) >> 2;
	}
	return true;
}

// Compute binary threshold and band gap images
static bool
thresholdImage(ABitMap2 *bm, ABitMap2 *bg, const ImgStruct *im, int thresh)
{
	bool	srcIsRGB;
	if (PmatchColorSpace(im->csp, &ICS_sRGB, PICMptype))
		srcIsRGB = true;
	else if (PmatchColorSpace(im->csp, &ICS_Y8, PICMptype))
		srcIsRGB = false;
	else
		return false;
	if (!bm->NewBitMap(im->xres, im->yres, false) ||
			!bg->NewBitMap(im->xres, im->yres, true))
		return false;
	const int	gaplower = thresh - 3;
	const int	gapupper = thresh + 3;
	const int	pstep = (srcIsRGB ? 3 : 1);
	int		x, y;
	for (y = 0; y < im->yres; y++) {
		const uby8 *	pp = ProwPtr(im, y);
		if (srcIsRGB) pp += P_CHAN;
		for (x = 0; x < im->xres; x++, pp += pstep)
			if (*pp > thresh) {
				bm->Set(x, y);
				if (*pp <= gapupper)
					bg->Reset(x, y);
			} else if (*pp >= gaplower)
				bg->Reset(x, y);
	}
	return true;
}

// Struct and callback for sorting bitmap differences
struct BMdiff {
	int		err;		// # bits that disagree
	int		shft[2];	// (x,y) shift
};

static int
cmpBMdiff(const void *p1, const void *p2)
{
	return (*(const BMdiff *)p1).err - (*(const BMdiff *)p2).err;
}

// Compute combined half-pixel shift to align this image at the given resolution
static bool
shiftImage(int shft[2], const ImgStruct *thisImg, const ImgStruct *thatImg,
		short thisThresh, short thatThresh, int shiftLim)
{
	shft[0] = shft[1] = 0;		// compute next coarser shift
	if ((shiftLim > 0) & (thisImg->xres >= 16) & (thisImg->yres >= 16)) {
		ImgStruct	thisShrunk, thatShrunk;
		bool		ok;
		thisShrunk.img = thatShrunk.img = NULL;
		ok = shrinkImage(&thisShrunk, thisImg) &&
				shrinkImage(&thatShrunk, thatImg);
		if (ok)
			ok = shiftImage(shft, &thisShrunk, &thatShrunk,
					thisThresh, thatThresh, shiftLim-1);
		PfreeImage(&thisShrunk);
		PfreeImage(&thatShrunk);
		if (!ok)
			return false;
	}
					// compute threshold bitmaps
	ABitMap2	thisBitmap, thatBitmap;
	ABitMap2	thisBandgap, thatBandgap;
	if (!thresholdImage(&thisBitmap, &thisBandgap, thisImg, thisThresh))
		return false;
	if (!thresholdImage(&thatBitmap, &thatBandgap, thatImg, thatThresh))
		return false;
					// compute errors for 21 possible shifts
	BMdiff		pDiff[21];
	int		n = 0;
	int		i, j;
	for (i = -2; i <= 2; i++)
		for (j = -2; j <= 2; j++) {
			if (i*i + j*j >= 8)     // skip corners
				continue;
			ABitMap2	cBitmap = thisBitmap;
			ABitMap2	cBandgap = thisBandgap;
			pDiff[n].shft[0] = shft[0] + i;
			pDiff[n].shft[1] = shft[1] + j;
			cBitmap.Shift(pDiff[n].shft[0], pDiff[n].shft[1], -1);
			cBandgap.Shift(pDiff[n].shft[0], pDiff[n].shft[1], 0);
			cBitmap ^= thatBitmap;
			cBitmap &= thatBandgap;
			cBitmap &= cBandgap;
			/*      the band gap image takes care of this
			cBitmap.ClearRect( xs<0 ? cBitmap.Width()+xs : 0,
					ys>0 ? ys : 0,
					abs(xs), cBitmap.Height()-abs(ys) );
			cBitmap.ClearRect( 0, ys<0 ? cBitmap.Height()+ys : 0,
					cBitmap.Width(), abs(ys) );
			*/
			pDiff[n++].err = cBitmap.SumTotal();
		}
					// sort based on thrshold differences
	qsort(pDiff, n, sizeof(BMdiff), cmpBMdiff);
	shft[0] = pDiff[0].shft[0] << 1;
	shft[1] = pDiff[0].shft[1] << 1;
					// average in shifts within threshold
	const int       errThresh = thisImg->xres * thisImg->yres / 4000;
	for (i = 1; i < n; i++) {
		if (pDiff[i].err - pDiff[0].err > errThresh)
			break;
		shft[0] += pDiff[i].shft[0] << 1;
		shft[1] += pDiff[i].shft[1] << 1;
	}
	shft[0] /= i;
	shft[1] /= i;
#ifdef phda_debug
	fprintf(stderr, "Error rate at level %d: %d/%d (%.1f%%) -> shift (%d/2,%d/2)\n",
			shiftLim, pDiff[0].err,
			thisBitmap.Width()*thisBitmap.Height(),
			pDiff[0].err*100./(thisBitmap.Width()*thisBitmap.Height()),
			shft[0], shft[1]);
	{
		static int      nimg = 0;
		++nimg;
		char		fname[64];
		ABitMap2	cBitmap = thisBitmap;
		ABitMap2	cBandgap = thisBandgap;
		cBitmap.Shift(shft[0], shft[1], -1);
		cBandgap.Shift(shft[0], shft[1], 0);
		cBitmap ^= thatBitmap;
		sprintf(fname, "a%02debmapA.bmp", nimg);
		WriteBitMap2(thisBitmap, fname);
		sprintf(fname, "a%02debmapB.bmp", nimg);
		WriteBitMap2(thatBitmap, fname);
		sprintf(fname, "a%02dbandgap.bmp", nimg);
		WriteBitMap2(thisBandgap, fname);
		sprintf(fname, "a%02dbefore.bmp", nimg);
		WriteBitMap2(cBitmap, fname);
		cBitmap &= thatBandgap;
		cBitmap &= cBandgap;
		sprintf(fname, "a%02dafter.bmp", nimg);
		WriteBitMap2(cBitmap, fname);
	}
#endif
	return true;
}

// Estimate rotation using the given image region and its mirror
static bool
getRotation(float *angp, const ImgStruct *thisImg, const ImgStruct *thatImg,
		short thisThresh, short thatThresh, int shiftLim,
		const ImgRect *b0p)
{
	ImgRect		b1;
	ImgStruct	imgA, imgB;
	int		shft0[2], shft1[2];
	int		dx, dy;
	double		dxr, dyr;
	bool		ok;
	
	*angp = 0;
					// compute first shift
	if (!PlinkSubimage(&imgA, thisImg, b0p) ||
			!PlinkSubimage(&imgB, thatImg, b0p)) {
		DMESG(DMCassert, "PlinkSubimage failed!");
		return false;
	}
	ok = shiftImage(shft0, &imgA, &imgB, thisThresh, thatThresh, shiftLim);
	PfreeImage(&imgA);
	PfreeImage(&imgB);
	if (!ok)
		return false;
					// get mirror position
	b1.xleft = thisImg->xres - b0p->xright;
	b1.xright = thisImg->xres - b0p->xleft;
	b1.ytop = thisImg->yres - b0p->ybottom;
	b1.ybottom = thisImg->yres - b0p->ytop;
					// compute second shift
	if (!PlinkSubimage(&imgA, thisImg, &b1) ||
			!PlinkSubimage(&imgB, thatImg, &b1)) {
		DMESG(DMCassert, "PlinkSubimage failed!");
		return false;
	}
	ok = shiftImage(shft1, &imgA, &imgB, thisThresh, thatThresh, shiftLim);
	PfreeImage(&imgA);
	PfreeImage(&imgB);
	if (!ok)
		return false;
					// get rotation from shift pair
	dx = b1.xleft - b0p->xleft;
	dy = b1.ytop - b0p->ytop;
	dxr = dx + .5*(shft1[0] - shft0[0]);	// half-pixel shifts
	dyr = dy + .5*(shft1[1] - shft0[1]);
	if (fabs(sqrt((dxr*dxr + dyr*dyr)/(dx*dx + dy*dy)) - 1.) > P_MAXDIVERGENCE)
		return false;
	*angp = (float)(180./M_PI*(atan2(dyr,dxr) - atan2((double)dy,(double)dx)));
	return true;
}

// Check angles to determine consensus average
static bool
avgAngle(float *angp, const float angl[], int nang)
{
	*angp = 0;
	if (nang < 1)
		return false;
	int	npos = 0, nneg = 0;
	int	j;
	for (j = nang; j--; ) {
		npos += (angl[j] > P_ROTHRESH);
		nneg += (angl[j] < -P_ROTHRESH);
	}
	if (npos && nneg)
		return false;
	for (j = nang; j--; )
		*angp += angl[j];
	*angp /= (float)nang;
	return true;
}

// Compute half-pixel shift and rotation between images
static bool
shiftRotImage(int shft[2], float *angp,
		const ImgStruct *thisImg, const ImgStruct *thatImg,
		short thisThresh, short thatThresh, int shiftLim)
{
	ImgRect		barbell;
	ImgStruct	rotImg;
	float		angl[3];
	int		nang;
	bool		ok;
	
	*angp = 0;
	if (shiftLim <= 2)
		return shiftImage(shft, thisImg, thatImg,
					thisThresh, thatThresh, shiftLim);
	nang = 0;			// get 3 barbell rotations
	if (thisImg->xres > thisImg->yres) {
		barbell.xleft = 0;
		barbell.xright = thisImg->xres/4;
		barbell.ytop = 0;
		barbell.ybottom = thisImg->yres/3;
		nang += getRotation(&angl[nang], thisImg, thatImg,
				thisThresh, thatThresh, shiftLim, &barbell);
		barbell.ytop = thisImg->yres/3;
		barbell.ybottom = thisImg->yres*2/3;
		nang += getRotation(&angl[nang], thisImg, thatImg,
				thisThresh, thatThresh, shiftLim, &barbell);
		barbell.ytop = thisImg->yres*2/3;
		barbell.ybottom = thisImg->yres;
		nang += getRotation(&angl[nang], thisImg, thatImg,
				thisThresh, thatThresh, shiftLim, &barbell);
	} else {
		barbell.ytop = 0;
		barbell.ybottom = thisImg->yres/4;
		barbell.xleft = 0;
		barbell.xright = thisImg->xres/3;
		nang += getRotation(&angl[nang], thisImg, thatImg,
				thisThresh, thatThresh, shiftLim, &barbell);
		barbell.xleft = thisImg->xres/3;
		barbell.xright = thisImg->xres*2/3;
		nang += getRotation(&angl[nang], thisImg, thatImg,
				thisThresh, thatThresh, shiftLim, &barbell);
		barbell.xleft = thisImg->xres*2/3;
		barbell.xright = thisImg->xres;
		nang += getRotation(&angl[nang], thisImg, thatImg,
				thisThresh, thatThresh, shiftLim, &barbell);
	}
	rotImg.img = NULL;		// rotate image
	if (!avgAngle(angp, angl, nang))
		DMESG(DMCwarning, "Rotational alignment failed");
	if (fabs(*angp) > P_ROTHRESH)
		ok = fastRotate(&rotImg, thisImg, *angp);
	else
		ok = PlinkImage(&rotImg, thisImg);
	if (!ok)
		return false;
					// shift-align (rotated) image
	ok = shiftImage(shft, &rotImg, thatImg,
				thisThresh, thatThresh, shiftLim);
	PfreeImage(&rotImg);
	return ok;
}

// Limit image offset based on pixel value limits and threshold
static bool
limitShift(int *ms, int vmin, int vthrsh, int vmax)
{
	const float	tooClose = .015f;
	const float	increment = .005f;
	float		dist;
	int		shiftLim;
	if (vmin >= vmax)
		return false;
	dist = float( (vthrsh-vmin < vmax-vthrsh) ? vthrsh-vmin : vmax-vthrsh )
			/ float(vmax - vmin);
	if (dist <= tooClose)
		return false;
	shiftLim = int((dist - tooClose)/increment);
	if (*ms > shiftLim)
		*ms = shiftLim;
	return true;
}

// Compute position and rotation offset for alignment
// Assumes ComputeUsability() already called
bool
PHDExposure::ComputeOffset(PHDExposure *that)
{
	int	maxShift = P_MAXSHIFT;
	int	shft[2];
	float	rot;

	xoff = that->xoff; yoff = that->yoff;
	degCW = that->degCW;
	if ((imr.xres != that->imr.xres) | (imr.yres != that->imr.yres)) {
		DMESG(DMCdata, "Exposure resolution mismatch");
		return false;
	}
					// pick best threshold percentile
	int	thresh0 = median[P_CHAN];
	int	thresh1 = that->median[P_CHAN];
	int	headroom, legroom, headabove, legbelow;
	if (thresh0 > thresh1) {
		headroom = maxv[P_CHAN] - thresh0;
		headabove = maxv[P_CHAN] - pct83[P_CHAN];
		legroom = thresh1 - that->minv[P_CHAN];
		legbelow = that->pct17[P_CHAN] - that->minv[P_CHAN];
	} else {
		headroom = that->maxv[P_CHAN] - thresh1;
		headabove = that->maxv[P_CHAN] - that->pct83[P_CHAN];
		legroom = thresh0 - minv[P_CHAN];
		legbelow = pct17[P_CHAN] - minv[P_CHAN];;
	}
	if (headroom < legbelow) {
		thresh0 = pct17[P_CHAN];
		thresh1 = that->pct17[P_CHAN];
	} else if (legroom < headabove) {
		thresh0 = pct83[P_CHAN];
		thresh1 = that->pct83[P_CHAN];
	}
	if (!limitShift(&maxShift, minv[P_CHAN],
				thresh0, maxv[P_CHAN]) ||
			!limitShift(&maxShift, that->minv[P_CHAN],
				thresh1, that->maxv[P_CHAN])) {
		DMESG(DMCwarning, "Exposure too extreme for alignment");
		return true;
	}
	if (!GetReady())
		return false;
	if (!that->GetReady()) {
		Release();
		return false;
	}
					// compute half-pixel shift & rotation
	if (!shiftRotImage(shft, &rot, &imr, &that->imr, thresh0, thresh1, maxShift)) {
		DMESG(DMCdata, "Cannot compute image offset");
		Release(); that->Release();
		return false;
	}
	Release(); that->Release();
	xoff -= shft[0] >> 1;
	yoff -= shft[1] >> 1;
	degCW -= rot;
#ifdef phda_debug
	fprintf(stderr, "Shift limit set to %d\n", maxShift);
	fprintf(stderr, "Offset to previous image: (%d/2,%d/2) (%f degCW)\n",
			-shft[0], -shft[1], -rot);
	fprintf(stderr, "Offset to first image: (%d,%d) (%f degCW)\n",
			xoff, yoff, degCW);
#endif
	return true;
}

