/*
 *  pdraw.cpp
 *  panlib
 *
 *  Drawing and stamping function implementations.
 *
 *  Created by Greg Ward on 1/24/19.
 *  Copyright 2019 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "pdraw.h"
#include "dmessage.h"
#include "astack.h"

static inline int
iround(double x)
{
	return int(x + .5) - int(x < -.5);
}

// Draw (partial) circle into bitmap
int
BdrawCircle(ABitMap2 *map, int xc, int yc, int rad, int sc, bool val)
{
	if (!map | !sc) return 0;
	if (rad <= 0)
		return map->TestAndSet(xc, yc, val);
	int	np = 0;
	int	x = rad;
	int	y = 0;
	int	dx = 1;
	int	dy = 1;
	int	err = dx - ((rad+1)<<1);

	while (x >= y) {
		if (sc & QSlowerRight) {
			np += map->TestAndSet(xc+x, yc+y, val);
			np += map->TestAndSet(xc+y, yc+x, val);
		}
		if (sc & QSlowerLeft) {
			np += map->TestAndSet(xc-x, yc+y, val);
			np += map->TestAndSet(xc-y, yc+x, val);
		}
		if (sc & QSupperLeft) {
			np += map->TestAndSet(xc-x, yc-y, val);
			np += map->TestAndSet(xc-y, yc-x, val);
		}
		if (sc & QSupperRight) {
			np += map->TestAndSet(xc+x, yc-y, val);
			np += map->TestAndSet(xc+y, yc-x, val);
		}
		if (err <= 0) {
			y++;
			err += dy;
			dy += 2;
		}
		if (err > 0) {
			x--;
			dx += 2;
			err += dx - ((rad+1) << 1);
		}
	}
	return np;
}

// Draw a line into bitmap
int
BdrawLine(ABitMap2 *map, int x0, int y0, int x1, int y1, bool val)
{
	if (!map) return 0;

	int	np = 0;
	double	xstep, ystep;
	int	n;

	if (abs(y1 - y0) > abs(x1 - x0)) {
		n = abs(y1 - y0) + 1;
		ystep = 2*(y1>y0) - 1;
		xstep = (x1 - x0)/(n-1.);
	} else {
		n = abs(x1 - x0) + 1;
		xstep = 2*(x1>x0) - 1;
		ystep = (y1 - y0)/(n-0.99999999999999);
	}
	while (n--)
		np += map->TestAndSet(iround(x0 + n*xstep), iround(y0 + n*ystep), val);
	return np;
}

// Draw a fat line into bitmap
int
BfillLine(ABitMap2 *map, int x0, int y0, int x1, int y1, int rad, bool val)
{
	if (rad <= 0)			// no width to line??
		return BdrawLine(map, x0, y0, x1, y1, val);

	const double	len = sqrt((double)(x1-x0)*(x1-x0) +
					(double)(y1-y0)*(y1-y0));
	int	np = 0;
	if (len <= .75) {		// just a fat point?
		np += BdrawCircle(map, x0, y0, rad, QSfull, val);
		np += BfloodFill(map, x0, y0, val);
		return np;
	}
	const int	dx = iround(rad*(y0-y1)/len);
	const int	dy = iround(rad*(x1-x0)/len);
					// draw outer edges
	np += BdrawLine(map, x0+dx, y0+dy, x1+dx, y1+dy, val);
	np += BdrawLine(map, x0-dx, y0-dy, x1-dx, y1-dy, val);
	if (dx > 0) {			// draw appropriate 3/4 circle end-caps
		if (dy > 0) {
			np += BdrawCircle(map, x0, y0, rad, QSleft+QSlowerRight, val);
			np += BdrawCircle(map, x1, y1, rad, QSright+QSupperLeft, val);
		} else {
			np += BdrawCircle(map, x0, y0, rad, QSright+QSlowerLeft, val);
			np += BdrawCircle(map, x1, y1, rad, QSleft+QSupperRight, val);
		}
	} else {
		if (dy > 0) {
			np += BdrawCircle(map, x0, y0, rad, QSleft+QSupperRight, val);
			np += BdrawCircle(map, x1, y1, rad, QSright+QSlowerLeft, val);
		} else {
			np += BdrawCircle(map, x0, y0, rad, QSright+QSupperLeft, val);
			np += BdrawCircle(map, x1, y1, rad, QSleft+QSlowerRight, val);
		}
	}
					// flood-fill from likely seed points
	np += BfloodFill(map, (x0+x1)>>1, (y0+y1)>>1, val);
	np += BfloodFill(map, x0, y0, val);
	np += BfloodFill(map, x1, y1, val);
	return np;
}

// Draw a (rounded) rectangle into bitmap
int
BdrawRectangle(ABitMap2 *map, int x0, int y0, int x1, int y1, int rad, bool val)
{
	if (!map) return 0;

	int	i;
	if (x0 > x1) { i = x0; x0 = x1; x1 = i; }
	if (y0 > y1) { i = y0; y0 = y1; y1 = i; }
	rad *= (rad > 0);
	if (x1 - x0 < 2*rad)
		rad = (x1 - x0)/2;
	if (y1 - y0 < 2*rad)
		rad = (y1 - y0)/2;
	int	np = 0;
	np += BdrawLine(map, x0+rad, y0, x1-rad, y0, val);
	np += BdrawLine(map, x0+rad, y1, x1-rad, y1, val);
	np += BdrawLine(map, x0, y0+rad, x0, y1-rad, val);
	np += BdrawLine(map, x1, y0+rad, x1, y1-rad, val);
	if (rad) {
		np += BdrawCircle(map, x0+rad, y0+rad, rad, QSupperLeft, val);
		np += BdrawCircle(map, x1-rad, y0+rad, rad, QSupperRight, val);
		np += BdrawCircle(map, x0+rad, y1-rad, rad, QSlowerLeft, val);
		np += BdrawCircle(map, x1-rad, y1-rad, rad, QSlowerRight, val);
	}
	return np;
}

// Flood-fill area starting from seed point which must be inside closed region
int
BfloodFill(ABitMap2 *map, int xs, int ys, int val)
{
	static const short	hdg[4][2] = {
		{0,1}, {1,0}, {0,-1}, {-1,0},
	};
	if ((xs < 0) | (ys < 0))
		return 0;
	if (!map || (xs >= map->Width()) | (ys >= map->Height()))
		return 0;
	const bool	stop = !map->Check(xs, ys);
	if ((val >= 0) & (stop != (val > 0)))	// check we're going right way
		return 0;
	map->Toggle(xs, ys);			// start of flood
	int		cnt = 1;
	IxyPoint	pt(xs, ys);
	IxyStack	stk;
	do {					// loop on search points...
		int		h = 4;
		while (h--) {			// recruit next door neighbor
			xs = pt.x + hdg[h][0];
			ys = pt.y + hdg[h][1];
			if (map->OffBitMap(xs, ys) || map->Check(xs, ys) == stop)
				continue;
			map->Toggle(xs, ys);	// fill this pixel
			++cnt;			// make him a new search point
			stk.Push(IxyPoint(xs,ys));
		}
	} while (stk.Pop(&pt));

	return cnt;				// all done
}

// Draw a polygon as a (thick) outline
int
BdrawPolygon(ABitMap2 *map, const int vxy[][2], int nv, int rad, bool val)
{
	int		np = 0;
	const int *	lastp = vxy[0];

	while (nv-- > 0) {
		np += BfillLine(map, lastp[0], lastp[1],
					vxy[nv][0], vxy[nv][1], rad, val);
		lastp = vxy[nv];
	}
	return np;
}

// Fill a polygon (may have holes connected by seams)
int
BfillPolygon(ABitMap2 *map, const int vxy[][2], int nv, bool val)
{
	int	xymin[2] = {map->Width(), map->Height()};
	int	xymax[2] = {0, 0};
	int	i = nv;
	while (i-- > 0) {			// get extrema
		if (vxy[i][0] < xymin[0]) xymin[0] = vxy[i][0];
		if (vxy[i][0] > xymax[0]) xymax[0] = vxy[i][0];
		if (vxy[i][1] < xymin[1]) xymin[1] = vxy[i][1];
		if (vxy[i][1] > xymax[1]) xymax[1] = vxy[i][1];
	}
	if ((xymin[0] >= xymax[0]) | (xymin[1] >= xymax[1]))
		return 0;
	ABitMap2	pMask(xymax[0]-xymin[0]+1, xymax[1]-xymin[1]+1);
	const int *	lastp = vxy[0];
	for (i = nv; i--; ) {			// draw edges vertically
		int	y = lastp[1];
		int	rise = vxy[i][1] - y;
		if (rise) {
			int	x = lastp[0];
			int	run = vxy[i][0] - x;
			if (rise < 0) {
				x += run; y += rise;
				run = -run; rise = -rise;
			}
			int	xstep = 1;
			if (run < 0) {
				xstep = -1;
				run = -run;
			}
			int	run2=0, rise2=0;
			int	n = rise;
			while (n)
				if (rise2 >= run2) {
					pMask.Toggle(x-xymin[0], y-xymin[1]);
					++y;
					run2 += run;
					--n;
				} else {
					x += xstep;
					rise2 += rise;
				}
		}
		lastp = vxy[i];
	}
	int	np = 0;				// bits toggle fill on each scan
	for (int yy = 0; yy < pMask.Height(); yy++) {
		bool	inPoly = false;
		for (int xx = 0; xx < pMask.Width(); xx++)
			if ((inPoly ^= pMask.Check(xx, yy)))
				np += map->TestAndSet(xymin[0]+xx, xymin[1]+yy, val);
	}
	return np;
}

// Stamp the given color according to a bitmap with optional offset and downsample
int
PstampInk(PanImage *rimp, const ABitMap2 &stencil, PixelVal pv, int xleft, int ytop, const int overs)
{
	if (!rimp | !stencil.Width() | (overs <= 0))
		return -1;
	PImgCond		cnd = rimp->Ready();
	if (cnd == PICread) {
		DMESG(DMCparameter, "Cannot stamp ink on locked image");
		return -1;
	}
						// determine output color space
	const ImgColorSpace *	spc = (cnd >= PICspace) ? rimp->GetCS() :
					pv.csp ? pv.csp : &ICS_Y8;
	pv = PconvertPixel(pv, spc);
	if (cnd != PICready) {			// need to initialize image?
		if (cnd == PICfree) {
			if (!rimp->Init(PICzero))
				return -1;
		} else if (!rimp->Init(xleft+(stencil.Width()+overs-1)/overs,
					ytop+(stencil.Height()+overs-1)/overs,
					*spc, PICzero))
			return -1;
		cnd = PICzero;
	}
	if ((xleft*overs + stencil.Width() <= 0) | (ytop*overs + stencil.Height() <= 0))
		return 0;
	PanImage	myImg;			// working on subimage?
	if (xleft | ytop || (rimp->Width()*overs != stencil.Width()) |
				(rimp->Height()*overs != stencil.Height())) {
		ImgRect	rct;
		rct.xleft = xleft*(xleft > 0);
		rct.ytop = ytop*(ytop > 0);
		rct.xright = xleft + (stencil.Width()+overs-1)/overs;
		if (rct.xright > rimp->Width()) rct.xright = rimp->Width();
		rct.ybottom = ytop + (stencil.Height()+overs-1)/overs;
		if (rct.ybottom > rimp->Height()) rct.ybottom = rimp->Height();
		if (!myImg.Link(rimp, &rct))
			return -1;
		rimp = &myImg;
		xleft *= (xleft < 0);
		ytop *= (ytop < 0);
	}
	if (overs > 1) {			// anti-aliasing bitmap?
		PanImage	overIm(rimp->Width()*overs, rimp->Height()*overs,
								*rimp->GetCS());
		if (!(cnd == PICzero ? overIm.Init(PICzero) :
				PsizeImage(overIm.Img(), rimp->GetImg(), PSnearest)))
			return -1;
						// oversampled stamp (1-level recursion)
		int	np = PstampInk(&overIm, stencil, pv, xleft*overs, ytop*overs);
		if (np <= 0)
			return np;
						// downsample for final result
		if (!PsizeImage(rimp->Img(), overIm.GetImg(), PSbox))
			return -1;
						// all done!
		return (np + overs*overs - 1)/(overs*overs);
	}
	bool	ink = true;			// else 1-1 stencil
	if (cnd == PICzero) {			// zeroed new image?
		if (!pv.csp)			// black on black?
			return (int)stencil.SumTotal();
		if (stencil.SumTotal()*2L > (uint32)rimp->Width()*rimp->Height()) {
			*rimp = pv;		// quicker to reset black bits
			pv = Pblack;
			ink = false;
		}
	}
	int	np = 0;				// draw our pixels
	int	x, y;
	if (!pv.csp & !xleft & !ytop) {		// use faster zeroing method?
		const int	psiz = rimp->PSize();
		for (x = y = 0; stencil.Find(&x, &y, ink) && y < rimp->Height(); ) {
			const int	xstart = x++;
			if (xstart >= rimp->Width()) {
				y++; x = 0;
				continue;
			}
			while (x < rimp->Width() && stencil.Check(x, y) == ink)
				x++;
			memset((*rimp)[y]+xstart*psiz, 0, (x-xstart)*psiz);
			np += x-xstart;
		}
	} else {				// regular pixel writing loop
		for (x = 0, y = -ytop; stencil.Find(&x, &y, ink) && y+ytop < rimp->Height(); x++)
			np += rimp->SetPixel(x+xleft, y+ytop, pv);
	}
	return ink ? np : rimp->Width()*rimp->Height() - np;
}

// Convert image to bitmap(s) based on threshold color
bool
PthreshMap(ABitMap2 map[], const PanImage &im, PixelVal pv, PHistoChan chan)
{
	if (!map | (im.Ready() < PICread))
		return false;
	const int	nc = im.NComp();
	if (nc < 3) {
		if (chan <= PHCblue) {
			static const char	chN[3][6] = {"red","green","blue"};
			DMESGF(DMCparameter, "Cannot access %s channel of gray image",
					chN[chan]);
			return false;
		}
		chan = PHCluminance;
	}
	if (!pv.csp && im.GetCS()->dtype != IDTfloat) {
		DMESG(DMCwarning, "Zero threshold - setting bitmap(s) to true");
		if (!map->NewBitMap(im.Width(), im.Height(), true))
			return false;
		if (chan != PHCrgb3)
			return true;
		return map[1].NewBitMap(im.Width(), im.Height(), true) &&
			map[2].NewBitMap(im.Width(), im.Height(), true);
	}
	if ((chan == PHCluminance) & (nc != 1)) {
		ImgColorSpace	mySpace;	// convert input to Y
		im.GetCS(&mySpace);
		mySpace.format = IPFy;
		if (pv.csp->dtype < mySpace.dtype)
			mySpace.dtype = pv.csp->dtype;
		PanImage	myImg(mySpace);
		return myImg.Load(im) && PthreshMap(map, myImg, pv, PHCluminance);
	}
	if (!map->NewBitMap(im.Width(), im.Height()))
		return false;
	if (chan == PHCrgb3 && 
			(!map[1].NewBitMap(im.Width(), im.Height()) ||
				!map[2].NewBitMap(im.Width(), im.Height())))
		return false;
	pv = PconvertPixel(pv, im.GetCS());
	int	x, y;
	switch (im.GetCS()->dtype) {
	case IDTubyte:
		for (y = 0; y < im.Height(); y++) {
			const uby8 *	bp = im.GetRowByte(y);
			for (x = 0; x < im.Width(); x++, bp += nc) {
				if (chan == PHCrgb) {
					if (bp[0] + bp[1] + bp[2] >=
						pv.v.b[0] + pv.v.b[1] + pv.v.b[2])
						map->Set(x, y);
					continue;
				}
				if (chan <= PHCblue) {
					if (bp[chan] >= pv.v.b[chan])
						map->Set(x, y);
					continue;
				}
				if (bp[0] >= pv.v.b[0])
					map->Set(x, y);
				if (chan == PHCrgb3) {
					if (bp[1] >= pv.v.b[1])
						map[1].Set(x, y);
					if (bp[2] >= pv.v.b[2])
						map[2].Set(x, y);
				}
			}
		}
		break;
	case IDTfloat:
		for (y = 0; y < im.Height(); y++) {
			const float *	fp = im.GetRowFloat(y);
			for (x = 0; x < im.Width(); x++, fp += nc) {
				if (chan == PHCrgb) {
					if (fp[0] + fp[1] + fp[2] >=
						pv.v.f[0] + pv.v.f[1] + pv.v.f[2])
						map->Set(x, y);
					continue;
				}
				if (chan <= PHCblue) {
					if (fp[chan] >= pv.v.f[chan])
						map->Set(x, y);
					continue;
				}
				if (fp[0] >= pv.v.f[0])
					map->Set(x, y);
				if (chan == PHCrgb3) {
					if (fp[1] >= pv.v.f[1])
						map[1].Set(x, y);
					if (fp[2] >= pv.v.f[2])
						map[2].Set(x, y);
				}
			}
		}
		break;
	case IDTushort:
		for (y = 0; y < im.Height(); y++) {
			const unsigned short *	sp = im.GetRowShort(y);
			for (x = 0; x < im.Width(); x++, sp += nc) {
				if (chan == PHCrgb) {
					if (sp[0] + sp[1] + sp[2] >=
						pv.v.s[0] + pv.v.s[1] + pv.v.s[2])
						map->Set(x, y);
					continue;
				}
				if (chan <= PHCblue) {
					if (sp[chan] >= pv.v.s[chan])
						map->Set(x, y);
					continue;
				}
				if (sp[0] >= pv.v.s[0])
					map->Set(x, y);
				if (chan == PHCrgb3) {
					if (sp[1] >= pv.v.s[1])
						map[1].Set(x, y);
					if (sp[2] >= pv.v.s[2])
						map[2].Set(x, y);
				}
			}
		}
		break;
	default:
		DMESG(DMCparameter, "Unsupported pixel data type");
		return false;
	}
	return true;
}

// Stencil operator to overlay selected image pixels, no offset
bool
PstencilImage(PanImage *rimp, const ABitMap2 &stencil, const PanImage &sim)
{
	if (!rimp | (sim.Ready() < PICread))
		return false;
	if ((stencil.Width() != sim.Width()) | (stencil.Height() != sim.Height())) {
		DMESG(DMCparameter, "Stencil and source image dimensions must match");
		return false;
	}
	int		x=0, y=0;		// check for empty bitmap
	const bool	empty = !stencil.Find(&x, &y);

	switch (rimp->Ready()) {		// initialize / finish easy cases
	case PICready:
		if (rimp->Compat(sim, PICMok)) {
			DMESG(DMCparameter,
				"Source and destination image dimensions must match");
			return false;
		}
		break;
	case PICfree:				// force dimensions to agree
		rimp->Init(sim.Width(), sim.Height(), *rimp->GetCS(), PICfree);
		// fall through...
	case PICspace:
		if (empty || stencil.SumTotal()*2 <=
				(uint32)stencil.Width()*stencil.Height()) {
			if (!rimp->Init(sim.Width(), sim.Height(),
						*rimp->GetCS(), PICzero))
				return false;
			break;			// save on conversion time...
		}
		// falling through...
	case PICvoid:
		if (empty)			// empty in source color space?
			return rimp->Init(sim.Width(), sim.Height(),
						*sim.GetCS(), PICzero);
		if (!rimp->Load(sim))		// else load original...
			return false;
		*rimp *= stencil;		// and set "false" to black
		return true;
	case PICread:
		DMESG(DMCparameter, "Cannot apply stencil to locked image");
		return false;
	default:				// XXX should never get here
		return false;
	}
	if (empty)				// nothing to stencil?
		return true;
						// stencil/convert pixels
	if (rimp->Compat(sim, PICMall, false)) {
		const int	psiz = sim.PSize();
		for (x = y = 0; stencil.Find(&x, &y); ) {
			const int	xstart = x++;
			while (stencil.Check(x, y))
				x++;
			memcpy((*rimp)[y]+xstart*psiz, sim[y]+xstart*psiz,
					(x-xstart)*psiz);
		}
	} else {				// convert color as we go
		const void *	ccvp = PgetColorConv(rimp->GetCS(), sim.GetCS());
		if (!ccvp)
			return false;
		for (x = y = 0; stencil.Find(&x, &y); ) {
			const int	xstart = x++;
			while (stencil.Check(x, y))
				x++;
			if (!PmapPixels(rimp->P(xstart,y), sim.GetP(xstart,y),
					x-xstart, ccvp))
				return false;
		}
	}
	return true;
}

// Perform operation on selected image pixels; compare with PblendImageOpCB()
bool
PstencilImageOpCB(PanImage *idstp, const ABitMap2 mask, PimageOp *op, void *udp)
{
	const int	fmarg = 4;		// margin for filtering

	if (!idstp | !op | (mask.Width() <= 0) | (mask.Height() <= 0))
		return false;
	if (idstp->Ready() != PICready) {
		DMESGF(DMCparameter, "Cannot operate on %s image",
			idstp->ReadOnly() ? "read-only" : "uninitialized");
		return false;
	}
	if ((idstp->Width() != mask.Width()) | (idstp->Height() != mask.Height())) {
		DMESG(DMCparameter, "Mismatched image dimensions");
		return false;
	}
	/*
	 * Optimize operation by trimming mask.
	 * This way, a small patch somewhere in the image only
	 * operates on that region.
	 */
	int	xymin[2], wh[2];
	if (!mask.GetBoundRect(xymin, wh))
		return true;			// nothing to do, but OK
	ImgRect		subRect;		// else link subregion
	subRect.xleft = xymin[0] - fmarg;
	subRect.xleft *= (subRect.xleft > 0);
	subRect.ytop = xymin[1] - fmarg;
	subRect.ytop *= (subRect.ytop > 0);
	if ((subRect.xright = xymin[0] + wh[0] + fmarg) > idstp->Width())
		subRect.xright = idstp->Width();
	if ((subRect.ybottom = xymin[1] + wh[1] + fmarg) > idstp->Height())
		subRect.ybottom = idstp->Height();
	PanImage	iSub, iWork;
	if (!iSub.Link(idstp, &subRect))
		return false;			// should never happen
	double		frac = PrectArea(&subRect) /
					double(idstp->Width()*idstp->Height());
	DTESTF((frac < .999), DMCtrace, "Operating on %.1f%% submimage", 100.*frac);
						// op initializes destination
	return ((*op)(iWork.Img(), iSub.GetImg(), udp) > 0 &&
			PstencilImage(&iSub,
				mask.GetRect(subRect.xleft, subRect.ytop,
						subRect.xright - subRect.xleft,
						subRect.ybottom - subRect.ytop),
							iWork));
}

// Resize a bitmap -- use a lower threshold to be more conservative with 1's
bool
BresizeBitmap(ABitMap2 *map, const ABitMap2 &src, float thresh)
{
	if (!map || !map->Width() | !map->Height())
		return false;
	if (!src.Width() | !src.Height())
		return false;
	if ((map->Width() == src.Width()) & (map->Height() == src.Height())) {
		*map = src;
		return true;
	}
	PixelVal	pv;
	pv.csp = &ICS_Y8;
	pv.v.b[0] = 255;
	PanImage	origIm;
	int		np = PstampInk(&origIm, src, pv);
	if (np <= 0) {
		map->ClearBitMap(false);
		return (np == 0);
	}
	if (np >= src.Width()*src.Height()) {
		map->ClearBitMap(true);
		return true;
	}
	PanImage	resIm(map->Width(), map->Height(), *pv.csp);
	if (!PsizeImage(resIm.Img(), origIm.GetImg(), PSlinear))
		return false;

	pv.v.b[0] = (thresh < 1.5f/256.f) ? 1 :
		(thresh > 254.5f/256.f) ? 254 :
				int(thresh * 256.f);

	return PthreshMap(map, resIm, pv);
}
