/*
 *  pdequant.cpp
 *  pan
 *
 *  Dequantize image using nearest greater/smaller valued pixel primary
 *  Interpolate linearly between to generate smooth gradients
 *
 *  Created by Greg Ward on 3/2/15.
 *  Copyright 2015 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "dmessage.h"
#include "panimage.h"
#include "abitmap.h"
#include "astack.h"

// return maximum component value for image color space
static inline int
maxIntComponent(const PanImage &im)
{
	if (im.GetCS()->dtype == IDTubyte)
		return 0xff;

	return 0xffff;
}

// return integer component value (no checks!)
static inline int
getIntComponent(const PanImage &im, int x, int y, int c)
{
	const int		nc = im.NComp();
	
	if (im.GetCS()->dtype == IDTubyte)
		return im[y][x*nc + c];

	return ((unsigned short *)im[y])[x*nc + c];
}

// Candidate
struct Candidate {
	short		gx, gy;		// candidate grid position
};

// Class to find nearest neighbor with component above or below reference
class NeighborFinder {
	bool		above;		// looking for next larger value?
	PanImage	origIm;		// read-only link to image
	int		comp;		// component we're interested in
	PanImage	exIm;		// extremum image
	int		subsmp;		// subdivision rate
	Candidate *	ca;		// our candidate list
	int		ncands;		// number of valid candidates
	long		curmaxD2;	// distance^2 with guaranteed outlier
	int		candLen;	// current candidate array length
	int		rx, ry;		// reference position
	int		rval;		// reference value
	int		vlimit;		// min or max over entire image
protected:
	long		GetMinimum2(int gx, int gy) const {
				gx -= rx/subsmp;
				if (gx < 0) gx = -gx-1;
				else if (gx > 0) --gx;
				gy -= ry/subsmp;
				if (gy < 0) gy = -gy-1;
				else if (gy > 0) --gy;
				return (long)subsmp*subsmp*(gx*gx + gy*gy);
			}
	long		GetMaximum2(int gx, int gy) const {
				gx -= rx/subsmp;
				if (gx < 0) gx = 1-gx;
				else ++gx;
				gy -= ry/subsmp;
				if (gy < 0) gy = 1-gy;
				else ++gy;
				return (long)subsmp*subsmp*(gx*gx + gy*gy);
			}
	bool		ExceedsAt(int gx, int gy) const {
				int	v;
				if (exIm.GetCS()->dtype == IDTubyte)
					v = exIm.GetRowByte(gy)[gx];
				else
					v = exIm.GetRowShort(gy)[gx];
				if (above) return (v > rval);
				return (v < rval);
			}
	bool		CorralCandidates(int rad);
	int		CullCandidates();
	int		AddCandidate(int gx, int gy);
	void		CandidateRect(ImgRect *r, int n) const {
				r->xright = r->xleft = ca[n].gx*subsmp;
				if ((r->xright += subsmp) > origIm.Width())
					r->xright = origIm.Width();
				r->ybottom = r->ytop = ca[n].gy*subsmp;
				if ((r->ybottom += subsmp) > origIm.Height())
					r->ybottom = origIm.Height();
			}
	long		FindInCandidate(int xy[2], int n) const;
public:
			NeighborFinder(const PanImage &orig, int c,
						bool upwards, int div = 5);
			~NeighborFinder() {
				if (ca) free(ca);
			}
	const PanImage &OrigImage() const {
				return origIm;
			}
	bool		LookingUp() const {
				return above;
			}
	int		Component() const {
				return comp;
			}
	int		Decimation() const {
				return subsmp;
			}
	long		NewReference(int x, int y, int v = -1) {
				if (v < 0) v = getIntComponent(origIm,x,y,comp);
				if (v == vlimit) return 0;
				if (v == rval && x/subsmp == rx/subsmp &&
						y/subsmp == ry/subsmp) {
					rx = x; ry = y;
					return curmaxD2*(ncands > 0);
				}
				rx = x; ry = y; rval = v;
				ncands = 0;
				curmaxD2 = ((long)origIm.Width()*origIm.Width() +
					    (long)origIm.Height()*origIm.Height());
				int	r = 0;
				while (CorralCandidates(r)) r++;
				CullCandidates();
				return curmaxD2*(ncands > 0);
			}
	int		GetReference(int xy[2]) const {
				xy[0] = rx; xy[1] = ry;
				return rval;
			}
	long		FindNearest2(int nxy[2] = 0) const {
				long	nearest2 = curmaxD2 + 1;
				for (int i = ncands; i-- > 0; ) {
					int	txy[2];
					long	test2 = FindInCandidate(txy, i);
					if (test2 < 0) return -1;
					if (!test2 | (test2 >= nearest2))
						continue;
					if (nxy) {
						nxy[0] = txy[0]; nxy[1] = txy[1];
					}
					nearest2 = test2;
				}
				return nearest2*(nearest2 <= curmaxD2);
			}
};

// Element of sorted search array for block
struct BlockPixel {
	uby8		bx, by;		// offset within block
	unsigned short	v;		// pixel component value
};

// Search neighborhood in efficient order for NeighborFinder object
class BlockSearch {
	PanImage	origIm;		// read-only link to original image
	int		comp;		// which component do we care about?
	int		bsiz;		// block side length (decimation)
	int		gx, gy;		// current block's grid position
	BlockPixel *	bpa;		// block pixel array (sorted)
	int		bw, bh;		// block width and height
	int		ndx;		// current index in block
public:
			BlockSearch(const NeighborFinder &nf) {
				origIm.LinkRO(nf.OrigImage());
				comp = nf.Component();
				bsiz = nf.Decimation();
				bpa = new BlockPixel [bsiz*bsiz];
				InitSearch();
			}
			~BlockSearch() {
				delete [] bpa;
			}
	bool		InitSearch() {
				if (origIm.Ready() < PICread) return false;
				gy = 0; gx = -1;
				bw = bh = 0; ndx = 0;
				return true;
			}
	int		NextPixel(int xy[2]);
};

// Set up neighbor finder for nearest above/below value
NeighborFinder::NeighborFinder(const PanImage &orig, int c,
					bool upwards, int div)
{
	ca = NULL; ncands = candLen = 0;
	rx = ry = -1; rval = -1;
	if (!origIm.LinkRO(orig))
		return;
	DASSERT(orig.GetCS()->dtype == IDTubyte || orig.GetCS()->dtype == IDTushort);
	origIm.ClearInfo();
	const int	nc = orig.NComp();
	DASSERT((0 <= c) & (c < nc));
	DASSERT(div > 0);
	above = upwards;
	comp = c;
	subsmp = div;
	ImgColorSpace	mySpace;
	orig.GetCS(&mySpace);
	mySpace.format = IPFy;			 // compute local min or max
	if (!exIm.Init((orig.Width()+subsmp-1)/subsmp, (orig.Height()+subsmp-1)/subsmp, mySpace))
		return;
	int		x, y;			// do min or max search
	if (above) {
		exIm = .0f;
		vlimit = 0;
	} else {
		exIm = 1.f;
		vlimit = maxIntComponent(origIm);
	}
	if (orig.GetCS()->dtype == IDTubyte) {
		for (y = orig.Height(); y--; ) {
			const uby8 *	obp = orig.GetRowByte(y) + comp;
			uby8 *		gbp = exIm.RowByte(y/subsmp);
			for (x = 0; x < orig.Width(); obp += nc,
							gbp += !(++x % subsmp))
				if ((*obp < *gbp) ^ above)
					*gbp = *obp;
		}
		for (y = exIm.Height(); y--; ) {
			const uby8 *	gbp = exIm.GetRowByte(y);
			for (x = 0; x < exIm.Width(); x++, gbp++)
				if ((*gbp < vlimit) ^ above)
					vlimit = *gbp;
		}
	} else {
		for (y = orig.Height(); y--; ) {
			const unsigned short *	osp = orig.GetRowShort(y) + comp;
			unsigned short *	gsp = exIm.RowShort(y/subsmp);
			for (x = 0; x < orig.Width(); osp += nc,
							gsp += !(++x % subsmp))
				if ((*osp < *gsp) ^ above)
					*gsp = *osp;
		}
		for (y = exIm.Height(); y--; ) {
			const unsigned short *	gsp = exIm.GetRowShort(y);
			for (x = 0; x < exIm.Width(); x++, gsp++)
				if ((*gsp < vlimit) ^ above)
					vlimit = *gsp;
		}
	}
	exIm.SetReadOnly();			// now tied to component image
}

// Add a new candidate to our list
int
NeighborFinder::AddCandidate(int gx, int gy)
{
	if (ca == NULL)
		ca = (Candidate *)malloc(sizeof(Candidate)*(candLen = 64));
	else if (ncands >= candLen)
		ca = (Candidate *)realloc(ca, sizeof(Candidate)*(candLen += 64));
	if (ca == NULL) {
		DMESG(DMCmemory, "malloc() failed in AddCandidate()");
		return ncands = candLen = 0;
	}
	ca[ncands].gx = gx; ca[ncands].gy = gy;
	return ++ncands;
}

// Find candidates at the given square radius
bool
NeighborFinder::CorralCandidates(int rad)
{
	if (exIm.Ready() < PICread)
		return false;

	const int	gxc = rx/subsmp;
	const int	gyc = ry/subsmp;
	bool		inRange = false;
	bool		possible;
	int		gx, gy;
	long		max2;

	for (gy = gyc-rad; gy <= gyc+rad; gy++) {
	    if (gy < 0) continue;
	    if (gy >= exIm.Height()) break;
	    const int	gxstep = (gyc-rad < gy) & (gy < gyc+rad) ? 2*rad : 1;
	    for (gx = gxc-rad; gx <= gxc+rad; gx += gxstep) {
		if (gx < 0) continue;
		if (gx >= exIm.Width()) break;
		inRange |= possible = (GetMinimum2(gx, gy) <= curmaxD2);
		if (possible && ExceedsAt(gx, gy)) {
			if (!AddCandidate(gx, gy))
				return false;
			if ((max2 = GetMaximum2(gx, gy)) < curmaxD2)
				curmaxD2 = max2;
		}
	    }
	}
	return inRange;		// return true if there might still be more
}

// Remove candidates that are clearly not going to be nearest the reference
int
NeighborFinder::CullCandidates()
{
	int	ngood = 0;
	int	n;

	for (n = 0; n < ncands; n++) {
		if (GetMinimum2(ca[n].gx, ca[n].gy) > curmaxD2)
			continue;	// not useful
		if (ngood < n)
			ca[ngood] = ca[n];
		++ngood;
	}
	n -= ngood;
	ncands = ngood;
	return n;			// return number culled
}

// Find nearest exceedant in the given candidate block
long
NeighborFinder::FindInCandidate(int xy[2], int n) const
{
	ImgRect		myRect;
	PanImage	subIm;

	CandidateRect(&myRect, n);
	if (!subIm.LinkRO(origIm, &myRect))
		return -1;
	long		nearest2 = curmaxD2 + 1;
	for (int y = subIm.Height(); y--; ) {
	    const int		dy = (myRect.ytop + y) - ry;
	    for (int x = 0; x < subIm.Width(); x++) {
		const int	dx = (myRect.xleft + x) - rx;
		const long	d2 = (long)dx*dx + (long)dy*dy;
		if (d2 >= nearest2)
			continue;
		int		v = getIntComponent(subIm, x, y, comp);
		if ((v != rval) & ((v < rval) ^ above)) {
			xy[0] = myRect.xleft + x;
			xy[1] = myRect.ytop + y;
			nearest2 = d2;
		}
	    }
	}
	return nearest2*(nearest2 <= curmaxD2);
}

// Sort block position by component value
static int
cmpBP(const void *p1, const void *p2)
{
	return (int)(*(const BlockPixel *)p1).v -
			(int)(*(const BlockPixel *)p2).v;
}

// Go to next pixel in our neighborhood search
int
BlockSearch::NextPixel(int xy[2])
{
	if (origIm.Ready() < PICread)
		return -1;
	if (ndx < bw*bh) {		// still in block?
		xy[0] = gx*bsiz + bpa[ndx].bx;
		xy[1] = gy*bsiz + bpa[ndx].by;
		return bpa[ndx++].v;	// advance
	}
					// else advance to & sort next block
	if (++gx*bsiz >= origIm.Width()) {
		if (++gy*bsiz >= origIm.Height())
			return -1;	// all done!
		gx = 0;
	}
	bw = bsiz;			// determine block dimensions
	if (gx*bsiz + bw > origIm.Width())
		bw = origIm.Width() - gx*bsiz;
	bh = bsiz;
	if (gy*bsiz + bh > origIm.Height())
		bh = origIm.Height() - gy*bsiz;
	for (int j = bh; j--; )		// assign pixels
		for (int i = bw; i--; ) {
			BlockPixel &	blk = bpa[j*bw + i];
			blk.bx = i;
			blk.by = j;
			blk.v = getIntComponent(origIm, gx*bsiz+i, gy*bsiz+j, comp);
		}
					// sort by value to optimize queries
	qsort(bpa, bw*bh, sizeof(BlockPixel), cmpBP);
	ndx = 0;			// return first entry & advance
	xy[0] = gx*bsiz + bpa[0].bx;
	xy[1] = gy*bsiz + bpa[0].by;
	return bpa[ndx++].v;
}

// dequantize image component
static bool
dequantizeComponent(PanImage *dstp, const PanImage &isrc, const int c,
			PquantMethod *qm, void *udp)
{
	const int	bsiz = 5;
	const int	nc = dstp->NComp();
	const float	normf = 1.f/(maxIntComponent(isrc) + 1.f);
	NeighborFinder	belowNF(isrc, c, false, bsiz);
	NeighborFinder	aboveNF(isrc, c, true, bsiz);
	BlockSearch	traverse(belowNF);
	int		pixy[2], refv;

	while ((refv = traverse.NextPixel(pixy)) >= 0) {
		float *	dp = dstp->RowFloat(pixy[1]) + pixy[0]*nc + c;
		*dp = normf*(refv + .5f);	// default value
		int		quant = 1;
		if (qm != NULL) {		// quantization call-back
			quant = (*qm)(pixy[0], pixy[1], c, udp);
			if (quant <= 0)
				continue;
		}
		long		aboveR2 = aboveNF.NewReference(pixy[0], pixy[1], refv);
		if (!aboveR2) {			// maximum plateau?
			*dp -= .5f*quant;	// match surrounding potential
			continue;
		}
		long		belowR2 = belowNF.NewReference(pixy[0], pixy[1], refv);
		if (!belowR2) {			// minimum valley?
			*dp += .5f*quant;	// match surrounding potential
			continue;
		}
		if ((belowR2 <= bsiz*bsiz) & (aboveR2 <= bsiz*bsiz))
			continue;		// don't bother to correct
		int	belowPos[2], abovePos[2];
		belowR2 = belowNF.FindNearest2(belowPos);
		aboveR2 = aboveNF.FindNearest2(abovePos);
		if ((belowR2 <= bsiz*bsiz) & (aboveR2 <= bsiz*bsiz))
			continue;		// don't bother to adjust
		if (quant > 1) {		// else check quantization
			int	q = abs(refv - (belowR2 <= aboveR2 ?
					getIntComponent(isrc,
						belowPos[0],belowPos[1],c) :
					getIntComponent(isrc,
						abovePos[0],abovePos[1],c)) );
			if (q < quant)
				quant = q;	// neighbor is nearer in value
		}
						// interpolate linearly
		const float	belowR = (float)sqrt((double)belowR2);
		const float	aboveR = (float)sqrt((double)aboveR2);
		*dp = normf*( ((refv-.5f*quant)*aboveR + (refv+.5f*quant)*belowR)
					/ (aboveR + belowR) + .5f );
	}
	return true;
}

// Callback using integer image to look up quantization value
int
PQmethodImage(int x, int y, int c, void *udp)
{
	const ImgStruct *	ip = (const ImgStruct *)udp;

	if (ip == NULL || (ip->csp == NULL) | (ip->img == NULL))
		return 1;
	if ((x < 0) | (x >= ip->xres) | (y < 0) | (y >= ip->yres))
		return 1;
	const int	nc = ImgPixelLen[ip->csp->format];
	if (c >= nc)
		c = 0;			// allows grayscale for RGB input
	switch (ip->csp->dtype) {
	case IDTubyte:
		return ProwPtr(ip,y)[x*nc + c];
	case IDTushort:
		return ((unsigned short *)ProwPtr(ip,y))[x*nc + c];
	}
	return 1;
}

// Simple call-back to return constant value from integer's address
int
PQmethodConstant(int x, int y, int c, void *udp)
{
	if (udp == NULL)
		return 1;
	return *(int *)udp;
}

// Dequantize image (using optional quantization from call-back at each pixel)
bool
Pdequantize(PanImage *dstp, const PanImage &isrc, float sf,
			PquantMethod *qm, void *udp)
{
	if ((dstp == NULL) || !dstp->Ready() | (isrc.Ready() < PICread))
		return false;
	const bool	inp8bit = (isrc.GetCS()->dtype == IDTubyte);
	if (!inp8bit && isrc.GetCS()->dtype != IDTushort) {
		DMESG(DMCparameter, "Input must be 8- or 16-bit to dequantize");
		return false;
	}
	ImgColorSpace	destCS;
	PanImage	inpIm;
	ImgColorSpace	ourCS;
	int		c;
	dstp->GetCS(&destCS);
	inpIm.LinkRO(isrc);
	if (PimagesOverlap(dstp->GetImg(), isrc.GetImg())) {
		inpIm.GetCS(&destCS);
		PrealCS(&destCS);		// assumed destination space
	} else if (destCS.dtype != IDTfloat) {
		DMESG(DMCparameter, "Destination image must be float");
		return false;
	} else if (destCS.format != inpIm.GetCS()->format) {
		DMESG(DMCparameter, "Mismatched image formats");
		return false;
	}
	inpIm.GetCS(&ourCS);			// prepare dequantized holder
	ourCS.dtype = IDTfloat;
	if (!dstp->Init(inpIm.Width(), inpIm.Height(), ourCS))
		return false;
	if (inpIm.GetInfo())			// copy metadata
		inpIm.GetInfo(dstp->Info());
	for (c = inpIm.NComp(); c--; )		// dequantize each component
		if (!dequantizeComponent(dstp, inpIm, c, qm, udp))
			return false;
						// destination color conversion
	return dstp->ConvertCS(destCS, sf);
}

// Same-valued flood-fill (no recursion)
static unsigned long
floodEqual(ABitMap2 *nwr, const PanImage &isrc, int x, int y, const int c)
{
	static const short	hdg[4][2] = {
		{0,1}, {1,0}, {0,-1}, {-1,0},
	};
	const int	rval = getIntComponent(isrc, x, y, c);
						// initialize destination
	nwr->NewBitMap(isrc.Width(), isrc.Height());
	nwr->Set(x, y);				// one pixel to start
	unsigned long	np = 1;
	IxyPoint	pt(x,y);
	IxyStack	stk;
	do {					// loop on candidates...
		int		h = 4;
		while (h--) {			// recruit next door neighbor
			x = pt.x + hdg[h][0];
			y = pt.y + hdg[h][1];
			if (nwr->OffBitMap(x, y) || nwr->Check(x, y))
				continue;
			if (getIntComponent(isrc, x, y, c) != rval)
				continue;
			nwr->Set(x, y);		// add this pixel
			++np;			// make him a new candidate
			stk.Push(IxyPoint(x,y));
		}
	} while (stk.Pop(&pt));

	return np;				// all done
}

// Fill the given region in quantization image (no checks for speed)
static void
fillRegion(PanImage *dstp, const ABitMap2 &fillMap, const int c, const int val)
{
	ImgStruct *	ip = dstp->Img();
	const int	nc = dstp->NComp();
	int		x, y;

	if (dstp->GetCS()->dtype == IDTubyte) {
		for (x = y = 0; fillMap.Find(&x, &y); x++)
			ProwPtr(ip,y)[x*nc + c] = val;
		return;
	}
	for (x = y = 0; fillMap.Find(&x, &y); x++)
		((unsigned short *)ProwPtr(ip,y))[x*nc + c] = val;
}

// Estimate quantization map for image based on neighbor values
bool
Pquantization(PanImage *dstp, const PanImage &isrc)
{
	static const int	oxy[4][2] = {-1,0, 1,0, 0,-1, 0,1};
	if ((dstp == NULL) | (isrc.Ready() < PICread))
		return false;
	if (isrc.GetCS()->dtype != IDTubyte && isrc.GetCS()->dtype != IDTushort) {
		DMESG(DMCparameter, "Input must be 8- or 16-bit for quantization");
		return false;
	}
	if (!dstp->Init(isrc.Width(), isrc.Height(), *isrc.GetCS(), PICzero))
		return false;
	for (int c = isrc.NComp(); c--; ) {	// each component independent
	    ABitMap2	doneMap(isrc.Width(), isrc.Height());
	    int		xs, ys;
	    for (xs = ys = 0; doneMap.Find(&xs, &ys, false); xs++) {
		const int	rval = getIntComponent(isrc, xs, ys, c);
		ABitMap2	fillMap;
		unsigned long	nfilled = floodEqual(&fillMap, isrc, xs, ys, c);
		doneMap |= fillMap;
		if (nfilled < 50)
		    continue;			// below threshold of interest
		unsigned long	nperim = 0;	// else check perimeter
		int		closestAbove = 0x10000;
		int		closestBelow = -1;
		int		xt, yt;		// find nearest-valued neighbors
		for (xt = yt = 0; fillMap.Find(&xt, &yt); xt++) {
		    for (int i = 4; i--; ) {
			const int	bx = xt + oxy[i][0];
			const int	by = yt + oxy[i][1];
			if (fillMap.OffBitMap(bx, by) || fillMap.Check(bx, by))
			    continue;
			const int	v = getIntComponent(isrc, bx, by, c);
			if (v > rval) {
				if (v < closestAbove) closestAbove = v;
			} else {
				if (v > closestBelow) closestBelow = v;
			}
			++nperim;
		    }
		}
		if (nfilled-nperim <= nperim)
		    continue;			// not enough interior
						// else take larger gap
		int	gap = closestAbove - rval;
		if (closestAbove > 0xffff) {
		    if (closestBelow < 0)
			continue;		// single-valued image?
		    gap = rval - closestBelow;
		} else if (closestBelow >= 0 &&
			    rval - closestBelow > gap &&
			    rval - closestBelow <= 2*gap + 1)
		    gap = rval - closestBelow;
						// fill with measured quantum
		fillRegion(dstp, fillMap, c, gap);
	    }
	}
	return true;				// quantization map ready
}

#if 0			// original (slow) method -- takes about 50 times longer

// fill the circle from the given center to circumference point
static void
fillBitCircle(ABitMap2 *bp, const int xc, const int yc, int rad2)
{
	const int	rad = int(sqrt(rad2) + .5);

	rad2 = rad*rad;				// ensure perfect square
	for (int y = yc-rad; y <= yc+rad; y++) {
		const int	yr2 = (y-yc)*(y-yc);
		const int	spanRad = int(sqrt(rad2 - yr2) + .99);
		for (int x = xc-spanRad; x <= xc+spanRad; x++)
			if ((x-xc)*(x-xc) + yr2 <= rad2)
				bp->Set(x, y);
	}
}

// find nearest pixel position whose value is above/below the given reference
static int
getNearest(int xyNear[2], const PanImage &isrc, const int rx, const int ry,
			const int c, const bool above)
{
	const int	refv = getIntComponent(isrc, rx, ry, c);
	const int	maxr = isrc.Width() <= isrc.Height() ?
					isrc.Width()>>1 : isrc.Height()>>1 ;
	int		r, r2, lastr2=0;
	int		x, y;
						// widening search radius
	for (r = 1; r < maxr; r++) {
		r2 = r*r;
		for (y = ry-r; y <= ry+r; y++) {
			if (y < 0) continue;
			if (y >= isrc.Height()) break;
			const int	yr2 = (y-ry)*(y-ry);
			const int	spanr = int(sqrt(r2 - yr2) + .99);
			for (x = rx-spanr; x <= rx+spanr; x++) {
				if (x < 0) continue;
				if (x >= isrc.Width()) break;
				const int	curr2 = (x-rx)*(x-rx) + yr2;
				if (curr2 <= lastr2) {
					if (x < rx) x = rx + (rx - x);
					continue;
				}
				if (curr2 > r2)
					continue;
				int	v = getIntComponent(isrc, x, y, c);
				if ((v != refv) & ((v < refv) ^ above)) {
					xyNear[0] = x; xyNear[1] = y;
					return curr2;
				}
			}
		}
		lastr2 = r2;
	}
	return 0;				// failure!
}

// move active reference to indicated position and update nearest
static int
updateNearest(int xyNear[2], ABitMap2 *bp, const PanImage &isrc,
		const int rx, const int ry, const int c, const bool above)
{
	const int	begr2 = (xyNear[0] - rx)*(xyNear[0] - rx) +
					(xyNear[1] - ry)*(xyNear[1] - ry);
	if (begr2 <= 2) {
		bp->NewBitMap(0,0);		// below mapping threshold
		return 1;
	}
	const int	refv = getIntComponent(isrc, rx, ry, c);
	ABitMap2	checkMap = *bp;		// remember previous map
	const ABitMap2 *cmp = &checkMap;
	if (!bp->NewBitMap(isrc.Width(), isrc.Height()))
		return 0;
	fillBitCircle(bp, rx, ry, begr2);
	if (checkMap.Width()) {			// update to previous map?
		checkMap.Invert();
		checkMap &= *bp;
	} else					// else check all neighbors
		cmp = bp;
						// search for closer
	int	curr2 = begr2;
	int	x, y;
	for (x = y = 0; cmp->Find(&x, &y); x++) {
		const int	tstr2 = (x-rx)*(x-rx) + (y-ry)*(y-ry);
		if (tstr2 >= curr2)
			continue;
		int	v = getIntComponent(isrc, x, y, c);
		if ((v == refv) | ((v > refv) ^ above))
			continue;
		xyNear[0] = x;
		xyNear[1] = y;
		curr2 = tstr2;
	}
	if (curr2 <= 2) {			// below mapping threshold?
		bp->NewBitMap(0,0);
	} else if (curr2 < begr2) {		// else need to update map?
		bp->ClearBitMap();
		fillBitCircle(bp, rx, ry, curr2);
	}
	return curr2;
}

// dequantize image component
static bool
dequantizeComponent(PanImage *dstp, const PanImage &isrc, const int c,
					const unsigned short extrema[2])
{
	const int	nc = dstp->NComp();
	const float	normf = 1.f/(maxIntComponent(isrc) + 1.f);
	ABitMap2	belowMap, aboveMap;
	int		nearBelow[2], belowR2, nearAbove[2], aboveR2;

	for (int y = 0; y < isrc.Height(); y++) {
	    float *	dp = dstp->RowFloat(y) + c;
	    int		last_refv = -1;
	    for (int x = 0; x < isrc.Width(); x++, dp += nc) {
		const int	refv = getIntComponent(isrc, x, y, c);
		*dp = normf*(refv + 0.5f);	// in case we skip this one
		if (refv == last_refv) {
			belowR2 = updateNearest(nearBelow, &belowMap,
							isrc, x, y, c, false);
			aboveR2 = updateNearest(nearAbove, &aboveMap,
							isrc, x, y, c, true);
		} else {
			last_refv = -1;		// starting search over
			if ((refv == extrema[0]) | (refv == extrema[1]))
				continue;
			belowR2 = getNearest(nearBelow, isrc, x, y, c, false);
			if (!belowR2)
				continue;
			aboveR2 = getNearest(nearAbove, isrc, x, y, c, true);
			if (!aboveR2)
				continue;
			if (belowR2 > 2) {
				if (!belowMap.NewBitMap(isrc.Width(), isrc.Height()))
					return false;
				fillBitCircle(&belowMap, x, y, belowR2);
			} else
				belowMap.NewBitMap(0,0);
			if (aboveR2 > 2) {
				if (!aboveMap.NewBitMap(isrc.Width(), isrc.Height()))
					return false;
				fillBitCircle(&aboveMap, x, y, aboveR2);
			} else
				aboveMap.NewBitMap(0,0);
			last_refv = refv;
		}
		if ((belowR2 <= 2) & (aboveR2 <= 2))
			continue;		// close enough
						// else interpolate linearly
		const float	belowR = (float)sqrt((double)belowR2);
		const float	aboveR = (float)sqrt((double)aboveR2);
		*dp = normf*( ((refv-.5f)*aboveR + (refv+.5f)*belowR)
					/ (aboveR + belowR) + .5f );
	    }
	}
	return true;
}

// Dequantize image
bool
Pdequantize(PanImage *dstp, const PanImage &isrc, float sf)
{
	if ((dstp == NULL) || !dstp->Ready() | (isrc.Ready() < PICread))
		return false;
	const bool	inp8bit = (isrc.GetCS()->dtype == IDTubyte);
	if (!inp8bit && isrc.GetCS()->dtype != IDTushort) {
		DMESG(DMCparameter, "Input must be 8- or 16-bit to dequantize");
		return false;
	}
	uby8		b_extrema[3][2];
	unsigned short	s_extrema[3][2];
	ImgColorSpace	destCS = *dstp->GetCS();
	PanImage	inpIm;
	ImgColorSpace	ourCS;
	int		c;
	inpIm.LinkRO(isrc);
	if (PimagesOverlap(dstp->GetImg(), isrc.GetImg())) {
		inpIm = isrc;			// need our own copy
		PcopyCS(&destCS, inpIm.GetCS());
		PrealCS(&destCS);		// assumed destination space
	} else if (destCS.dtype != IDTfloat) {
		DMESG(DMCparameter, "Destination image must be float");
		return false;
	} else if (destCS.format != inpIm.GetCS()->format) {
		DMESG(DMCparameter, "Mismatched image formats");
		return false;
	}
						// get extrema
	memset(b_extrema, 0, sizeof(b_extrema));
	memset(s_extrema, 0, sizeof(s_extrema));
	PcomputeHisto(inp8bit ? (void *)b_extrema : (void *)s_extrema,
					NULL, 0, inpIm.GetImg(), PHCrgb3);
	for (c = inpIm.NComp()*inp8bit; c--; ) {
		s_extrema[c][0] = b_extrema[c][0];
		s_extrema[c][1] = b_extrema[c][1];
	}
	PcopyCS(&ourCS, inpIm.GetCS());		// prepare dequantized holder
	ourCS.dtype = IDTfloat;
	if (!dstp->Init(inpIm.Width(), inpIm.Height(), ourCS))
		return false;
	for (c = inpIm.NComp(); c--; )	// dequantize each component
		if (!dequantizeComponent(dstp, inpIm, c, s_extrema[c]))
			return false;
						// destination color conversion
	return dstp->ConvertCS(destCS, sf);
}

#endif		// end of slower method
