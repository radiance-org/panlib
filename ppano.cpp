/*
 *  ppano.cpp
 *  Photosphere
 *
 *  Class to handle building of panoramas from LDR or HDR input images.
 *
 *  Created by Greg Ward on 5/4/06.
 *  Copyright 2006 Anyhere Software. All rights reserved.
 *
 */

/* #define PblendPano	PblendPan2 */

#include <math.h>
#include "pancine.h"
#include "ppano.h"

#ifndef MAXANCHORS
#define MAXANCHORS	6	/* maximum anchors for single image */
#endif
#ifndef MAXANCHORPCT
#define MAXANCHORPCT	8	/* maximum anchor size (% image size) */
#endif
#ifndef MAXMATCHERR
#define MAXMATCHERR	0.70	/* maximum RMS error for matches */
#endif

#define	scaleWhole(i,sf)	int(float(i)*(sf) + .5f)

static const char	linkName[] = "PanoLink";

// Can we (re)set the given origin for anchor?
bool
PPanorama::CheckOrig(const DBRecord &irec, const ImgRect &anchr) const
{
	if (irec.GetFieldInfo() == NULL)
		return false;
					// check rectangle
	int	xsiz = PDBgetField(irec,PDBFxsize).GetInt();
	int	ysiz = PDBgetField(irec,PDBFysize).GetInt();
	
	if (!PlegalRect(&anchr, xsiz, ysiz))
		return false;
	if ((anchr.xright - anchr.xleft)*100 > xsiz*MAXANCHORPCT)
		return false;
	if ((anchr.ybottom - anchr.ytop)*100 > ysiz*MAXANCHORPCT)
		return false;
					// check if it's allowed
	if (irl.GetSize() <= 0)
		return true;
	int		lri = FindImage(irec);
	if (lri < 0)			// must be in list already
		return false;
	if (irl.Get(lri)[linkID].GetNV() >= MAXANCHORS*3)
		return false;
	return true;
}

// Assign image with link origin (opening)
bool
PPanorama::LinkOrig(const DBRecord &irec, const ImgRect &anchr)
{
	if (!CheckOrig(irec, anchr)) {
		DMESG(DMCinput, "Bad panorama origin anchor");
		return false;
	}
					// find/insert image record
	anchor = InsertImage(irec);
					// assign anchor box
	anchRect = anchr;
	return true;
}

// Verify link destination and check for loops
bool
PPanorama::CheckDest(const DBRecord &irec, const ImgRect &anchr) const
{
	if (anchor < 0)			// is anchor even open?
		return false;
	if (irec.GetFieldInfo() == NULL)
		return false;
					// check densities
	if (PDBgetField(irec,PDBFxdensity) !=
			PDBgetField(irl.Get(anchor),PDBFxdensity) ||
			PDBgetField(irec,PDBFydensity) !=
			    PDBgetField(irl.Get(anchor),PDBFydensity))
		return false;
					// check rectangle
	int	xsiz = PDBgetField(irec,PDBFxsize).GetInt();
	int	ysiz = PDBgetField(irec,PDBFysize).GetInt();
	
	if (!PlegalRect(&anchr, xsiz, ysiz))
		return false;
	if ((anchr.xright - anchr.xleft)*100 > xsiz*MAXANCHORPCT)
		return false;
	if ((anchr.ybottom - anchr.ytop)*100 > ysiz*MAXANCHORPCT)
		return false;
					// check prior links
	int		lri = FindImage(irec);
	if (lri >= 0) {			// replacing existing anchor?
		int32	linka[MAXANCHORS*3];
		int	la = irl.Get(lri)[linkID].Get(linka, MAXANCHORS*3) / 3;
		while (la--)
			if (linka[3*la+2] == anchor)
				break;
		if (la < 0)
			return false;	// not linked to current anchor
	}
	return true;
}

// Assign image with link destination (closing)
bool
PPanorama::LinkDest(const DBRecord &irec, const ImgRect &anchr)
{
	if (!CheckDest(irec, anchr)) {
		DMESG(DMCinput, "Bad panorama destination anchor");
		return false;
	}
	int	matchPoint[4];		// match up anchor regions
	if (!AlignAnchor(matchPoint, irec, anchr))
		return false;
	float	xcorr, ycorr;		// get resolution scaling factors
	if (PDBgetField(irec,PDBFxdensity).Get(&xcorr) &&
			PDBgetField(irec,PDBFydensity).Get(&ycorr)) {
		if (xcorr < ycorr) {
			ycorr = xcorr/ycorr;
			xcorr = 1.f;
		} else {
			xcorr = ycorr/xcorr;
			ycorr = 1.f;
		}
	} else
		xcorr = ycorr = 1.f;
	int32	linka[MAXANCHORS*3];	// assign/replace bidirectional links
	int	i, nv;
	int	lri = InsertImage(irec);
	if (lri < 0)
		return false;
	nv = irl.Get(anchor)[linkID].Get(linka, MAXANCHORS*3);
	for (i = nv/3; i--; )
		if (linka[3*i+2] == lri) break;
	if (i < 0) i = nv;
	else i *= 3;
	linka[i++] = scaleWhole(matchPoint[0], xcorr);
	linka[i++] = scaleWhole(matchPoint[1], ycorr);
	linka[i++] = lri;
	if (i > nv) nv = i;
	irl[anchor].SetField(linkID, linka, nv);
	nv = irl.Get(lri)[linkID].Get(linka, MAXANCHORS*3);
	for (i = nv/3; i--; )
		if (linka[3*i+2] == anchor) break;
	if (i < 0) i = nv;
	else i *= 3;
	linka[i++] = scaleWhole(matchPoint[2], xcorr);
	linka[i++] = scaleWhole(matchPoint[3], ycorr);
	linka[i++] = anchor;
	if (i > nv) nv = i;
	irl[lri].SetField(linkID, linka, nv);
	anchor = -1;			// anchor link is now closed
	return true;
}

// Check to make sure we're ready & get size
const char *
PPanorama::CheckPano(ImgColorSpace *dcsp, int res[2], int org0[2]) const
{
	if (irl.GetSize() < 2)
		return "Too few images for panorama";

	if (dcsp != NULL)		// assign color space
		PcopyCS(dcsp, &cspc);
					// check links & get bounds
	const bool	bound = (res != NULL) | (org0 != NULL);
	int		siz[2];
	ImgRect		rect0;
	if (bound && !GetImageSize(siz, 0))
		return "Bad image size";
	rect0.xleft = rect0.ytop = 0;
	rect0.xright = siz[0];
	rect0.ybottom = siz[1];
	for (int i = irl.GetSize(); --i > 0; ) {
		if (!PmatchOrient(irl.Get(i), irl.Get(0)))
			return "Images do not match orientation";
		int		offset[2];
		ABitMap		bmap(irl.GetSize());
		if (!Linked(0, i, &bmap, bound ? offset : (int *)NULL))
			return "Image set not fully joined";
		if (!bound)
			continue;
		if (!GetImageSize(siz, i))
			return "Bad image size";
		if (offset[0] < rect0.xleft)
			rect0.xleft = offset[0];
		if (offset[0] + siz[0] > rect0.xright)
			rect0.xright = offset[0] + siz[0];
		if (offset[1] < rect0.ytop)
			rect0.ytop = offset[1];
		if (offset[1] + siz[1] > rect0.ybottom)
			rect0.ybottom = offset[1] + siz[1];
	}
	if (res != NULL) {
		res[0] = rect0.xright - rect0.xleft;
		res[1] = rect0.ybottom - rect0.ytop;
	}
	if (org0 != NULL) {
		org0[0] = -rect0.xleft;
		org0[1] = -rect0.ytop;
	}
	return NULL;			// we're ready
}

// Render panorama from current image set
bool
PPanorama::RenderPano(ImgStruct *pimg, ABitMap2 *coverage)
{
	int		org0[2], pres[2];
	const char *	emsg = CheckPano(NULL, pres, org0);
	if (emsg != NULL) {
		DMESG(DMCinput, emsg);
		return false;
	}
	if (reportProgress && !(*reportProgress)("Building Panorama", 0)) {
		DMESG(DMCtrace, "User cancelled Panorama builder");
		return false;
	}
	const bool	freshStart = (pimg->img == NULL);
	if (!PnewImage(pimg, double(pres[1])/double(pres[0])))
		return false;
	if (freshStart) {
		sprintf(dmessage_buf, "Clearing %dx%d canvas",
					pimg->xres, pimg->yres);
		DMESG(DMCtrace, dmessage_buf);
		PclearImage(pimg, NULL);
	}
	if (coverage != NULL)
		coverage->NewBitMap(pimg->xres, pimg->yres);
	const float	isca = float(pimg->xres) / float(pres[0]);
	int		siz[2];
	ImgRect		crect;
	ImgStruct	img0;
	ImgReader *	ir;
					// render starting image
	GetImageSize(siz, 0);
	crect.xleft = scaleWhole(org0[0], isca);
	if (crect.xleft < 0) crect.xleft = 0;
	crect.ytop = scaleWhole(org0[1], isca);
	if (crect.ytop < 0) crect.ytop = 0;
	crect.xright = crect.xleft + scaleWhole(siz[0], isca);
	if (crect.xright > pimg->xres) crect.xright = pimg->xres;
	crect.ybottom = crect.ytop + scaleWhole(siz[1], isca);
	if (crect.ybottom > pimg->yres) crect.ybottom = pimg->yres;
	if (!PlinkSubimage(&img0, pimg, &crect)) {
		DMESG(DMCparameter, "PlinkSubimage failed");
		return false;
	}
	ir = PopenImageD(irl.Get(0), false);
	if (ir == NULL) {
		PfreeImage(pimg);
		return false;
	}
	DMESGF(DMCtrace, "Loading image '%s'", ir->file);
	if (!PrenderImageR(&img0, ir, false)) {
		IRclose(ir);
		PfreeImage(pimg);
		return false;
	}
	IRclose(ir);
	ABitMap		doneMap(irl.GetSize());
					// render linked images
	if (!RenderLinked(pimg, isca, 0, &img0, &doneMap, coverage)) {
		PfreeImage(&img0);
		PfreeImage(pimg);
		return false;
	}
	PfreeImage(&img0);		// unlink first child
	return true;			// success!
}

// Adjust image exposure
static void
adjustExposure(ImgStruct *ip, const float sf)
{
	int	x, y;

	if (ip->csp->dtype != IDTfloat) {
		DMESGF(DMCwarning, "Uncorrected exposure difference (%f)", sf);
		return;
	}
	for (y = 0; y < ip->yres; y++) {
		float *	fp = (float *)ProwPtr(ip, y);
		for (x = ip->xres*ImgPixelLen[ip->csp->format]; x--; )
			*fp++ *= sf;
	}
}

// Render linked children into panorama (private)
bool
PPanorama::RenderLinked(ImgStruct *pimg, const float isca,
		int i1, const ImgStruct *ipar, ABitMap *bm, ABitMap2 *cvg)
{
	int		xo1, yo1;	// get parental origin
	if (!PpixPos(&xo1, &yo1, pimg, ipar->img)) {
		DMESG(DMCassert, "PpixPos failed!");
		return false;
	}
	if (cvg != NULL)		// parent's footprint in mask
		cvg->ClearRect(xo1, yo1, ipar->xres, ipar->yres, true);
	bm->Set(i1);
	if (reportProgress && !(*reportProgress)(NULL,
					100*bm->SumTotal()/irl.GetSize())) {
		DMESG(DMCtrace, "User cancelled panorama builder");
		PfreeImage(pimg);
		return false;
	}
					// read and blend children
	float		stonits = GetStoNits();
	int32		linka[MAXANCHORS*3];
	int		la = irl.Get(i1)[linkID].Get(linka, MAXANCHORS*3) / 3;
	while (la-- > 0) {
		const int	kid = linka[3*la+2];
		if (bm->Check(kid))
			continue;	// already done
		ImgReader *	kidRdr;
		ImgRect		kidRect, blndRect;
		ImgStruct	kidImg, blndImg;
		float		sf;
		int		parAnch[2], kidAnch[2];
		{			// find child's anchor position
			int32		linkb[MAXANCHORS*3];
			int		lb = irl.Get(kid)[linkID].Get(linkb,
							MAXANCHORS*3) / 3;
			while (lb--)
				if (linkb[3*lb+2] == i1)
					break;
			if (lb < 0) {
				DMESG(DMCassert, "Missing link-back");
				return false;
			}
			kidAnch[0] = scaleWhole(linkb[3*lb], isca);
			kidAnch[1] = scaleWhole(linkb[3*lb+1], isca);
		}
		parAnch[0] = scaleWhole(linka[3*la], isca);
		parAnch[1] = scaleWhole(linka[3*la+1], isca);
		kidRdr = PopenImageD(irl.Get(kid), false);
		if (kidRdr == NULL)
			return false;
		DMESGF(DMCtrace, "Loading image '%s'", kidRdr->file);
		kidImg.csp = pimg->csp;
		kidImg.xres = scaleWhole(kidRdr->xres, isca);
		kidImg.yres = scaleWhole(kidRdr->yres, isca);
		kidImg.img = NULL;
		if (!PrenderImageR(&kidImg, kidRdr, false)) {
			IRclose(kidRdr);
			return false;
		}
		IRclose(kidRdr);
		if (stonits > 0 && PDBgetField(irl.Get(kid),
						PDBFstonits).Get(&sf) &&
				((sf /= stonits) < .98f || sf > 1.01f)) {
			DMESGF(DMCtrace, "Adjusting exposure by %f", sf);
			adjustExposure(&kidImg, sf);
		}
		DMESG(DMCtrace, "Computing blend region");
		blndRect.xleft = parAnch[0];
		blndRect.xright = blndRect.xleft - kidAnch[0] + kidImg.xres;
		blndRect.ytop = parAnch[1];
		blndRect.ybottom = blndRect.ytop - kidAnch[1] + kidImg.yres;
		blndImg.csp = pimg->csp;
		blndImg.img = NULL;
		if (!PblendPano(&blndImg, &blndRect, ipar, &kidImg))
			return false;
		DMESG(DMCtrace, "Compositing into panorama");
		kidRect.xleft = xo1 + parAnch[0] - kidAnch[0];
		kidRect.ytop = yo1 + parAnch[1] - kidAnch[1];
		kidRect.xright = kidRect.xleft + kidImg.xres;
		kidRect.ybottom = kidRect.ytop + kidImg.yres;
		if (!PcopyImage(pimg, &kidImg, kidRect.xleft, kidRect.ytop)) {
			PfreeImage(&kidImg);
			PfreeImage(&blndImg);
			return false;
		}
		PfreeImage(&kidImg);
		blndRect.xleft = xo1;
		if (parAnch[0] > kidAnch[0])
			blndRect.xleft += ipar->xres - blndImg.xres;
		blndRect.ytop = yo1;
		if (parAnch[1] > kidAnch[1])
			blndRect.ytop += ipar->yres - blndImg.yres;
		if (!PcopyImage(pimg, &blndImg, blndRect.xleft, blndRect.ytop)) {
			PfreeImage(&blndImg);
			return false;
		}
		PfreeImage(&blndImg);		// now, link to result
		if (kidRect.xleft < 0) kidRect.xleft = 0;
		if (kidRect.xright > pimg->xres) kidRect.xright = pimg->xres;
		if (kidRect.ytop < 0) kidRect.ytop = 0;
		if (kidRect.ybottom > pimg->yres) kidRect.ybottom = pimg->yres;
		if (!PlinkSubimage(&kidImg, pimg, &kidRect)) {
			DMESG(DMCparameter, "PlinkSubimage failed");
			return false;
		}
						// recurse on grandkids
		if (!RenderLinked(pimg, isca, kid, &kidImg, bm, cvg))
			return false;
		PfreeImage(&kidImg);		// done with this child
	}
	return true;
}

// Find a particular record in our image list (private)
int
PPanorama::FindImage(const DBRecord &irec) const
{
	DBQuery		dq;		// find in record list
	if (!PCsetQuery(&dq, irec, true)) {
		DMESG(DMCparameter, "PCsetQuery failed in FindImage");
		return -1;
	}
	return irl.FindFirst(&dq);
}

// Insert image record in our list if not present (private)
int
PPanorama::InsertImage(const DBRecord &irec)
{
	int	res;
	if (irl.GetSize() <= 0) {	// initialize our image record list
		fieldInfo = *irec.GetFieldInfo();
		linkID = fieldInfo.Field(linkName, DBDTint,
				"Panorama Link Rectangles",
				DBFFarray|DBFFhide|DBFFcustom,
				&DBformatInteger);
		if (linkID < 0) {
			DMESG(DMCparameter, "Out of fields in InsertImage");
			return -1;
		}
		irl.Init(&fieldInfo);
	} else if ((res = FindImage(irec)) >= 0)
		return res;
	res = irl.GetSize();		// need to add image
	irl.NewRecord() = irec;
	if ((PDBgetField(irl.Get(res),PDBFxsize).GetInt() <= 0) |
		(PDBgetField(irl.Get(res),PDBFysize).GetInt() <= 0)) {
		DMESG(DMCparameter, "Illegal image size in InsertImage");
		irl.Resize(res);
		return -1;
	}
	if (PDBgetField(irl.Get(res),PDBFcrop).GetNV() > 0) {
		DMESG(DMCwarning, "Ignoring crop for panorama");
		PDBclearField(&irl[res],PDBFcrop);
	}
	return res;
}

// Get rendered dimensions of ith image (private)
bool
PPanorama::GetImageSize(int siz[2], int i) const
{
	if (siz == NULL)
		return false;
	if ((i < 0) | (i >= irl.GetSize()))
		return false;
					// get original image size
	const DBRecord &	irec = irl.Get(i);
	int32			xres, yres;
	if (!PDBgetField(irec,PDBFxsize).Get(&xres) ||
			!PDBgetField(irec,PDBFysize).Get(&yres))
		return false;
					// adjust for pixel aspect ratio
	float			xdens, ydens;
	if (PDBgetField(irec,PDBFxdensity).Get(&xdens) &&
			PDBgetField(irec,PDBFydensity).Get(&ydens)) {
		if (xdens < ydens*.99f)
			yres = scaleWhole(yres, xdens/ydens);
		else if (ydens < xdens*.99f)
			xres = scaleWhole(xres, ydens/xdens);
	}
	siz[0] = xres;
	siz[1] = yres;
	return (siz[0] > 1) & (siz[1] > 1);
}

// Find existing link in image list (private)
bool
PPanorama::Linked(int i1, int i2, ABitMap *bm, int offset[2]) const
{
	if (i1 == i2) {
		if (offset != NULL)
			offset[0] = offset[1] = 0;
		return true;
	}
	if (!bm->TestAndSet(i1))	// already checked this image?
		return false;
	int32	linka[MAXANCHORS*3];
	int	la = irl.Get(i1)[linkID].Get(linka, MAXANCHORS*3) / 3;
	while (la-- > 0)		// search our link list
		if (Linked(linka[3*la+2], i2, bm, offset)) {
			if (offset == NULL)
				return true;
			int32	linkb[MAXANCHORS*3];
			int	lb = irl.Get(linka[3*la+2])[linkID].Get(linkb,
							MAXANCHORS*3) / 3;
			while (lb--)
				if (linkb[3*lb+2] == i1)
					break;
			if (lb < 0) {
				DMESG(DMCassert, "Missing link-back");
				return false;
			}
			offset[0] += linka[3*la] - linkb[3*lb];
			offset[1] += linka[3*la+1] - linkb[3*lb+1];
			return true;
		}
	return false;
}

// Compute normalized RMS difference between two images
static double
imageDiff(const ImgStruct *img1, const ImgStruct *img2)
{
	DASSERT((img1->xres == img2->xres) & (img1->yres == img2->yres));
	DASSERT((img1->csp->dtype == IDTfloat) & (img2->csp == img1->csp));
	double		diff2 = .0;
	int		nv = 0;
	for (int y = img1->yres; y--; ) {
		const float *	fp1 = (const float *)ProwPtr(img1,y);
		const float *	fp2 = (const float *)ProwPtr(img2,y);
		for (int i = img1->xres*ImgPixelLen[img1->csp->format]; i--; ) {
			double	d = double(*fp1 + *fp2);
			if (d <= .0) { fp1++; fp2++; continue; }
			d = double(*fp1++ - *fp2++)/d;
			diff2 += d*d;
			++nv;
		}
	}
	if (nv < 1) return MAXMATCHERR;
	return 2.*sqrt( diff2 / double(nv) );
}

// Compute dot product of two images
static double
imageProd(const ImgStruct *img1, const ImgStruct *img2)
{
	DASSERT((img1->xres == img2->xres) & (img1->yres == img2->yres));
	DASSERT((img1->csp->dtype == IDTfloat) & (img2->csp == img1->csp));
	double		prod = .0;
	for (int y = img1->yres; y--; ) {
		const float *	fp1 = (const float *)ProwPtr(img1,y);
		const float *	fp2 = (const float *)ProwPtr(img2,y);
		for (int i = img1->xres*ImgPixelLen[img1->csp->format]; i--; )
			prod += double(*fp1++)*double(*fp2++);
	}
	return prod;
}

// Find prominent (corner) point in image
static void
findPromonotory(int pt[2], const ImgStruct *img)
{
	const int		krad = 10;
	static ImgStruct	kern;
	double			maxResp = .0;
	int			xo, yo;
	
	if (kern.img == NULL) {		// initialize corner kernel (saddle)
		kern.csp = &ICS_Y;
		kern.xres = 2*krad + 1;
		kern.yres = 2*krad + 1;
		if (!PsetImage(&kern, Pblack)) return;
		for (yo = -krad; yo <= krad; yo++) {
			const int	xr = int(sqrt(double(krad*krad - yo*yo)) + .5);
			float *		kp = (float *)PpixPtr(&kern,krad-xr,krad+yo);
			for (xo = -xr; xo <= xr; xo++)
				*kp++ = sqrt(double(xo*xo + yo*yo)) *
					sin(2.*atan2(double(yo),double(xo)));
		}
	}
	DASSERT(img->csp == kern.csp);
					// find maximum corner response
	pt[0] = img->xres/2;
	pt[1] = img->yres/2;
	for (yo = -(img->yres - kern.yres); yo <= img->yres - kern.yres; yo++)
	    for (xo = -(img->xres - kern.xres); xo <= img->xres - kern.xres; xo++) {
		ImgRect		rct;
		ImgStruct	isub;
		double		resp;
		rct.xleft = (img->xres - kern.xres + xo)/2;
		rct.xright = rct.xleft + kern.xres;
		rct.ytop = (img->yres - kern.yres + yo)/2;
		rct.ybottom = rct.ytop + kern.yres;
		if (!PlinkSubimage(&isub, img, &rct)) {
			DMESG(DMCparameter, "PlinkSubimage failed!");
			return;
		}
		resp = imageProd(&isub, &kern);
		PfreeImage(&isub);
		resp *= resp;
		if (resp <= maxResp)
			continue;
		pt[0] = (img->xres + xo)/2;
		pt[1] = (img->yres + yo)/2;
		maxResp = resp;
	    }
}

// Find best matching position for image pair
static bool
matchPoint(int match[4], const ImgStruct *img1, const ImgStruct *img2)
{
	DASSERT((img1->xres == img2->xres) & (img1->yres == img2->yres));
	const int	xsdia = img1->xres/2;
	const int	ysdia = img1->yres/2;
	const int	xsubres = img1->xres - xsdia;
	const int	ysubres = img1->yres - ysdia;
	double		minErr = 100.*MAXMATCHERR;
	int		xoffs, yoffs;
	ImgRect		rect1, rect2;
	ImgStruct	isub1, isub2;
	double		err;
					// symmetric search for best offset
	for (yoffs = -ysdia; yoffs < ysdia; yoffs++) {
		rect1.ytop = (img1->yres - ysubres + yoffs + 1)/2;
		rect1.ybottom = rect1.ytop + ysubres;
		rect2.ytop = (img2->yres - ysubres - yoffs)/2;
		rect2.ybottom = rect2.ytop + ysubres;
		for (xoffs = -xsdia; xoffs < xsdia; xoffs++) {
			rect1.xleft = (img1->xres - xsubres + xoffs + 1)/2;
			rect1.xright = rect1.xleft + xsubres;
			rect2.xleft = (img2->xres - xsubres - xoffs)/2;
			rect2.xright = rect2.xleft + xsubres;
			if (!PlinkSubimage(&isub1, img1, &rect1) ||
					!PlinkSubimage(&isub2, img2, &rect2)) {
				DMESG(DMCparameter, "PlinkSubimage failed!");
				return false;
			}
			err = imageDiff(&isub1, &isub2);
			PfreeImage(&isub1); PfreeImage(&isub2);
			if (err >= minErr)
				continue;
			match[0] = xoffs;
			match[1] = yoffs;
			minErr = err;
		}
	}
	if (minErr > MAXMATCHERR) {
		DMESGF(DMCdata, "Cannot align image features (%.0f%% RMS error)",
				100.*minErr);
		return false;
	}
	xoffs = match[0];		// find promonotory
	yoffs = match[1];
	rect1.xleft = (img1->xres - xsubres + xoffs + 1)/2;
	rect1.xright = rect1.xleft + xsubres;
	rect1.ytop = (img1->yres - ysubres + yoffs + 1)/2;
	rect1.ybottom = rect1.ytop + ysubres;
	PlinkSubimage(&isub1, img1, &rect1);
	findPromonotory(match, &isub1);
	PfreeImage(&isub1);
	match[2] = match[0] + (img2->xres - xsubres - xoffs)/2;
	match[3] = match[1] + (img2->yres - ysubres - yoffs)/2;
	match[0] += rect1.xleft;
	match[1] += rect1.ytop;
	return true;
}

// Align our anchor link and find promonotory (private)
bool
PPanorama::AlignAnchor(int matchpt[4], const DBRecord &irec2,
			const ImgRect &anchr2) const
{
	if (anchor < 0)
		return false;
	const DBRecord &	irec1 = irl.Get(anchor);
	bool			do_reporting = (reportProgress != NULL);
					// expand search areas
	ImgRect			rect1, rect2;
	int			siz, res;
	siz = PrectArea(&anchRect);
	res = PrectArea(&anchr2);
	if (res > siz)
		siz = res;
	siz = int(sqrt(double(siz) + .5));
	siz += siz/2 + 16;
	for ( ; ; ) {			// make square areas fit within images
		if ((rect1.xleft = (anchRect.xleft + anchRect.xright - siz)/2) < 0) {
			siz = anchRect.xleft + anchRect.xright;
			continue;
		}
		res = PDBgetField(irec1,PDBFxsize).GetInt();
		if ((rect1.xright = rect1.xleft + siz) > res) {
			siz = 2*res - 1 - anchRect.xleft - anchRect.xright;
			continue;
		}
		if ((rect1.ytop = (anchRect.ytop + anchRect.ybottom - siz)/2) < 0) {
			siz = anchRect.ytop + anchRect.ybottom;
			continue;
		}
		res = PDBgetField(irec1,PDBFysize).GetInt();
		if (res < 1150)
			do_reporting = false;
		if ((rect1.ybottom = rect1.ytop + siz) > res) {
			siz = 2*res - 1 - anchRect.ytop - anchRect.ybottom;
			continue;
		}
		if ((rect2.xleft = (anchr2.xleft + anchr2.xright - siz)/2) < 0) {
			siz = anchr2.xleft + anchr2.xright;
			continue;
		}
		res = PDBgetField(irec2,PDBFxsize).GetInt();
		if ((rect2.xright = rect2.xleft + siz) > res) {
			siz = 2*res - 1 - anchr2.xleft - anchr2.xright;
			continue;
		}
		if ((rect2.ytop = (anchr2.ytop + anchr2.ybottom - siz)/2) < 0) {
			siz = anchr2.ytop + anchr2.ybottom;
			continue;
		}
		res = PDBgetField(irec2,PDBFysize).GetInt();
		if (res < 1150)
			do_reporting = false;
		if ((rect2.ybottom = rect2.ytop + siz) > res) {
			siz = 2*res - 1 - anchr2.ytop - anchr2.ybottom;
			continue;
		}
		break;
	}
					// load associated image regions
	ImgStruct		img1, img2, itmp;
	ImgReader *		ir;
	ImgReadBuf		rbuf;
	ir = PopenImageD(irec1, false);
	if (ir == NULL)
		return false;
	PcopyCS(&cspc, &ir->cs);
	img1.csp = &ir->cs;
	img1.xres = rect1.xright - rect1.xleft;
	img1.yres = rect1.ybottom - rect1.ytop;
	img1.img = NULL;
	if (!PnewImage(&img1, .0))
		return false;
	rbuf.cs = ir->cs;
	rbuf.r = rect1;
	rbuf.subsample = 1;
	rbuf.buf = img1.img;
	sprintf(dmessage_buf, "Reading feature from '%s'",
				PgetFilename(ir->file));
	if (do_reporting && !(*reportProgress)(dmessage_buf, 0))
		return false;
	DMESG(DMCtrace, dmessage_buf);
	if (IRreadRec(ir, &rbuf) != IREnone) {
		PreportReaderErr(ir);
		IRclose(ir);
		PfreeImage(&img1);
		return false;
	}
	IRclose(ir);
	if (rbuf.buf != img1.img) {
		DMESG(DMCparameter, "Read buffer problem in AlignAnchor");
		free(rbuf.buf);
		PfreeImage(&img1);
		return false;
	}
	if (!PconvertColorSpace(&img1, &ICS_Y, 1.f))
		return false;
	ir = PopenImageD(irec2, false);
	if (ir == NULL)
		return false;
	if (!PmatchColorSpace(&ir->cs, &cspc, PICMall))
		DMESG(DMCwarning, "Mismatched image color spaces");
	img2.csp = &ir->cs;
	img2.xres = rect2.xright - rect2.xleft;
	img2.yres = rect2.ybottom - rect2.ytop;
	img2.img = NULL;
	if (!PnewImage(&img2, .0))
		return false;
	rbuf.cs = ir->cs;
	rbuf.r = rect2;
	rbuf.subsample = 1;
	rbuf.buf = img2.img;
	sprintf(dmessage_buf, "Reading feature from '%s'",
				PgetFilename(ir->file));
	if (do_reporting && !(*reportProgress)(dmessage_buf, 30))
		return false;
	DMESG(DMCtrace, dmessage_buf);
	if (IRreadRec(ir, &rbuf) != IREnone) {
		PreportReaderErr(ir);
		IRclose(ir);
		PfreeImage(&img2);
		PfreeImage(&img1);
		return false;
	}
	IRclose(ir);
	if (rbuf.buf != img2.img) {
		DMESG(DMCparameter, "Read buffer problem in AlignAnchor");
		free(rbuf.buf);
		PfreeImage(&img2);
		PfreeImage(&img1);
		return false;
	}
	float	stonits1, stonits2;
	if (!PDBgetField(irec1,PDBFstonits).Get(&stonits1) ||
			!PDBgetField(irec2,PDBFstonits).Get(&stonits2))
		stonits1 = stonits2 = 1.f;
	if (!PconvertColorSpace(&img2, &ICS_Y, stonits2/stonits1))
		return false;
	if (do_reporting && !(*reportProgress)("Matching features", 60))
		return false;
	DMESG(DMCtrace, "Attempting to match image features");
	itmp.csp = img1.csp;			// blur to reduce noise
	itmp.xres = img1.xres;
	itmp.yres = img1.yres;
	itmp.img = NULL;
	if (!PblurImage(&itmp, &img1, 1.7f)) {
		PfreeImage(&img1);
		PfreeImage(&img2);
		return false;
	}
	PfreeImage(&img1);
	img1 = itmp;
	itmp.csp = img2.csp;
	itmp.xres = img2.xres;
	itmp.yres = img2.yres;
	itmp.img = NULL;
	if (!PblurImage(&itmp, &img2, 1.7f)) {
		PfreeImage(&img1);
		PfreeImage(&img2);
		return false;
	}
	PfreeImage(&img2);
	img2 = itmp;
					// find matching position
	if (!matchPoint(matchpt, &img1, &img2)) {
		PfreeImage(&img2);
		PfreeImage(&img1);
		return false;
	}
	PfreeImage(&img2);
	PfreeImage(&img1);
	matchpt[0] += rect1.xleft;
	matchpt[1] += rect1.ytop;
	matchpt[2] += rect2.xleft;
	matchpt[3] += rect2.ytop;
	if (do_reporting)
		(*reportProgress)("Done", 100);
	return true;
}
