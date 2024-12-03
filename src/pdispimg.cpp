/*
 *  pdispimg.cpp
 *  panlib
 *
 *  Class implementation to facilitate HDR image dispaly.
 *
 *  Created by gward on Tue Jan 08 2002.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#include "pancine.h"
#include "pdispimg.h"
#include "resolu.h"
#include "tiffio.h"
#include "tmaptiff.h"
#include <math.h>

int		PEMDefault = PEMauto;		// default exposure mode

// Enable/change false color mapping for HDR image
void
PImageRing::FalseColorOn(double maxv, double minv)
{
	if (!IsHDR())
		return;
	bool	fcChanged = false;
	if (fcs == NULL) {
		fcs = fcInit(NULL);
		if (fcs == NULL)
			return;
		fcChanged = true;
	}
	if (!(fcManual = (maxv > 1e-7))) {	// auto scaling?
		changed = true;			// force update
		return;
	}
	float	sfcorr;				// apply reverse correction
	if (GetRec()[PDBFstonits].Get(&sfcorr)) {
		sfcorr = GetSF() / sfcorr;
		minv *= sfcorr;
		maxv *= sfcorr;
	}
	fcChanged |= (fcs->mbrmax != tmCvLuminance(maxv));
	if ((1e-7 < minv) & (minv < maxv)) {	// log mapping
		fcChanged |= (fcs->mbrmin != tmCvLuminance(minv));
		if (fcChanged)
			fcFixedLog(fcs, minv, maxv);
	} else {				// linear mapping
		fcChanged |= (fcIsLogMap(fcs) != 0);
		if (fcChanged)
			fcFixedLinear(fcs, maxv);
	}
	changed |= fcChanged;			// update changed flag
}

// Get false color limits, return true if manual mode
bool
PImageRing::GetFCLimits(double minmax[2]) const
{
	if (fcs == NULL || fcs->lumap == NULL) {
		minmax[0] = minmax[1] = .0;
		return false;
	}
	minmax[1] = tmLuminance(fcs->mbrmax);
	if (fcIsLogMap(fcs) > 0)
		minmax[0] = tmLuminance(fcs->mbrmin);
	else
		minmax[0] = .0;
	float	sfcorr;				// apply correction
	if (GetRec()[PDBFstonits].Get(&sfcorr)) {
		sfcorr /= GetSF();
		minmax[0] *= sfcorr;
		minmax[1] *= sfcorr;
	}
	return (fcManual > 0);
}

// Check our current settings and update
void
PImageRing::CheckSettings()
{
	// used to return if !changed, but this caused problems for fixed exp.
	int	tmflags = 0;
	if (!(expMode & PEMauto))
		tmflags |= TM_F_LINEAR;
	if (expMode & PEMhuman)
		tmflags |= (TM_F_HCONTR|TM_F_MESOPIC);
	SetTMFlags(tmflags);
	if ((froffs != 0) | (frsmode == IRSabs)) {
		if (!PDisplayImage::SeekFrame(froffs,frsmode) & (frsmode==IRSrel))
			if (froffs < 0)
				PDisplayImage::SeekFrame(0);
			else
				PDisplayImage::SeekFrame(GetNFrames()-1);
	}
	froffs = 0; frsmode = IRSrel;
	if (!IsHDR()) {			// calls CheckLoad (important!)
		expMode = 0;
		FalseColorDone();
	} else if (PEMisManual(expMode) && fcManual < 0)
		SetFixedLinear(expMult);
	changed = false;
}

// Render the current image into the given window (original orientation)
bool		
PImageRing::RenderImage(ImgStruct *ims, const void *pf, const ImgRect *sir)
{
	if (!Ready())
		return false;
	if (ims == NULL || (ims->xres <= 0) | (ims->yres <= 0))
		return false;
	CheckSettings();			// synchronize settings
	bool	hasbord = false;		// compute subimage source/dest.
	ImgRect	rrect, srect;
	rWidth = ims->xres;
	rHeight = ims->yres;
	if (!IsZoomToFit()) {			// fixed zoom or zoom-to-fit?
		curZNum = zoomNum;
		curZDenom = zoomDenom;
	} else if (ims->xres*ior.yres > ims->yres*ior.xres) {
		curZNum = ims->yres;
		curZDenom = ior.yres;
	} else {
		curZNum = ims->xres;
		curZDenom = ior.xres;
	}
	int		ileft = 0;		// get our image rectangle
	int		itop = 0;
	ImgStruct	sidest;
	bool		ok;
	if (sir != NULL) {			// write to subsection of image?
		SetPurgeable(false);		// avoid purge
		ok = PnewImage(ims, .0);	// make sure it's allocated
		SetPurgeable(true);
		if (!ok)
			return false;
		if (!PlinkSubimage(&sidest, ims, sir)) {
			DMESG(DMCparameter, "Bad subrectangle in RenderImage");
			return false;
		}
		if (sidest.mbase != NULL)	// don't retain subimage
			MobjRelease(sidest.mbase);
		ileft = sir->xleft;
		itop = sir->ytop;
		ims = &sidest;			// render into subimage
	}
	MapToSource(&rrect.xleft, &rrect.ytop, ileft, itop);
	rrect.xright = rrect.xleft + ims->xres*curZDenom/curZNum;
	rrect.ybottom = rrect.ytop + ims->yres*curZDenom/curZNum;
	rrect.xright += (rrect.xright == rrect.xleft);
	rrect.ybottom += (rrect.ybottom == rrect.ytop);
						// adjust for borders
	srect.xleft = srect.ytop = 0;
	srect.xright = ims->xres; srect.ybottom = ims->yres;
	if (rrect.xleft < 0) {
		srect.xleft = -rrect.xleft * curZNum / curZDenom;
		rrect.xleft = 0;
		hasbord = true;
	}
	if (rrect.xright > ior.xres) {
		srect.xright -= (rrect.xright - ior.xres) * curZNum / curZDenom;
		rrect.xright = ior.xres;
		hasbord = true;
	}
	if (rrect.ytop < 0) {
		srect.ytop = -rrect.ytop * curZNum / curZDenom;
		rrect.ytop = 0;
		hasbord = true;
	}
	if (rrect.ybottom > ior.yres) {
		srect.ybottom -= (rrect.ybottom - ior.yres) * curZNum / curZDenom;
		rrect.ybottom = ior.yres;
		hasbord = true;
	}
	ImgStruct *	rimg = ims;
	ImgStruct	subimg;
						// destination outside source?
	if (hasbord) {
		SetPurgeable(false);		// avoid purge
		ok = PnewImage(ims, .0);
		SetPurgeable(true);
		if (!ok)
			return false;
		PclearImage(ims, pf);		// create background
		if (!PlinkSubimage(&subimg, ims, &srect)) {
			DMESG(DMCparameter, "Border botch in RenderImage");
			return false;
		}
		rimg = &subimg;			// render into subimage
	}
						// render (sub)image
	if (expMode & PEMlocal)			// use local exposure?
		ok = GetSubimage(rimg, &rrect, &expMult);
	else					// else global exposure
		ok = GetSubimage(rimg, &rrect);
	if (rimg != ims)			// free rendering subimage
		PfreeImage(rimg);
	if (!ok)				// error return
		return false;
	if (!(expMode & PEMlocal)) {		// compute global exposure
		double	em = GetMultiplier();
		if (expMult <= 0. || fabs(em/expMult - 1.) > .05)
			expMult = em;		// don't accrue errors
	}
	return true;
}

// Render the current image into the given window (final orientation)
bool		
PImageRing::MapImage(ImgStruct *ims, const void *pf, const ImgRect *sir)
{
	if (!Ready())
		return false;
	if (ior.IsNoop())			// no need to reorient
		return RenderImage(ims, pf, sir);
	if (ims == NULL || (ims->xres <= 0) | (ims->yres <= 0))
		return false;
	POrient	rOr = ior;			// get size & orientation
	if (XYswapped()) {
		rWidth = ims->yres;
		rHeight = ims->xres;
	} else {
		rWidth = ims->xres;
		rHeight = ims->yres;
	}
	rOr.SetOrigSize(rWidth, rHeight);
	ImgStruct	sidest, orgimg;
	bool		ok;
	SetPurgeable(false);			// avoid purge
	ok = PnewImage(ims, .0);		// make sure we're allocated
	SetPurgeable(true);
	if (!ok)
		return false;
	orgimg.csp = ims->csp;			// render (sub)image
	orgimg.xres = rWidth;
	orgimg.yres = rHeight;
	orgimg.img = NULL;
	if (sir != NULL) {			// rendering subsection?
		int		x0, y0, x1, y1;
		ImgRect		orgrect;
		if (!PlinkSubimage(&sidest, ims, sir)) {
			DMESG(DMCparameter, "Bad subrectangle in MapImage");
			return false;
		}
		if (sidest.mbase != NULL)	// don't retain subimage
			MobjRelease(sidest.mbase);
		ims = &sidest;
		rOr.MapToOrig(&x0, &y0, sir->xleft, sir->ytop);
		rOr.MapToOrig(&x1, &y1, sir->xright-1, sir->ybottom-1);
		if (x0 <= x1) {
			orgrect.xleft = x0;
			orgrect.xright = x1 + 1;
		} else {
			orgrect.xleft = x1;
			orgrect.xright = x0 + 1;
		}
		if (y0 <= y1) {
			orgrect.ytop = y0;
			orgrect.ybottom = y1 + 1;
		} else {
			orgrect.ytop = y1;
			orgrect.ybottom = y0 + 1;
		}
		// XXX very wasteful to allocate whole image for subsection
		ok = RenderImage(&orgimg, pf, &orgrect) &&
				PlinkSubimage(&orgimg, &orgimg, &orgrect);
		rOr.SetOrigSize(orgimg.xres, orgimg.yres);
	} else {
		ok =  RenderImage(&orgimg, pf);
	}
	if (ok)					// reorient result
		ok = rOr.OrientImage(ims, &orgimg);
	PfreeImage(&orgimg);			// free temp. image
	return ok;
}

// Get the currect selection in rendered (window pane) coordinates
bool
PImageRing::GetSelection(ImgRect *rp) const
{
	if ((this == NULL) | (rp == NULL))
		return false;
	if (!PlegalRect(&srcrect, ior.xres, ior.yres))
		return false;
	int	x0, y0, x1, y1;			// map to rendered rectangle
	POrient	rOr = ior;
	rOr.SetOrigSize(rWidth, rHeight);
	MapToDest(&x0, &y0, srcrect.xleft, srcrect.ytop);
	MapToDest(&x1, &y1, srcrect.xright, srcrect.ybottom);
	rOr.MapToCorrect(&x0, &y0, x0, y0);
	rOr.MapToCorrect(&x1, &y1, x1, y1);
	if (x0 < x1) {
		rp->xleft = x0;
		rp->xright = x1;
	} else if (x0 > x1) {
		rp->xleft = x1;
		rp->xright = x0;
	} else
		return false;
	if (y0 < y1) {
		rp->ytop = y0;
		rp->ybottom = y1;
	} else if (y0 > y1) {
		rp->ytop = y1;
		rp->ybottom = y0;
	} else
		return false;
	return true;
}

// Set selection rectangle based on rendered point pair
void
PImageRing::SelectRect(int x0, int y0, int x1, int y1)
{
	POrient	rOr = ior;
						// compute source rectangle
	rOr.SetOrigSize(rWidth, rHeight);
	rOr.MapToOrig(&x0, &y0, x0, y0);
	rOr.MapToOrig(&x1, &y1, x1, y1);
	MapToSource(&x0, &y0, x0, y0);
	MapToSource(&x1, &y1, x1, y1);
	if (x0 <= x1) {
		srcrect.xleft = x0;
		srcrect.xright = x0<x1 ? x1 : x0+1;
	} else {
		srcrect.xleft = x1;
		srcrect.xright = x0;
	}
	if (y0 <= y1) {
		srcrect.ytop = y0;
		srcrect.ybottom = y0<y1 ? y1 : y0+1;
	} else {
		srcrect.ytop = y1;
		srcrect.ybottom = y0;
	}
						// make it legal
	if (srcrect.xleft < 0)
		srcrect.xleft = 0;
	if (srcrect.xright > ior.xres)
		srcrect.xright = ior.xres;
	if (srcrect.ytop < 0)
		srcrect.ytop = 0;
	if (srcrect.ybottom > ior.yres)
		srcrect.ybottom = ior.yres;
}

// Set selection based on original image rectangle
void
PImageRing::SetOrigSelection(ImgRect *orp)
{
	if (orp == NULL)
		return;
						// convert from original coords
	GetViewRect(&srcrect, orp);
						// make it legal
	if (srcrect.xleft < 0)
		srcrect.xleft = 0;
	if (srcrect.xright > ior.xres)
		srcrect.xright = ior.xres;
	if (srcrect.ytop < 0)
		srcrect.ytop = 0;
	if (srcrect.ybottom > ior.yres)
		srcrect.ybottom = ior.yres;
}

// Get the (approximate) multiplier corresponding to the given tone-mapping
double
PDisplayList::GetMultiplier(TMstruct *tm)
{
	if (tm == NULL || tm->lumap == NULL)
		return 1.;
	int	li = tm->mbrmax;		// compute at maximum
	while (li > tm->mbrmin && tm->lumap[li-tm->mbrmin] >= TM_BRES)
		--li;
	return pow((tm->lumap[li-tm->mbrmin]+.5)*(1./TM_BRES),tm->mongam) *
			tm->inpsf / tmLuminance(li);
}

// Store a new image record along with its field information
void
PDisplayImage::SetRecord(const DBRecord &pr)
{
	if (!pr.GetNAlloc())
		return;
	if (&pr == &ircd)
		return;
	if (pr.GetFieldInfo() == &finf) {	// avoid munging source field info
		ircd.Init(&finf, true);
		ircd = pr;
		return;
	}
	int	idMap[DB_MAXFIELD];
	if (!PDBFInfo.MapIDs(idMap, pr.GetFieldInfo())) {
		ircd.Init(&PDBFInfo, true);
		ircd = pr;
		return;
	}
	finf = PDBFInfo;			// add new fields at end
	for (int i = 0; i < pr.GetNFields(); i++)
		if (idMap[i] < 0)
			idMap[i] = finf.Field(pr.GetFieldName(i), pr.GetFieldDetails(i));
	ircd.Init(&finf, true);
	ircd = pr;
}

// Set our display object to use the given image file
bool
PDisplayImage::SetImage(const char *fpath)
{
	CloseImage();				// close previous image
	if (fpath == NULL)
		return false;
	ircd.Init(&PDBFInfo, true);
	PgrokFile(&ircd, fpath);
	return SetImage(ircd);
}

// Set our display object to use the given image record
bool
PDisplayImage::SetImage(const DBRecord &pr)
{
	if (&ircd != &pr) {
		CloseImage();			// close previous image
		SetRecord(pr);
	}
	if (ircd.GetField(PDBFfile) == NULL ||
			ircd.GetField(PDBFlocation) == NULL ||
			!ircd[PDBFxsize].Get(&oxres) ||
			!ircd[PDBFysize].Get(&oyres)) {
		DMESG(DMCdata, "Cannot parse image file");
		ircd.Init();
		oxres = oyres = 0;
		return false;
	}
	need2reload.ClearBitMap();
	UpdateRec(ircd);
	return true;
}

// Open the given frame number
bool
PDisplayImage::SetFrame(int fn)
{
	if (!Ready())
		return false;
	if (fn < 0)
		return false;
	if ((nframes > 0) & (fn >= nframes))
		return false;
	if (dlhead == NULL)			// need to start list?
		dlhead = dlp = new PDisplayList(this);
	else if (fn < dlp->frame)
		dlp = dlhead;
	if (fn == 0)				// just looking for frame 0?
		return dlp->CheckLoad();
	if (ir == NULL) {			// open image
		ir = PopenImageD(ircd, quiet);
		if (ir == NULL)
			return false;
		if ((ir->nframes > 0) & (fn >= ir->nframes))
			goto fail;
	}
						// precheck frame exists
	if (nframes == 0 && IRseekFrame(ir, fn, IRSabs) != IREnone)
		goto fail;
	while (fn > dlp->frame) {		// skip to frame
		if (dlp->next == NULL)		// append new frame to list
			dlp->next = new PDisplayList(this, dlp->frame+1);
		dlp = dlp->next;
	}
	return dlp->CheckLoad();
fail:
	if (ir != NULL) {
		nframes = ir->nframes;
		IRclose(ir); ir = NULL;
	}
	return false;
}

// Seek to the given frame in a sequence
bool
PDisplayImage::SeekFrame(int offs, ImgSeekMode sm)
{
	if (sm == IRSabs)
		return SetFrame(offs);
	if (!CheckLoad())
		return false;
	if (offs == 0)
		return true;
	if (nframes == 1)			// nowhere to go
		return false;
	DASSERT(dlp != NULL);
	int	goal = dlp->frame + offs;
	if (nframes == 0) {			// need to search...
		if ((sm == IRSadv) & (goal < 0))
			if (frameType == IRFbackforth)
				goal = -goal;
			else if (frameType == IRFloop)
				goal = MAXFRAMES;
		if (SetFrame(goal))
			return true;
		if (sm == IRSrel)
			return false;
		if (nframes == 0) {
			DMESG(DMCparameter, "Cannot determine number of frames");
			return false;
		}
	}
	goal = dlp->frame;			// convert to absolute seek
	if (PabsFrame(&goal, nframes, frameType, offs, sm))
		return SetFrame(goal);
	return false;
}

// Update our image record and recompute orientation/size
bool
PDisplayImage::UpdateRec(const DBRecord &pr)
{
	if (!pr.GetNAlloc())
		return false;
	bool	changed = ( &ircd != &pr &&
			( !PmatchOrient(ircd, pr) ||
			ircd[PDBFcrop] != PDBgetField(pr,PDBFcrop) ||
			ircd[PDBFspotexp] != PDBgetField(pr,PDBFspotexp) ) );
	PDisplayList *	dli;
	if (&ircd != &pr &&
		(ircd[PDBFredeye].GetNV() > PDBgetField(pr,PDBFredeye).GetNV() ||
			((dlp == NULL || dlp->tms == NULL) &&
			  ircd[PDBFspotexp] != PDBgetField(pr,PDBFspotexp)) ||
			ircd[PDBFhashval] != PDBgetField(pr,PDBFhashval) ||
			ircd[PDBFnbytes] != PDBgetField(pr,PDBFnbytes))) {
		Reload();			// red-eye or image changed on us
		PDBgetField(pr,PDBFxsize).Get(&oxres);
		PDBgetField(pr,PDBFysize).Get(&oyres);
		for (dli = dlhead; dli != NULL; dli = dli->next) {
			if (dli->tms != NULL) {
				tmDone(dli->tms);
				dli->tms = NULL;
			}
			dli->oxres = oxres;
			dli->oyres = oyres;
		}
		changed = true;
	}
	if (&ircd != &pr &&			// remove new red-eye(s)
		ircd[PDBFredeye].GetNV() < PDBgetField(pr,PDBFredeye).GetNV() &&
			dlp->GetPointers()) {
		dlp->RemoveRedEye(&pr);
		dlp->FreePointers();
		changed = true;
	}
	SetRecord(pr);
	ImgRect	cr;				// reassign ior regardless
	float	xdens, ydens;
	if (ircd[PDBFxdensity].Get(&xdens) && ircd[PDBFydensity].Get(&ydens))
		pixAspect = ydens/xdens;
	else
		pixAspect = 1.f;
	CropRect(&cr);
	int	xr = cr.xright - cr.xleft;
	int	yr = cr.ybottom - cr.ytop;
	CorrectAspect(&xr, &yr);
	ior.SetRecord(&ircd);
	ior.SetOrigSize(xr, yr);
	if (changed)				// need new histograms
		for (dli = dlhead; dli != NULL; dli = dli->next)
			if (dli->tms != NULL)
				tmClearHisto(dli->tms);
	return changed;
}

// Close the current image if one is open
bool
PDisplayImage::CloseImage(bool cancancel)
{
	if (!Ready())
		return true;
	if (cloz != NULL && !(*cloz)(this, cancancel))
		return false;		// callback cancelled operation
	cloz = NULL;
	IRclose(ir); ir = NULL;
	if (dlhead != NULL) {
		dlhead->HolderRelease();
		dlhead = dlp = NULL;
	}
	FalseColorDone();
	ior.SetRecord(NULL);
	ircd.Init();
	return true;
}

// Callback to delete temporary image on closing -- never cancels
bool
PDeleteTempImage(PDisplayImage *di, bool)
{
	char	ipath[1024];
	if (PcurrentFile(ipath, sizeof(ipath), di->GetRec()) == NULL)
		return true;
	if (remove(ipath) < 0)
		DMESGF(DMCwarning, "Cannot delete temp image '%s'", ipath);
	return true;
}

// Check if our image needs to be (re)loaded
bool
PDisplayList::CheckLoad()
{
	if (parent == NULL)
		return false;
						// check if we need to reload
	if (parent->need2reload.Check(frame) && !InUse()) {
		Free();
		parent->need2reload.Reset(frame);
	}
	if (objSize <= 0) {			// (re)load image
		ReleaseObject(GetCacheObject());
		if (objSize <= 0)
			return false;
	}
	parent->UpdateRec(parent->ircd);	// sync image dimensions
	return true;
}

// Get subimage region from cache image
bool
PDisplayList::GetSubimage(ImgStruct *ims, const ImgRect *r, double *localexp)
{
	if ((ims == NULL) | (parent == NULL))
		return false;
	if (!PmatchColorSpace(ims->csp, c_img.csp, PICMptype)) {
		if (!parent->quiet)
			DMESG(DMCparameter, "Request for illegal pixel type");
		return false;
	}
	if (!parent->quiet && !PmatchColorSpace(ims->csp, c_img.csp, PICMall))
		DMESG(DMCwarning, "Non-matching color space request");

	if (!GetPointers())			// check/load/convert image
		return false;
						// set up rendering region
	bool	dolocal = (localexp != NULL);
	bool	ok;
	ImgRect	rdflt;
	if (r == NULL) {
		rdflt.xleft = 0; rdflt.xright = parent->ior.xres;
		rdflt.ytop = 0; rdflt.ybottom = parent->ior.yres;
		r = &rdflt;
		dolocal = false;
	}
	if (tms != NULL) {			// HDR image?
		if (dolocal)
			ok = MapSubimage(ims, r, localexp);
		else {
			ok = MapSubimage(ims, r);
			if (localexp != NULL)
				*localexp = GetMultiplier();
		}
		FreePointers();
		return ok;
	}
	if ((ims->mbase != NULL) & (ims->mbase == c_img.mbase))
		PfreeImage(ims);		// don't overwrite ourselves

	if (ims->img == NULL &&			// just use reference?
			ims->xres == r->xright - r->xleft &&
			ims->yres == r->ybottom - r->ytop) {
		ok = PlinkSubimage(ims, &c_img, r);
	} else {
		ImgStruct	isub;		// else render from subimage
		isub.img = NULL;
		ok = PlinkSubimage(&isub, &c_img, r) &&
				PsizeImage(ims, &isub, PSbest);
		PfreeImage(&isub);
	}
	FreePointers();				// done with cache
	if (localexp != NULL)
		*localexp = 1.;
	return ok;
}

// Approximate normalized luminance (0-1) for 24-bit RGB
static inline double
rgb_lum(const uby8 rgb[3], const uby8 *imap = NULL)
{
	uby8    orgb[3];
	if (imap != NULL) {
		orgb[0] = imap[rgb[0]];
		orgb[1] = imap[rgb[1]];
		orgb[2] = imap[rgb[2]];
		rgb = orgb;
	}
	return (54*rgb[0]*rgb[0] + 183*rgb[1]*rgb[1] + 18*rgb[2]*rgb[2])*(1./16581375.);
}

// Compute normalized luminance over a 24-bit RGB image
static double
lumAvg(const ImgStruct *ims, const uby8 *imap = NULL)
{
	if (!PmatchColorSpace(ims->csp, &ICS_sRGB, PICMptype))
		return .0;
	double		lumSum = .0;
	for (int y = ims->yres; y-- > 0; ) {
		const uby8 *	bp = ProwPtr(ims, y);
		for (int x = ims->xres; x-- > 0; bp += 3)
			lumSum += rgb_lum(bp, imap);
	}
	return lumSum/((double)ims->xres*ims->yres);
}

// Average luminance over the given rectangle
double
PDisplayList::AvgLuminance(const ImgRect *r)
{
	if (r == NULL)
		return .0;
	if (!GetPointers())
		return .0;
	float		stonits;
	if (l_img.img == NULL) {		// low dynamic-range approx.
		ImgStruct	subimg;
		double		bravg;
		if (!PlinkSubimage(&subimg, &c_img, r)) {
			FreePointers();
			return .0;
		}
		if (!parent->ircd[PDBFstonits].Get(&stonits))
			stonits = 1.f;
		bravg = lumAvg(&subimg, invMap);
		PfreeImage(&subimg);
		FreePointers();
		return stonits*bravg;
	}
	if (!PlegalRect(r, parent->ior.xres, parent->ior.yres)) {
		FreePointers();
		return .0;
	}
	double	lumsum = .0;			// use HDR luminance
	ImgRect	orect = *r;			// get original rectangle
	OrigAspect(&orect.xleft, &orect.ytop);
	OrigAspect(&orect.xright, &orect.ybottom);
	for (int y = orect.ytop; y < orect.ybottom; y++) {
		TMbright *	lp = (TMbright *)PpixPtr(&l_img,orect.xleft,y);
		for (int i = orect.xright - orect.xleft; i-- > 0; )
			lumsum += tmLuminance(*lp++);
	}
	if (parent->ircd[PDBFstonits].Get(&stonits))
		lumsum *= stonits / tms->inpsf; // apply correction
	FreePointers();
	return lumsum / ((double)(orect.xright-orect.xleft) *
				(orect.ybottom-orect.ytop));
}

// Compute histogram over the given rectangle
// Logarithmic divisions for HDR luminance
bool
PDisplayList::GetHistogram(double minmax[2], unsigned long *hist, int hlen,
				PHistoChan hc, const ImgRect *r)
{
	if ((minmax == NULL) | (hist == NULL) | (hlen <= 1))
		return false;
	int     i, j;
	minmax[0] = minmax[1] = .0;
	memset(hist, 0, hlen*sizeof(*hist));
	if (!GetPointers())
		return false;
	ImgRect		myrect;
	ImgStruct       subimg;
	float		ifact;
	if (hc == PHCluminance && l_img.img != NULL) {
						// high dynamic range luminance
		if (r != NULL) {
			myrect = *r;
			OrigAspect(&myrect.xleft, &myrect.ytop);
			OrigAspect(&myrect.xright, &myrect.ybottom);
			r = &myrect;
		}
		if (!PlinkSubimage(&subimg, &l_img, r)) {
			FreePointers();
			return false;
		}
						// compute extrema
		const TMbright  LtooSmall = -50*TM_BRTSCALE;
		TMbright	lmin=16000, lmax=LtooSmall;
		TMbright *      lp;
		for (j = 0; j < subimg.yres; j++) {
			lp = (TMbright *)ProwPtr(&subimg, j);
			for (i = subimg.xres; i--; lp++) {
				if (lmin > *lp && *lp > LtooSmall)
					lmin = *lp;
				if (lmax < *lp)
					lmax = *lp;
			}
		}
		if (lmin >= lmax) {
			PfreeImage(&subimg);
			FreePointers();
			return false;
		}
		minmax[0] = tmLuminance(lmin);
		minmax[1] = tmLuminance(lmax);
		float   stonits;			// make corrections
		if (parent->ircd[PDBFstonits].Get(&stonits)) {
			minmax[0] *= stonits / tms->inpsf;
			minmax[1] *= stonits / tms->inpsf;
		}
		ifact = (hlen - .001)/(lmax - lmin);    // logarithmic histogram
		for (j = 0; j < subimg.yres; j++) {
			lp = (TMbright *)ProwPtr(&subimg, j);
			for (i = subimg.xres; i--; lp++)
				++hist[int(ifact*(*lp - lmin))];
		}
		PfreeImage(&subimg);
		FreePointers();
		return true;
	}
							// get SDR (sub)image
	if (l_img.img != NULL) {
		if (r == NULL) {
			myrect.xleft = myrect.ytop = 0;
			myrect.xright = parent->ior.xres;
			myrect.ybottom = parent->ior.yres;
			r = &myrect;
		}
		subimg.csp = &parent->cs;
		subimg.xres = r->xright - r->xleft;
		subimg.yres = r->ybottom - r->ytop;
		subimg.img = NULL;			// tone-map HDR
		if (!MapSubimage(&subimg, r)) {		// XXX wrong if PEMlocal
			FreePointers();
			return false;
		}
	} else if (!PlinkSubimage(&subimg, &c_img, r)) {
		FreePointers();
		return false;
	}
	int		offs = 0;			// linear histogram
	const uby8 *    bp;
	switch (hc) {
	case PHCluminance:
		minmax[1] = parent->ircd[PDBFstonits].GetFloat();
		if (minmax[1] <= 1e-16)
			minmax[1] = 1.;
		ifact = hlen - .001;
		for (j = 0; j < subimg.yres; j++) {
			bp = ProwPtr(&subimg, j);
			for (i = subimg.xres; i--; bp += 3)
				++hist[int(ifact*rgb_lum(bp,invMap))];
		}
		break;
	case PHCblue:
		++offs;
	    // fall through
	case PHCgreen:
		++offs;
	    // fall through
	case PHCred:
		minmax[1] = 255.;
		ifact = (hlen - .001)/255.;
		for (j = 0; j < subimg.yres; j++) {
			bp = ProwPtr(&subimg, j) + offs;
			for (i = subimg.xres; i--; bp += 3)
				++hist[int(ifact * *bp)];
		}
		break;
	case PHCrgb:
		minmax[1] = 255.;
		ifact = (hlen - .001)/255.;
		for (j = 0; j < subimg.yres; j++) {
			bp = ProwPtr(&subimg, j);
			for (i = 3*subimg.xres; i--; bp++)
				++hist[int(ifact * *bp)];
	}
		break;
	}
	PfreeImage(&subimg);
	FreePointers();
	return true;
}

// Assign/change the display color space
bool
PDisplayImage::SetDisplaySpace(const ImgColorSpace *dcs,
				double lddyn, double ldmax)
{
	if (dcs != NULL) {
		for (PDisplayList *dl = dlhead; dl != NULL; dl = dl->next) {
			if (!PmatchColorSpace(dcs, dl->c_img.csp, PICMall)) {
				if (!PmatchColorSpace(dcs, dl->c_img.csp, PICMptype))
					return false;
				need2reload.Set(dl->frame);
			}
			dl->c_img.csp = &cs;
		}
		cs = *dcs;
		cs.cstatic = true;		// true enough for our purposes
	}
	if (lddyn > 2.)
		dynrange = lddyn;
	if (dlmax > 1.)
		dlmax = ldmax;
	return true;
}

// Reset linear tone-mapping (internal -- assumes GetPointers() has been called)
bool
PDisplayList::ResetLinear(double expmult)
{
	if (expmult <= 0)
		return false;
	bool	ok = true;
	if (tms->histo == NULL)
		ok = RedoHisto();
	if (ok)
		ok = (tmFixedMapping(tms, expmult, c_img.csp->gamma,
					parent->dynrange) == TM_E_OK);
	return ok;
}

// Assign a fixed, linear tone-mapping
bool
PDisplayList::SetFixedLinear(double expmult)
{
	if (!GetPointers())
		return false;
	if (tms == NULL) {
		FreePointers();
		return false;		// only works for HDR images
	}
	tms->flags = parent->newTMflags;
	bool	ok = ResetLinear(expmult);
	FreePointers();
	if ((!ok & !parent->quiet))
		DMESG(DMCwarning, "Cannot set fixed linear tone-mapping");
	return ok;
}

// Get image data pointer(s), reconverting image if necessary
bool
PDisplayList::GetPointers()
{
	if (parent == NULL)
		return false;
	double	expmult = -1;
	if (tms != NULL && ((tms->flags ^ parent->newTMflags) & TM_F_MESOPIC) != 0)
		parent->need2reload.Set(frame);
	if (parent->need2reload.Check(frame))
		if (InUse())
			DMESG(DMCparameter, "Cannot reload image -- active reference");
		else {
			if (tms != NULL && (tms->flags & parent->newTMflags & TM_F_LINEAR))
				expmult = GetMultiplier(tms);
			Free();		// kicks out old object memory & backup
			parent->need2reload.Reset(frame);
		}
	l_img.xres = oxres; l_img.yres = oyres;
	l_img.img = NULL; l_img.mbase = NULL;
	c_img.xres = oxres; c_img.yres = oyres;
	c_img.img = NULL; c_img.mbase = NULL;
	if (GetCacheObject() == NULL)	// get image object
		return false;
	ImgRect	crect;			// get crop subregion
	bool	docrop = parent->CropRect(&crect);
	if (tms != NULL) {		// tone-mapped has TMbright image
		size_t	lumlen = sizeof(TMbright)*oxres*oyres;
		if (l_img.img == NULL) {
			l_img.rowsize = sizeof(TMbright)*oxres;
			l_img.img = (uby8 *)objCache;
			l_img.mbase = NULL;
		}
		if (docrop && !PlinkSubimage(&l_img, &l_img, &crect))
			DMESG(DMCassert, "PlinkSubimage failed!");
		if (objSize > lumlen) {
			if (c_img.img == NULL) {
				c_img.rowsize = sizeof(uby8)*3*oxres;
				c_img.img = (uby8 *)objCache + lumlen;
				c_img.mbase = NULL;
			}
			if (docrop && !PlinkSubimage(&c_img, &c_img, &crect))
				DMESG(DMCassert, "PlinkSubimage failed!");
		}
	} else {			// standard RGB maps to square pixels
		if (c_img.img == NULL) {
			CorrectAspect(&c_img.xres, &c_img.yres);
			c_img.rowsize = ImgPixelSize(c_img.csp)*c_img.xres;
			c_img.img = (uby8 *)objCache;
			c_img.mbase = this;	// only passes retention
		}
		if (docrop) {
			CorrectAspect(&crect.xleft, &crect.ytop);
			CorrectAspect(&crect.xright, &crect.ybottom);
			if (!PlinkSubimage(&c_img, &c_img, &crect))
				DMESG(DMCassert, "PlinkSubimage failed!");
		}
	}
	DASSERT(HolderCount() > 0);
	DASSERT(InUse() > 0);
	if (expmult > 0)		// reset linear exposure value
		ResetLinear(expmult);
	return true;
}

// Get valid crop region for this image -- true if it's not the whole image
bool
PDisplayImage::CropRect(ImgRect *crp) const
{
	ImgRect	dummy;
	int32	cropr[4];
	if (crp == NULL)		// just checking...
		crp = &dummy;
	if (ircd[PDBFcrop].Get(cropr,4) == 4) {
		bool    cropped = false;
		if ((crp->xleft = cropr[0]) > 0)
			cropped = true;
		else
			crp->xleft = 0;
		if ((crp->ytop = cropr[1]) > 0)
			cropped = true;
		else
			crp->ytop = 0;
		if ((crp->xright = cropr[2]) < GetOWidth())
			cropped = true;
		else
			crp->xright = GetOWidth();
		if ((crp->ybottom = cropr[3]) < GetOHeight())
			cropped = true;
		else
			crp->ybottom = GetOHeight();
		if ((crp->xleft < crp->xright) & (crp->ytop < crp->ybottom))
			return cropped;
	}
	crp->xleft = crp->ytop = 0;
	crp->xright = GetOWidth(); crp->ybottom = GetOHeight();
	return false;
}

// Get corresponding original rectangle
bool
PDisplayImage::GetOrigRect(ImgRect *orp, const ImgRect *r)
{
	if (orp == NULL)
		return false;
	if (!CheckLoad())
		return false;
	if (r == NULL) {
		CropRect(orp);
		return true;
	}
	ImgRect	cr;
	CropRect(&cr);
	*orp = *r;
	OrigAspect(&orp->xleft, &orp->ytop);
	orp->xleft += cr.xleft; orp->ytop += cr.ytop;
	OrigAspect(&orp->xright, &orp->ybottom);
	orp->xright += cr.xleft; orp->ybottom += cr.ytop;
	return true;
}

// Get corresponding view rectangle
bool
PDisplayImage::GetViewRect(ImgRect *vrp, const ImgRect *r)
{
	if (vrp == NULL)
		return false;
	if (!CheckLoad())
		return false;
	if (r == NULL) {
		vrp->xleft = vrp->ytop = 0;
		vrp->xright = ior.xres;
		vrp->ybottom = ior.yres;
		return true;
	}
	ImgRect	cr;
	CropRect(&cr);
	*vrp = *r;
	vrp->xleft -= cr.xleft; vrp->ytop -= cr.ytop;
	CorrectAspect(&vrp->xleft, &vrp->ytop);
	vrp->xright -= cr.xleft; vrp->ybottom -= cr.ytop;
	CorrectAspect(&vrp->xright, &vrp->ybottom);
	return true;
}

// Open image file and load into memory
bool
PDisplayList::NewOriginal()
{
	if (parent == NULL)
		return false;
	DASSERT(objCache == NULL);
	ImgReader *	ir = parent->ir;
	if (ir == NULL) {		// make sure we have open image
		ir = PopenImageD(parent->ircd, parent->quiet);
		if (ir == NULL)
			return false;
		if (frame > 0) {	// trying to seek to frame
			if (ir->nframes == 1) {
				IRclose(ir);
				return false;
			} else		// keep open if seeking frame
				parent->ir = ir;
		}
	}
					// update frame information
	parent->frameType = ir->frameType;
	parent->frameRate = ir->frameRate;
					// go to this frame number
	if (frame != ir->frame && IRseekFrame(ir, frame, IRSabs) != IREnone) {
		if (!parent->quiet)
			PreportReaderErr(ir);
		if ((parent->nframes = ir->nframes) == 1)
			parent->ir = NULL;
		if (ir != parent->ir)
			IRclose(ir);
		return false;
	}
					// get corrected image size
	oxres = ir->xres; oyres = ir->yres;
	pixAspect = ir->pixAspect;
	parent->UpdateRec(parent->ircd);
	if (invMap != NULL) {
		delete [] invMap;
		invMap = NULL;
	}
	bool	ok;
	switch (ir->cs.dtype) {
	case IDTubyte:			// 8-bit gray or 24-bit RGB
		if ((ir->cs.format != IPFrgb) & (ir->cs.format != IPFy))
			goto unsupported;
		c_img.xres = ir->xres; c_img.yres = ir->yres;
		CorrectAspect(&c_img.xres, &c_img.yres);
		if (!PsizeOK(c_img.xres, c_img.yres,
				c_img.rowsize = ImgPixelSize(c_img.csp))) {
			DMESG(DMCparameter, "Image too big for memory");
			return false;
		}
		objSize = (ssize_t)c_img.yres * (c_img.rowsize *= c_img.xres);
		objCache = Cmalloc(objSize);
		if (objCache == NULL) {
			objSize = 0;
			return false;
		}
		c_img.img = (uby8 *)objCache;
		c_img.mbase = NULL;
		ok = PrenderImageR(&c_img, ir, parent->quiet);
		if (c_img.img != (uby8 *)objCache) {
			PfreeImage(&c_img);
			ok = false;
		}
		if (ok)
			AdjustExposure();
		if (ok && ir->cs.format == IPFrgb)
			RemoveRedEye();
		lossy = PmatchSuffix("jpg", ir->ri->suffixList);
		if ((parent->nframes = ir->nframes) == 1)
			parent->ir = NULL;
		if (ir != parent->ir)
			IRclose(ir);
		if (!ok) {
			Cfree(objCache); 
			objCache = NULL; objSize = 0;
			return false;
		}
		c_img.mbase = this;	// only passes retention
		return true;
	case IDTfloat:			// gray or color HDR
	case IDTushort:
		if (parent->quiet)
			parent->newTMflags |= TM_F_NOSTDERR;
		else
			parent->newTMflags &= ~TM_F_NOSTDERR;
		switch (ir->cs.format) {
		case IPFrgb:
		case IPFxyz:
			parent->newTMflags &= ~TM_F_BW;
			break;
		case IPFy:
			parent->newTMflags |= TM_F_BW;
			break;
		default:
			goto unsupported;
		}
		ok = LoadHDR(ir);
		if (ok && ir->cs.format == IPFrgb)
			RemoveRedEye();
		lossy = true;		// should be false if IEEE float?
		if ((parent->nframes = ir->nframes) == 1)
			parent->ir = NULL;
		if (ir != parent->ir)
			IRclose(ir);
		return ok;
	default:
		break;
	}
unsupported:				// unsupported format
	if (!parent->quiet)
		DMESG(DMCdata, "Unsupported image data format");
	parent->nframes = ir->nframes;
	IRclose(ir);
	parent->ir = NULL;
	return false;
}

// Convert original rectangle to local one in image (internal)
bool
PDisplayList::LocalRect(ImgRect *rl, const int32 rect[4]) const
{
	int     xr, yr;
	ImgRect	crect;
					// get internal image res.
	if (c_img.img != NULL) {
		xr = c_img.xres; yr = c_img.yres;
	} else {
		xr = l_img.xres; yr = l_img.yres;
	}
	rl->xleft = rect[0];
	rl->ytop = rect[1];
	rl->xright = rect[2];
	rl->ybottom = rect[3];
					// adjust coordinates for crop
	if (parent->CropRect(&crect) &&
			crect.xright - crect.xleft == xr &&
			crect.ybottom - crect.ytop == yr) {
		if (!PclipRect(rl, &crect))
			return false;
		rl->xleft -= crect.xleft;
		rl->ytop -= crect.ytop;
		rl->xright -= crect.xleft;
		rl->ybottom -= crect.ytop;
	}
	if (tms == NULL) {		// mapped to square pixels
		CorrectAspect(&rl->xleft, &rl->ytop);
		CorrectAspect(&rl->xright, &rl->ybottom);
	}
	return PlegalRect(rl, xr, yr);
}

// Recompute histogram, accounting for spot-metered regions (internal)
bool
PDisplayList::RedoHisto(const ImgRect *ro)
{
	ImgRect	hrect;
	int     y;

	if (ro == NULL) {		// use entire cropped region?
		hrect.xleft = hrect.ytop = 0;
		hrect.xright = l_img.xres;
		hrect.ybottom = l_img.yres;
		ro = &hrect;
	}
	tmClearHisto(tms);		// start with overall
	for (y = ro->ytop; y < ro->ybottom; y++)
		if (tmAddHisto(tms, (TMbright *)PpixPtr(&l_img,ro->xleft,y),
				ro->xright - ro->xleft, 1) != TM_E_OK)
			return false;
					// add in spot exposure(s)
	int32   spotexpr[32];
	int     nv = parent->ircd[PDBFspotexp].Get(spotexpr, 32);
	if (!nv)
		return true;		// no spot meterings
	if (nv%4 != 0) {
		DMESG(DMCwarning, "Bad spot exposure parameter count");
		return true;
	}
	nv /= 4;
	const int       ngeneral = PrectArea(ro);
	int		nspot = 0;
	int		i;
	for (i = 0; i < nv; i++) {
		if (!LocalRect(&hrect, spotexpr+4*i))
			continue;
		if (!PclipRect(&hrect, ro))
			continue;
		nspot += PrectArea(&hrect);
	}
	if (!nspot)
		return true;		// empty or clipped
	int     sweight = ngeneral/nspot;
	if (sweight > 200)
		sweight = 200;
	for (i = 0; i < nv; i++) {      // weigh in spots
		if (!LocalRect(&hrect, spotexpr+4*i))
			continue;
		if (!PclipRect(&hrect, ro))
			continue;
		for (y = hrect.ytop; y < hrect.ybottom; y++)
			tmAddHisto(tms, (TMbright *)PpixPtr(&l_img,hrect.xleft,y),
					hrect.xright - hrect.xleft, sweight);
	}
	return true;
}

// Map pixels in a subimage using our tone-mapping operator (internal)
bool
PDisplayList::MapSubimage(ImgStruct *ims, const ImgRect *r, double *localexp)
{
	if ((parent == NULL) | (ims == NULL) | (l_img.img == NULL))
		return false;
	if (!PlegalRect(r, parent->ior.xres, parent->ior.yres))
		return false;
	if (tms == NULL)
		return false;
	ImgRect	orect = *r;			// get original rectangle
	OrigAspect(&orect.xleft, &orect.ytop);
	OrigAspect(&orect.xright, &orect.ybottom);
						// allocate destination
	/*
	if ( !PnewImage( ims, double(r->ybottom-r->ytop) /
					(r->xright-r->xleft) ) ) {
	*/
	if (!PnewImage(ims, .0))
		return false;
	ImgStruct *	oims = ims;
	ImgStruct	myImg;
	if ((ims->xres != orect.xright - orect.xleft) |
			(ims->yres != orect.ybottom - orect.ytop)) {
		myImg.csp = ims->csp;		// will need to resample
		myImg.xres = orect.xright - orect.xleft;
		myImg.yres = orect.ybottom - orect.ytop;
		myImg.img = NULL;
		if (!PnewImage(&myImg, .0)) {
			PfreeImage(ims);
			return false;
		}
		oims = &myImg;
	}
	int	y;
	bool	newhist = (tms->histo == NULL) | (localexp != NULL);
	if (newhist)				// recompute histogram?
		if (!RedoHisto(localexp!=NULL ? &orect : (ImgRect *)NULL))
			goto fail;
						// recompute tone-mapping?
	if (newhist | (tms->flags != parent->newTMflags) | (tms->lumap == NULL)) {
		tms->flags = parent->newTMflags;
		if (tmComputeMapping(tms, c_img.csp->gamma,
				parent->dynrange, parent->dlmax) != TM_E_OK)
			goto fail;
	}
						// automatic false color scale?
	if ((parent->fcs != NULL) & (parent->fcManual <= 0)) {
		parent->fcManual = false;
		if (tms->flags & TM_F_LINEAR)
			fcLinearMapping(parent->fcs, tms, 1);
		else
			fcLogMapping(parent->fcs, tms, 0.25);
	}
						// tone-map pixels
	for (y = orect.ytop; y < orect.ybottom; y++) {
		uby8 *	dptr = ProwPtr(oims, y-orect.ytop);
		if (parent->fcs != NULL) {
			if (fcMapPixels(parent->fcs, dptr,
					(TMbright *)PpixPtr(&l_img,orect.xleft,y),
					orect.xright-orect.xleft) != TM_E_OK)
				goto fail;
			continue;
		}
		if (tmMapPixels(tms, dptr,
				(TMbright *)PpixPtr(&l_img,orect.xleft,y),
				(c_img.img==NULL ? TM_NOCHROM
					: PpixPtr(&c_img,orect.xleft,y)),
				orect.xright-orect.xleft) != TM_E_OK)
			goto fail;
		if (c_img.img == NULL && oims->csp->format != IPFy)
			PgetRGB24fromY8(dptr, dptr, orect.xright-orect.xleft);
	}
	if (oims != ims) {			// resample final output
		if (!PsizeImage(ims, oims, PSbest))
			goto fail;
		PfreeImage(oims);
		oims = ims;
	}
	if (localexp != NULL) {
		*localexp = GetMultiplier(tms);
		tmClearHisto(tms);		// reset histogram if local
	}
	return true;
fail:						// land here on failure
	if (!parent->quiet)
		DMESG(DMCdata, "Cannot map HDR pixels in MapSubimage");
	if (localexp != NULL)
		tmClearHisto(tms);
	if (oims != ims)
		PfreeImage(oims);
	PfreeImage(ims);
	return false;
}

// Fix exposure on a 24-bit RGB image using spot metering
void
PDisplayList::AdjustExposure()
{
	if (tms != NULL || c_img.img == NULL ||
			!PmatchColorSpace(c_img.csp, &ICS_sRGB, PICMptype))
		return;				// can't handle these
	int32   spotexpr[32];
	int     nv = parent->ircd[PDBFspotexp].Get(spotexpr, 32);
	if (!nv)
		return;				// no spot metering
	if (nv%4 != 0) {
		DMESG(DMCwarning, "Bad spot exposure parameter count");
		return;
	}
	nv /= 4;
	uby8		extrem8[2] = {0, 255};
	unsigned long	pmap[256];
	int		i, x, y;
	uby8 *		bp;
						// compute histogram
	memset(pmap, 0, sizeof(pmap));
	PcomputeHisto(extrem8, pmap, 256, &c_img, PHCrgb);
	int		sum = 0;		// average over spots
	int		nsp = 0;
	for (i = 0; i < nv; i++) {
		ImgRect sr;
		if (!LocalRect(&sr, spotexpr+4*i))
			continue;
		for (y = sr.ytop; y < sr.ybottom; y++) {
			bp = PpixPtr(&c_img, sr.xleft, y);
			for (x = sr.xright-sr.xleft; x-- > 0; ) {
				sum += (54*bp[0] + 183*bp[1] + 18*bp[2]) >> 8;
				pmap[*bp++] += 50;
				pmap[*bp++] += 50;
				pmap[*bp++] += 50;
			}
		}
		nsp += PrectArea(&sr);
	}
	if (nsp <= 0)
		return;
	x = c_img.xres*c_img.yres*3/100;	// 1% allowed over max.
	for (i = 256; i--; )
		if ((x -= pmap[i]) <= 0)
			break;
	int		maxv = i;
	const int       avg = sum / nsp;
	if ((15 >= avg) | (avg >= 235))
		return;
	if (maxv <= avg)
		maxv = avg + 1;
	int		goal = 118 - (118 - avg)*5/12;
	if (goal < avg*255/maxv)
		goal = avg*255/maxv;
	/*
	 * Compute Reinhard global TMO such that
	 * pmap[avg] == goal & pmap[maxv] == 255.
	 */
	const double    K = ((double)goal/avg - 255.*avg/(maxv*maxv)) /
					(255. - goal);
	for (i = 0; i < maxv; i++)
		pmap[i] = int( .5 + 255.*K*i*(1. + i/(K*maxv*maxv)) /
						(1. + K*i) );
	for (i = maxv; i < 256; i++)
		pmap[i] = 255;
						// apply adjustment
	for (y = 0; y < c_img.yres; y++) {
		bp = ProwPtr(&c_img, y);
		for (x = 3*c_img.xres; x-- > 0; bp++)
			*bp = pmap[*bp];
	}
						// record inverse mapping
	if (invMap == NULL)
		invMap = new uby8 [256];
	memset(invMap, 0, sizeof(uby8)*256);
	for (i = 0; i < 256; i++)
		invMap[pmap[i]] = i;
	for (i = 1; i < 256; i++)
		if (!invMap[i])
			invMap[i] = invMap[i-1];
}

// Remove red-eye from color image
void
PDisplayList::RemoveRedEye(const DBRecord *drp)
{
	int32	redeyer[100];
	int	nv;
	if (drp == &parent->ircd)
		drp = NULL;
	int32 *	rer = redeyer;
	if (drp == NULL) {			// remove all
		nv = parent->ircd[PDBFredeye].Get(redeyer, 100);
	} else {				// else last rect(s) are new
		nv = parent->ircd[PDBFredeye].GetNV();
		rer += nv;
		nv = PDBgetField(*drp,PDBFredeye).Get(redeyer, 100) - nv;
	}
	if (nv%4 != 0) {
		DMESG(DMCwarning, "Bad red-eye parameter count");
		return;
	}
	for (nv /= 4; nv-- > 0; rer += 4)
		RemoveRedEye(rer);
}

// Heuristic check if pixel is characteristic red-eye color
static inline bool
isRedEye(const uby8 rgb[3])
{
	return (rgb[0] >= 110 && 3*rgb[0] >= 4*(rgb[1] + rgb[2]));
}

// Looser check if pixel is legal part of red-eye pupil
static inline bool
inRedEye(const uby8 rgb[3])
{
	if (rgb[0] >= 230)			// highlight check
		return true;
	return (rgb[0] >= 60 && rgb[0] >= (rgb[1] + rgb[2]));
}

// Check if pixel is legal part of eye white
static inline bool
inEyeWhite(const uby8 rgb[3])
{
	if (rgb[0] < 60 || rgb[1] < 60 || rgb[2] < 60)
		return false;
	float	rat;
	rat = (float)rgb[0] / (float)rgb[1];
	if (!(0.75 <= rat) & (rat <= 1.33))
		return false;
	rat = (float)rgb[2] / (float)rgb[1];
	return (0.75 <= rat) & (rat <= 1.33); 
}

// Perform flood fill of red-eye pixels
static int
floodRedEye(ABitMap2 *pbm, const ImgStruct *ims, int x, int y, int depth)
{
	if (depth <= 0)
		return 0;
	if ((x < 0) | (x >= ims->xres))
		return 0;
	if ((y < 0) | (y >= ims->yres))
		return 0;
	if (pbm->Check(x, y))
		return 0;			// already flooded
	if (!inRedEye(PpixPtr(ims, x, y)))
		return 0;			// not in red-eye pupil
	pbm->Set(x, y);				// else paint and flood neighbors
	--depth;
	return 1 + floodRedEye(pbm, ims, x+1, y, depth) +
			floodRedEye(pbm, ims, x, y+1, depth) +
			floodRedEye(pbm, ims, x-1, y, depth) +
			floodRedEye(pbm, ims, x, y-1, depth);
}

#undef rad2	// *$#! Windows
// Collect pupil around its center of mass and check for whites
static bool
goodPupil(ABitMap2 *pbm, const ImgStruct *ims)
{
	int     np = 0;				// compute center of mass
	int     x0 = 0, y0 = 0;
	int	i, j;
	for (i = j = 0; pbm->Find(&i, &j); i++)
		if (isRedEye(PpixPtr(ims, i, j))) {
			x0 += i; y0 += j;
			++np;
		}
	int	area = pbm->SumTotal();		// check for wild flood
	if (area >= 10*np)
		return false;
	x0 /= np; y0 /= np;
						// apply circular mask
	int		rad = (int)(1.4*sqrt(area/3.14) + 1.5);
	int		rad2 = rad*rad;
	for (j = 0; j < pbm->Height(); j++)
		for (i = 0; i < pbm->Width(); i++)
			if ((i-x0)*(i-x0) + (j-y0)*(j-y0) > rad2)
				pbm->Reset(i, j);
	pbm->Expand(1);				// smudge it around a bit
	area = pbm->SumTotal();			// revise pupil size
	if (!area)
		return false;
	rad = (int)(sqrt(area/3.14) + .5);
	rad2 = rad*rad;
						// check for nearby whites of eye
	const int	erad = 5*rad;		// count white pixels in 5*rad
	const int	erad2 = erad*erad;
	int		nwhite = 0;
	for (j = 0; j < ims->yres; j++) {
		if (j < y0-erad) j = y0-erad;
		if (j > y0+erad) break;
		for (i = 0; i < ims->xres; i++) {
			if (i < x0-erad) i = x0-erad;
			if (i > x0+erad) break;
			if ((i-x0)*(i-x0) + (j-y0)*(j-y0) > erad2)
				continue;
			if (pbm->Check(i, j))
				continue;
			nwhite += inEyeWhite(PpixPtr(ims, i, j));
		}
	}
	return (nwhite >= area/4);		// whites at least 1/4 pupil area
}

// Find red-eye pupils in subimage
static int
PfindRedEye(ABitMap2 *bm, const ImgStruct *ims)
{
						// only works for RGB
	DASSERT(PmatchColorSpace(ims->csp,&ICS_sRGB,PICMptype));
						// clear bitmap
	bm->NewBitMap(ims->xres, ims->yres);
						// assign min. pupil diameter
	int	mind = 3 * (ims->xres > ims->yres ?
					ims->xres : ims->yres) / 100;
	if (mind < 2) mind = 2;
						// find pupils
	for (int y = 0; y < ims->yres; y += mind/2)
		for (int x = 0; x < ims->xres; x += mind/2) {
			if (bm->Check(x, y))
				continue;	// already in pupil
			const uby8 *	rgb = PpixPtr(ims, x, y);
			if (!isRedEye(rgb))
				continue;	// not red enough
			ABitMap2	pbm(ims->xres, ims->yres);
			int		maxd = (ims->xres < ims->yres ?
						ims->xres : ims->yres);
			if (floodRedEye(&pbm, ims, x, y, maxd) < mind*mind)
				continue;
			if (!goodPupil(&pbm, ims))
				continue;
			if (!pbm.Check(x, y))
				continue;
			*bm |= pbm;		// add pupil to bitmap
		}
	return bm->SumTotal();
}

// Remove red-eye from rectangle in color image layer
void
PDisplayList::RemoveRedEye(const int32 rect[4])
{
	if (c_img.img == NULL)
		return;
	ImgStruct	subimg;
	ImgRect		rct;
	ABitMap2	rebm;
	if (!LocalRect(&rct, rect))		// map to image coord's
		return;
	if (5*5*PrectArea(&rct) > parent->ior.xres*parent->ior.yres) {
		DMESG(DMCwarning, "Region too large for red-eye removal");
		return;
	}
	if (!PlinkSubimage(&subimg, &c_img, &rct))
		return;
	if (!PfindRedEye(&rebm, &subimg)) {
		PfreeImage(&subimg);
		return;
	}
	for (int y = 0; y < subimg.yres; y++) {	// darkest gray in pupils
		uby8 *	bptr = ProwPtr(&subimg, y);
		for (int x = 0; x < subimg.xres; x++, bptr += 3) {
			if (!rebm.Check(x, y))
				continue;
			int	amt = 1 + rebm.Check(x-1,y) +
					rebm.Check(x,y-1) +
					rebm.Check(x+1,y) +
					rebm.Check(x,y+1);
			if (bptr[2] < bptr[1]) {
				bptr[0] = ((5-amt)*bptr[0] + amt*bptr[2])/5;
				bptr[1] = ((5-amt)*bptr[1] + amt*bptr[2])/5;
			} else {
				bptr[0] = ((5-amt)*bptr[0] + amt*bptr[1])/5;
				bptr[2] = ((5-amt)*bptr[2] + amt*bptr[1])/5;
			}
		}
	}
	PfreeImage(&subimg);
	FreeFile();				// disk cache is invalid
}

// Allocate and load a high dynamic-range image and set up tone-mapping
bool
PDisplayList::LoadHDR(ImgReader *ir)
{
	if (parent == NULL || ir == NULL ||
			(ir->cs.dtype != IDTfloat) & (ir->cs.dtype != IDTushort))
		return false;
	if (tms != NULL)			// free previous tone-mapping
		tmDone(tms);
						// set output color space
	DASSERT(c_img.csp != NULL);
	tms = tmInit(parent->newTMflags, const_cast<RGBPRIMP>(c_img.csp->chroma),
				c_img.csp->gamma);
	if (tms == NULL) {
		if (!parent->quiet)
			DMESGF(DMCdata, "Cannot set tone-mapping for '%s'",
					ir->file);
		return false;
	}
						// check if buffer will fit
	if (!PsizeOK(ir->xres, ir->yres, (tms->flags & TM_F_BW) ?
			sizeof(TMbright) : sizeof(TMbright)+sizeof(uby8)*3)) {
		DMESG(DMCparameter, "Image too big to fit in memory");
		return false;
	}
	const size_t	lumlen = sizeof(TMbright)*ir->xres*ir->yres;
	RGBPRIMP	inchroma = (RGBPRIMP)ir->cs.chroma;
	ImgReadBuf	irb;
	ImgInfo		iinfo;
	if (LoadHDRspecial(ir))			// try optimized HDR loader
		return true;			// success!
	iinfo.stonits = 1.f;
	if (!parent->ircd[PDBFstonits].Get(&iinfo.stonits) &&
			IRgetInfo(ir, &iinfo) == IREnone &&
			iinfo.flags & IIFstonits)
		parent->ircd.SetField(PDBFstonits, iinfo.stonits);
	if (ir->cs.format == IPFxyz)
		inchroma = TM_XYZPRIM;
	if (tmSetSpace(tms, inchroma, iinfo.stonits) != TM_E_OK) {
		if (!parent->quiet)
			DMESGF(DMCdata, "Cannot set color space for '%s'",
					ir->file);
		tmDone(tms);
		tms = NULL;
		return false;
	}
						// allocate storage
	if (tms->flags & TM_F_BW)
		objSize = lumlen;
	else
		objSize = lumlen + sizeof(uby8)*3*ir->xres*ir->yres;
	objCache = Cmalloc(objSize);		// allocate cache
	if (objCache == NULL) {
		objSize = 0;
		return false;
	}
						// read the image
	l_img.xres = ir->xres; l_img.yres = ir->yres;
	l_img.rowsize = sizeof(TMbright)*ir->xres;
	l_img.img = (uby8 *)objCache;
	l_img.mbase = NULL;
	if (objSize > lumlen) {
		c_img.xres = ir->xres; c_img.yres = ir->yres;
		c_img.rowsize = sizeof(uby8)*3*ir->xres;
		c_img.img = (uby8 *)objCache + lumlen;
		c_img.mbase = NULL;
	}
						// set up read buffer
	irb.cs = ir->cs;
	irb.r = ir->fr;				// start at the beginning
	irb.subsample = 1;
	irb.buf = NULL;
	do {					// read each strip
		if (IRreadRec(ir, &irb) != IREnone) {
			if (!parent->quiet)
				PreportReaderErr(ir);
			goto bigfailure;
		}
		for (int y = irb.r.ytop; y < irb.r.ybottom; y++) {
			size_t	offset = (ssize_t)y*ir->xres + irb.r.xleft;
			int	slen = irb.r.xright - irb.r.xleft;
			int	tmres = TM_E_CODERR1;
			if (irb.cs.dtype == IDTfloat) {
				if (irb.cs.format == IPFy)
					tmres = tmCvGrays(tms, (TMbright *)l_img.img + offset,
						(float *)irb.buf + size_t(y-irb.r.ytop)*slen,
						slen);
				else
					tmres = tmCvColors(tms, (TMbright *)l_img.img + offset,
						(c_img.img==NULL ? TM_NOCHROM
							: c_img.img + 3*offset),
						(COLOR *)irb.buf + size_t(y-irb.r.ytop)*slen,
						slen);
			} else /* irb.cs.dtype == IDTushort */ {
				if (irb.cs.format == IPFy)
					tmres = tmCvGray16(tms, (TMbright *)l_img.img + offset,
						(uint16 *)irb.buf + size_t(y-irb.r.ytop)*slen,
						slen, irb.cs.gamma);
				else
					tmres = tmCvRGB48(tms, (TMbright *)l_img.img + offset,
						(c_img.img==NULL ? TM_NOCHROM
							: c_img.img + 3*offset),
						(uint16 (*)[3])irb.buf + size_t(y-irb.r.ytop)*slen,
						slen, irb.cs.gamma);
			}
			if (tmres != TM_E_OK) {
				if (!parent->quiet)
					DMESGF(DMCdata,
					"Cannot map luminances in HDR image '%s'",
						ir->file);
				goto bigfailure;
			}
		}
		irb.r = ir->nr;			// advance to next strip
	} while (IRmoreRec(ir));
	free(irb.buf);				// image converted!
	irb.buf = NULL;
	return true;				// success!
bigfailure:					// land here on failure
	l_img.img = NULL;
	c_img.img = NULL;
	if (objCache != NULL)
		Cfree(objCache);
	objCache = NULL;
	objSize = 0;
	if (irb.buf != NULL)
		free(irb.buf);
	tmDone(tms);
	tms = NULL;
	return false;
}

// Allocate and load an HDR image using optimized library routine
bool
PDisplayList::LoadHDRspecial(ImgReader *ir)
{
	if (ir->frame)					// not for animations
		return false;
	if (ir->cs.logorig > 0)				// no general log input
		return false;
	const size_t	lumlen = sizeof(TMbright)*ir->xres*ir->yres;
	const size_t	chromlen = sizeof(uby8)*3*ir->xres*ir->yres;
	int		xr, yr;
	if (ir->ri == &IRInterfaceRad) {		// try Radiance routine
		int	sl, ns;
		FILE *	fp = fopen(ir->file, "rb");
		if (fp == NULL)
			return false;
		getheader(fp, NULL, NULL);		// check scan ordering
		if (!fscnresolu(&sl, &ns, fp)) {
			fclose(fp);			// non-standard order
			return false;
		}
		CacheMakeRoom(lumlen+chromlen);		// make space for it
		fseek(fp, 0L, 0);			// load/convert picture
		if (tmLoadPicture(tms, (TMbright **)&l_img.img, &c_img.img,
				&xr, &yr, ir->file, fp) != TM_E_OK) {
			fclose(fp);
			return false;
		}
		fclose(fp);
	} else if (ir->ri == &IRInterfaceTIFF &&
			((ir->cs.dtype == IDTushort && (ir->cs.gamma >= 2.1) &
					(ir->cs.gamma <= 2.3)) ||
			!strcmp(ir->encoding, "RLE-Log Luv") ||
			!strcmp(ir->encoding, "24-bit-Log Luv") ||
			!strcmp(ir->encoding, "RLE-Log luminance"))) {
		CacheMakeRoom(lumlen+chromlen);		// make space & load it
		if (tmLoadTIFF(tms, (TMbright **)&l_img.img, &c_img.img,
				&xr, &yr, ir->file, NULL) != TM_E_OK)
			return false;
	} else
		return false;
	DASSERT((xr==ir->xres) & (yr==ir->yres));
	l_img.xres = xr; l_img.yres = yr;
	l_img.rowsize = sizeof(TMbright)*xr;
	l_img.mbase = NULL;
	if (c_img.img == NULL) {			// just luminance
		objCache = (void *)l_img.img;
		objSize = lumlen;
		return true;
	}
							// else do chroma also
	objCache = Cmalloc(objSize = lumlen+chromlen);
	if (objCache == NULL) {
		free(c_img.img); c_img.img = NULL;
		free(l_img.img); l_img.img = NULL;
		objSize = 0;
		return false;
	}
	memcpy((uby8 *)objCache, l_img.img, lumlen);
	free(l_img.img);
	l_img.img = (uby8 *)objCache;
	memcpy((uby8 *)objCache+lumlen, c_img.img, chromlen);
	free(c_img.img);
	c_img.img = (uby8 *)objCache + lumlen;
	c_img.xres = xr; c_img.yres = yr;
	c_img.rowsize = sizeof(uby8)*3*xr;
	c_img.mbase = NULL;
	return true;				// success!
}
