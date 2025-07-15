/*
 *  pdispimg.h
 *  panlib
 *
 *  Convenience classes for high dynamic range image display.
 *
 *  Depends on pancine.h
 *
 *  Created by gward on Tue Jan 08 2002.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#ifndef _PDISPIMG_H_
#define _PDISPIMG_H_

#ifndef _PANCINE_H_
#include "pancine.h"
#endif
#include "color.h"
#include "tonemap.h"
#include "falsecolor.h"

#define MAXFRAMES		10000	// maximum number of frames

class PDisplayImage;			// forward declaration

// List class for holding image frames in a sequence (can be HDR)
class PDisplayList : public DiskCacheObject {
	friend class	PDisplayImage;
protected:
	ImgStruct	l_img;		// luminance image struct
	ImgStruct	c_img;		// RGB color image struct
	int		oxres, oyres;	// original image dimensions
	float		pixAspect;	// original pixel aspect ratio
	TMstruct *	tms;		// tone-mapping structure for HDR
	uby8 *		invMap;		// pixel inverse mapping
	bool		lossy;		// conversion to RGB loses information
	PDisplayImage *	parent;		// parent object
	int		frame;		// frame number
	PDisplayList *	next;		// next frame in our list
	static double	GetMultiplier(TMstruct *tm);
	bool		LocalRect(ImgRect *rl, const int32 rect[4]) const;
	bool		RedoHisto(const ImgRect *ro = NULL);
	bool		ResetLinear(double exposure);
	bool		LoadHDR(ImgReader *ir);
	bool		LoadHDRspecial(ImgReader *ir);
	bool		MapSubimage(ImgStruct *ims, const ImgRect *r,
					double *localexp = NULL);
	void		AdjustExposure();
	void		RemoveRedEye(const DBRecord *drp = NULL);
	void		RemoveRedEye(const int32 rect[4]);
	bool		GetPointers();
	void		FreePointers() {
				ReleaseObject(objCache);
			}
	void		Empty() {
				l_img.csp = &ICS_Y16;
				l_img.mbase = NULL; l_img.img = NULL;
				c_img.csp = &ICS_sRGB;
				c_img.mbase = NULL; c_img.img = NULL;
				oxres = oyres = 0; pixAspect = 1.f;
				tms = NULL; invMap = NULL; lossy = false;
				frame = 0; parent = NULL; next = NULL;
			}
	virtual bool	NewOriginal();
	virtual void	FreeOriginal() {
				Cfree(objCache); objCache = NULL;
				l_img.img = c_img.img = NULL;
			}
public:
			PDisplayList(PDisplayImage *par, int fno = 0);
	virtual		~PDisplayList() {
				if (next != NULL) next->HolderRelease();
				if (tms != NULL) tmDone(tms);
				delete [] invMap;
				Free();
			}
			// Get name of subclass
	virtual const char *
			ClassName() const {
				return "PDisplayList";
			}
			// Replacement for delete operator
	virtual void	HolderRelease() {
				if ((next != NULL) & (HolderCount() == 1))
					{ next->HolderRelease(); next = NULL; }
				DiskCacheObject::HolderRelease();
			}
			// Make sure we have the correct image information
	bool		CheckLoad();
			// Correct coordinates to square aspect pixels
	void		CorrectAspect(int *xp, int *yp) const {
				if (pixAspect < 0.995f)
					*xp = int(*xp*pixAspect + .5f);
				else if (pixAspect > 1.005f)
					*yp = int(*yp/pixAspect + .5f);
			}
			// Corresponding position in original image
	void		OrigAspect(int *xp, int *yp) const {
				if (pixAspect < 0.995f) {
					*xp = int(*xp/pixAspect + .5f);
					if (*xp > oxres) *xp = oxres;
				} else if (pixAspect > 1.005f) {
					*yp = int(*yp*pixAspect + .5f);
					if (*yp > oyres) *yp = oyres;
				}
			}
			// Set fixed linear mapping using the given multiplier
	bool		SetFixedLinear(double expmult = 1.);
			// Get the current multiplier (approx. if non-linear)
	double		GetMultiplier() {
				if (!CheckLoad()) return .0;
				return GetMultiplier(tms);
			}
			// These calls work in concert with PfreeImage()
	bool		GetSubimage(ImgStruct *ims, const ImgRect *r,
					double *localexp = NULL);
	bool		GetImage(ImgStruct *ims) {
				return GetSubimage(ims, NULL);
			}
			// Average luminance over the given rectangle
	double		AvgLuminance(const ImgRect *r);
			// Compute histogram over the given rectangle
			// Logarithmic divisions for HDR luminance
	bool		GetHistogram(double minmax[2],
					unsigned long *hist, int hlen,
					PHistoChan hc = PHCluminance,
					const ImgRect *r = NULL);
};

// Class for caching HDR image and converting to display RGB (only)
// Does not reorient image -- use ior to figure out how!
class PDisplayImage {
	friend class	PDisplayList;
private:
	DBFieldInfo	finf;		// record field info holder
	DBRecord	ircd;		// image database record
	int32		oxres, oyres;	// original image dimensions
	float		pixAspect;	// original pixel aspect ratio
	int		nframes;	// ir->nframes
	ImgFrameType	frameType;	// ir->frameType
	float		frameRate;	// ir->frameRate
	ImgReader *	ir;		// open image reader
	PDisplayList *	dlhead;		// display list head
	PDisplayList *	dlp;		// current display image
	ImgColorSpace	cs;		// display color space
	double		dynrange;	// display dynamic range
	double		dlmax;		// display maximum luminance
	int		newTMflags;	// flags to use in next tone-mapping
	ABitMap		need2reload;	// frame(s) need to be reloaded
	void		SetRecord(const DBRecord &pr);
	bool		SetFrame(int fn);
	void		Empty() {
				nframes = 1; frameType = IRFnone; frameRate = 0;
				ir = NULL;
				oxres = oyres = 0; pixAspect = 1.f;
				dlhead = dlp = NULL;
				cs = ICS_sRGB;
				dynrange = 100.; dlmax = 100.;
				newTMflags = TM_F_CAMERA;
				need2reload.NewBitMap(MAXFRAMES);
				fcs = NULL; fcManual = -1;
				quiet = false;
				cloz = NULL;
			}
public:
	FCstruct	*fcs;		// false color mapping
	int		fcManual;	// manual false color scale?
	POrient		ior;		// cropped image size and orientation
	bool		quiet;		// don't report errors
					// user callback to cleanup before closing
	bool		(*cloz)(PDisplayImage *, bool cancancel);
			PDisplayImage() {
				Empty();
			}
			PDisplayImage(const char *fname) {
				Empty();
				SetImage(fname);
			}
			PDisplayImage(const DBRecord &pr) {
				Empty();
				SetImage(pr);
			}
	virtual		~PDisplayImage() {
				CloseImage();
			}
			// Assign the display color space
	bool		SetDisplaySpace(const ImgColorSpace *dcs,
					double lddyn = .0, double ldmax = .0);
			// Get the current display color space
	void		GetDisplaySpace(ImgColorSpace *csret,
					double *lddynp = NULL,
					double *ldmaxp = NULL) const {
				if (csret != NULL) PcopyCS(csret, &cs);
				if (lddynp != NULL) *lddynp = dynrange;
				if (ldmaxp != NULL) *ldmaxp = dlmax;
			}
			// Close the current image and free resources
	bool		CloseImage(bool cancancel = false);
			// Use the given image file path
	bool		SetImage(const char *fpath);
			// Use the given image record
	bool		SetImage(const DBRecord &pr);
			// Are we good to go?
	bool		Ready() const {
				return ((oxres > 0) & (oyres > 0));
			}
			// Determine the number of frames (0 if unknown)
	int		GetNFrames() {
				if (!CheckLoad()) return -1;
				return nframes;
			}
			// Seek to the given frame in a sequence
	bool		SeekFrame(int offs, ImgSeekMode sm = IRSabs);
			// Say which frame we're on
	int		TellFrame() const {
				if (dlp == NULL) return -1;
				return dlp->frame;
			}
			// Get frame type (and rate)
	ImgFrameType	GetFrameType(double *fr = NULL) {
				if (fr != NULL) *fr = 0;
				if (!CheckLoad()) return IRFnone;
				if (fr != NULL) *fr = frameRate;
				return frameType;
			}
			// Reload our image
	void		Reload() {
				if (!Ready()) return;
				IRclose(ir); ir = NULL;
				need2reload.ClearBitMap(true);
			}
			// Check that the current image is up to date
	bool		CheckLoad() {
				if (dlp != NULL) return dlp->CheckLoad();
				return SetFrame(0);
			}
			// Control current image's purgeability (pair calls!)
	void		SetPurgeable(bool purgeOK) {
				if (dlp == NULL) return;
				if (!purgeOK) dlp->Retain();
				else if (dlp->RetainCount() > 0) dlp->Release();
			}
			// Get the current image record (type PDBFInfo++)
	const DBRecord &
			GetRec() const {
				return ircd;
			}
			// Update the image record
	bool		UpdateRec(const DBRecord &pr);
			// Set the tone-mapping flags indicating what to do
	void		SetTMFlags(int newflags) {
				newTMflags = newflags;
			}
			// Get the present tone-mapping flag settings
	int		GetTMFlags() const {
				return newTMflags;
			}
			// Is the image in high dynamic-range?
	bool		IsHDR() {
				if (!CheckLoad()) return false;
				return (dlp->tms != NULL);
			}
			// Is the conversion to RGB lossy?
	bool		IsLossy() {
				if (!CheckLoad()) return false;
				return dlp->lossy;
			}
			// Get sample to nits calibration factor
	float		GetSF() const {
				if (dlp != NULL && dlp->tms != NULL)
					return dlp->tms->inpsf;
				return ircd[PDBFstonits].GetFloat();
			}
			// Get original image width (before crop & orient)
	int		GetOWidth() const {
				if (dlp != NULL) return dlp->oxres;
				return oxres;
			}
			// Get original image height
	int		GetOHeight() const {
				if (dlp != NULL) return dlp->oyres;
				return oyres;
			}
			// Get crop rectangle, return true if not whole image
	bool		CropRect(ImgRect *crp = NULL) const;
			// Convert viewable rectangle to original coordinates
	bool		GetOrigRect(ImgRect *orp, const ImgRect *r = NULL);
			// Convert original coordinates to viewable rectangle
	bool		GetViewRect(ImgRect *vrp, const ImgRect *r = NULL);
			// Correct coordinates to square aspect pixels
	void		CorrectAspect(int *xp, int *yp) const {
				if (dlp != NULL)
					dlp->CorrectAspect(xp, yp);
				else if (pixAspect < 0.995f)
					*xp = int(*xp*pixAspect + .5f);
				else if (pixAspect > 1.005f)
					*yp = int(*yp/pixAspect + .5f);
			}
			// Corresponding position in original image
	void		OrigAspect(int *xp, int *yp) const {
				if (dlp != NULL) {
					dlp->OrigAspect(xp, yp);
				} else if (pixAspect < 0.995f) {
					*xp = int(*xp/pixAspect + .5f);
					if (*xp > oxres) *xp = oxres;
				} else if (pixAspect > 1.005f) {
					*yp = int(*yp*pixAspect + .5f);
					if (*yp > oyres) *yp = oyres;
				}
			}
			// Disable false color mapping
	void		FalseColorDone() {
				if (fcs == NULL) return;
				fcDone(fcs); fcs = NULL; fcManual = -1;
			}
			// Set fixed linear mapping using the given multiplier
	bool		SetFixedLinear(double expmult = 1.) {
				if (!CheckLoad()) return false;
				if (fcs != NULL) {
					fcFixedLinear(fcs, GetSF()/expmult);
					fcManual = true;
				}
				return dlp->SetFixedLinear(expmult);
			}
			// Get the current multiplier (approx. if non-linear)
	double		GetMultiplier() {
				if (dlp == NULL) return .0;
				return dlp->GetMultiplier();
			}
			// These calls work in concert with PfreeImage()
	bool		GetSubimage(ImgStruct *ims, const ImgRect *r,
					double *localexp = NULL) {
				if (!CheckLoad()) return false;
				return dlp->GetSubimage(ims, r, localexp);
			}
	bool		GetImage(ImgStruct *ims) {
				return GetSubimage(ims, NULL);
			}
			// Average luminance over the given rectangle
	double		AvgLuminance(const ImgRect *r) {
				if (!CheckLoad()) return .0;
				return dlp->AvgLuminance(r);
			}
			// Compute histogram over the given rectangle
	bool		GetHistogram(double minmax[2],
					unsigned long *hist, int hlen,
					PHistoChan hc = PHCluminance,
					const ImgRect *r = NULL) {
				if (!CheckLoad()) return false;
				return dlp->GetHistogram(minmax, hist, hlen, hc, r);
			}
};

// Initialize a display list node
inline
PDisplayList::PDisplayList(PDisplayImage *par, int fno)
{
	Empty();
	parent = par;
	frame = fno;
	if (par != NULL)
		c_img.csp = &par->cs;
}

// Callback to delete original image on closing
bool		PDeleteTempImage(PDisplayImage *di, bool cancancel);

// Exposure mode flags (modifies false color behavior as well)
enum PExposMode {
		PEMauto = 0x1,			// non-linear mode
		PEMlocal = 0x2,			// local adjustment
		PEMhuman = 0x4,			// human sensitivity & color
		PEMmask = 0x7			// all modes mask
};

#define PEMisManual(pem)	(!((pem) & (PEMauto|PEMlocal)))

inline int
iround(double x)
{
	return int(x + .5) - int(x < -.5);
}

extern int	PEMDefault;			// default exposure mode

// Class to form a ring of display images (knows how to reorient)
class PImageRing : public PDisplayImage {
	int		froffs;			// offset to new frame
	ImgSeekMode	frsmode;		// frame seek mode
	int		zoomNum, zoomDenom;	// assigned zoom factor = Num/Denom
	int		curZNum, curZDenom;	// last rendered zoom factor
	int		xCent, yCent;		// current offset from image center
	int		expMode;		// current exposure mode flags
	double		expMult;		// current exposure multiplier
	ImgRect		srcrect;		// box in source coords
	PImageRing *	prev;			// previous image in ring
	PImageRing *	next;			// next image in ring
	void		SetDefaults() {
				froffs = 0; frsmode = IRSrel;
				zoomNum = zoomDenom = 0;
				curZNum = curZDenom = 1;
				xCent = yCent = 0;
				expMode = PEMDefault; expMult = 1.;
				ClearSelection();
				rWidth = rHeight = 0;
				changed = true;
			}
	void		MapToSource(int *sxp, int *syp, int dx, int dy) const {
				*sxp = ior.xres/2 + xCent +
					(2*dx - rWidth)*(long)curZDenom/curZNum/2;
				*syp = ior.yres/2 + yCent +
					(2*dy - rHeight)*(long)curZDenom/curZNum/2;
			};
	void		MapToDest(int *dxp, int *dyp, int sx, int sy) const {
				*dxp = rWidth/2 +
				    (sx - ior.xres/2 - xCent)*(long)curZNum/curZDenom;
				*dyp = rHeight/2 +
				    (sy - ior.yres/2 - yCent)*(long)curZNum/curZDenom;
			}
	void		CheckSettings();
public:
	bool		changed;		// display image needs update?
	int		rWidth, rHeight;	// rendered image size
			PImageRing(const char *fn = NULL, PImageRing *prv = NULL) {
				next = prev = this;
				fcs = NULL;
				JoinRing(prv);
				SetImage(fn);
			}
			PImageRing(const DBRecord &dr, PImageRing *prv = NULL) {
				next = prev = this;
				fcs = NULL;
				JoinRing(prv);
				SetImage(dr);
			}
			~PImageRing() {
				Clean();
			}
			// Clear the entire ring of images
	void		Clean() {
				prev->next = NULL;
				if (next != NULL)
					delete next;
				prev = next = this;
				CloseImage();
			}
			// Join with another ring after the given point
	void		JoinRing(PImageRing *prv) {
				if ((prv == NULL) | (prv == this)) return;
				prv->next->prev = prev;
				prev->next = prv->next;
				prev = prv;
				prev->next = this;
			}
			// Unlink this object from its ring and return previous
	PImageRing *	Unlink() {
				if ((this == NULL) | (this == next)) return NULL;
				PImageRing *	prv = prev;
				next->prev = prev;
				prev->next = next;
				next = prev = this;
				return prv;
			}
			// Set this image to the given file
	bool		SetImage(const char *fn) {
				CloseImage(); SetDefaults();
				if (fn == NULL || !*fn) return false;
				return PDisplayImage::SetImage(fn);
			}
			// Set this image to the given record
	bool		SetImage(const DBRecord &dr) {
				CloseImage(); SetDefaults();
				return PDisplayImage::SetImage(dr);
			}
			// Are we ready to rock?
	bool		Ready() const {
				if (this == NULL) return false;
				return PDisplayImage::Ready();
			}
			// Reload the current image
	bool		Reload() {
				PDisplayImage::Reload();
				return changed = true;
			}
			// Are horizontal and vertical swapped during mapping?
	bool		XYswapped() const {
				return ior.xyswap;
			}
			// Update this image record (i.e., to crop or reorient)
	bool		UpdateRec(const DBRecord &dr) {
				if (PDBgetField(dr,PDBFcrop) != GetRec()[PDBFcrop])
					xCent = yCent = 0;
				if (PDisplayImage::UpdateRec(dr))
					return changed = true;
				return false;
			}
			// Get the next image in our ring
	PImageRing *	GetNext() const {
				if (this == NULL) return NULL;
				return next;
			}
			// Get the previous image in our ring
	PImageRing *	GetPrev() const {
				if (this == NULL) return NULL;
				return prev;
			}
			// Get the cropped, zoomed size and orientation
	bool		GetOrient(POrient *porp) {
				if (this == NULL) return false;
				if (porp == NULL) return false;
				*porp = ior;
				porp->SetOrigSize(ior.xres*curZNum/curZDenom,
							ior.yres*curZNum/curZDenom);
				return true;
			}
			// Get the currect selection in mapped coordinates
	bool		GetSelection(ImgRect *rp) const;
			// Get selection in original image coordinates
	bool		GetOrigSelection(ImgRect *orp) {
				if ((this == 0) | (orp == 0)) return false;
				return GetOrigRect(orp, &srcrect);
			}
			// Set selection rectangle using original image coordinates
	void		SetOrigSelection(ImgRect *orp);
			// Set selection rectangle based on mapped point pair
	void		SelectRect(int x0, int y0, int x1, int y1);
			// Select the entire viewable source image
	void		SelectAll() {
				srcrect.xleft = srcrect.ytop = 0;
				srcrect.xright = ior.xres;
				srcrect.ybottom = ior.yres;
			}
			// return 0 if none selected, 1 if some selected, -1 if all
	int		SomeSelected() const {
				int	w = srcrect.xright - srcrect.xleft;
				if (w <= 0) return 0;
				int	h = srcrect.ybottom - srcrect.ytop;
				if (h <= 0) return 0;
				if ((w < ior.xres) | (h < ior.yres)) return 1;
				return -1;
			}
			// Clear the selection rectangle
	void		ClearSelection() {
				srcrect.xleft = srcrect.xright = 0;
				srcrect.ytop = srcrect.ybottom = 0;
			}
			// Seek to the given frame in a sequence
	bool		SeekFrame(int offs, ImgSeekMode sm = IRSabs) {
				if ((offs == 0) & (sm != IRSabs)) return false;
				if ((frsmode = sm) == IRSabs) froffs = offs;
				else froffs += offs;
				return changed = true;
			}
			// Say which frame we're supposedly on
	int		TellFrame() {
				int	frame = PDisplayImage::TellFrame();
				if (frame < 0) return -1;
				if (!PabsFrame(&frame, GetNFrames(), GetFrameType(),
						froffs, frsmode)) return -1;
				return frame;
			}
			// Get the current scroll offset in mapped pixels
	void		GetOffset(int hvoff[2]) const {
				if (ior.xyswap) {
					hvoff[0] = yCent * curZNum / curZDenom;
					if (ior.vflip) hvoff[0] = -hvoff[0];
					hvoff[1] = xCent * curZNum / curZDenom;
					if (ior.hflip) hvoff[1] = -hvoff[1];
				} else {
					hvoff[0] = xCent * curZNum / curZDenom;
					if (ior.hflip) hvoff[0] = -hvoff[0];
					hvoff[1] = yCent * curZNum / curZDenom;
					if (ior.vflip) hvoff[1] = -hvoff[1];
				}
			}
			// Set the current scroll offset
	void		SetOffset(int hoff, int voff) {
				if (IsZoomToFit()) return;
				int	xoff, yoff;
				if (XYswapped()) {
					xoff = voff * zoomDenom / zoomNum;
					yoff = hoff * zoomDenom / zoomNum;
				} else {
					xoff = hoff * zoomDenom / zoomNum;
					yoff = voff * zoomDenom / zoomNum;
				}
				if (ior.hflip) xoff = -xoff;
				if (ior.vflip) yoff = -yoff;
				if ((xoff == xCent) & (yoff == yCent)) return;
				xCent = xoff; yCent = yoff;
				changed = true;
			}
			// Get the current zoom ratio
	double		GetZoom(int *numden = NULL) const {
				if (numden != NULL) {
					numden[0] = curZNum;
					numden[1] = curZDenom;
				}
				return double(curZNum)/double(curZDenom);
			}
			// Set the desired zoom ratio
	void		SetZoomRatio(int num, int denom) {
				if ((num <= 0) | (denom <= 0)) return;
				if ((num == zoomNum) & (denom == zoomDenom)) return;
				curZNum = zoomNum = num;
				curZDenom = zoomDenom = denom;
				changed = true;
			}
			// Return true if zoom-to-fit is set
	bool		IsZoomToFit() const {
				return ((zoomDenom <= 0) | (zoomNum <= 0));
			}
			// Set zoom-to-fit mode (centers image)
	void		SetZoomToFit() {
				if (IsZoomToFit()) return;
				zoomNum = zoomDenom = 0;
				curZNum = curZDenom = 1;
				xCent = yCent = 0;
				changed = true;
			}
			// Get the current exposure mode (and multiplier)
	int		GetExposure(double *multp = NULL) const {
				if (multp != NULL) *multp = expMult;
				return expMode;
			}
			// Set the exposure mode
	void		SetExposureMode(int newmode) {
				newmode &= PEMmask;
				if (newmode == expMode) return;
				expMode = newmode;
				fcManual = -1;		// triggers update
				changed = true;
			}
			// Set particular exposure mode flag(s) on or off
	void		SetExposureMode(int modefl, bool switchon) {
				SetExposureMode( switchon ? (expMode | modefl)
							: (expMode & ~modefl) );
			}
			// Set the exposure multiplier (turns off auto exposure)
	void		SetExposure(double mult) {
				if (PEMisManual(expMode) &&
					(expMult*.95 <= mult) &
					(mult <= expMult*1.05)) return;
				expMode &= ~(PEMauto|PEMlocal);
				expMult = mult;
				fcManual = -1;		// triggers update
				changed = true;
			}
			// Enable/change false color mapping for HDR image
	void		FalseColorOn(double maxv = .0, double minv = .0);
			// Get false color limits, if set
	bool		GetFCLimits(double minmax[2]) const;
			// Disable false color mapping
	void		FalseColorOff() {
				if (fcs == NULL) return;
				FalseColorDone();
				changed = true;
			}
			// Return true if no crop, resize (or tonemapping)
	bool		RenderIsNoop(bool ignoreHDR = false) {
				if (!Ready()) return true;
				if (IsZoomToFit()) return false;
				if (zoomNum != zoomDenom) return false;
				if (CropRect()) return false;
				if (GetRec().GetField(PDBFredeye) != NULL) return false;
				if (IsHDR()) { if (!ignoreHDR) return false; }
				else if (GetRec().GetField(PDBFspotexp) != NULL) return false;
				if ((ior.GetWidth() != GetOWidth()) |
						(ior.GetHeight() != GetOHeight()))
					return false;
				return true;
			}
			// Return true if no orient, crop, resize (or tonemapping)
	bool		MapIsNoop(bool ignoreHDR = false) {
				if (!Ready()) return true;
				if (!ior.IsNoop()) return false;
				return RenderIsNoop(ignoreHDR);
			}
			// Render the image into the given window (orig. orientation)
	bool		RenderImage(ImgStruct *ims, const void *pf = NULL,
						const ImgRect *sir = NULL);
			// Render to subsection of destination (orig. orientation)
	bool		RenderSubimage(ImgStruct *ims, const ImgRect *sir,
						const void *pf = NULL) {
				return RenderImage(ims, pf, sir);
			}
			// Render the image into the given window (reorients)
	bool		MapImage(ImgStruct *ims, const void *pf = NULL,
						const ImgRect *sir = NULL);
			// Render to subsection of destination (reorients)
	bool		MapSubimage(ImgStruct *ims, const ImgRect *sir,
						const void *pf = NULL) {
				return MapImage(ims, pf, sir);
			}
			// Average luminance over the selected region
	double		AvgLuminance(const ImgRect *rct = NULL) {
				if (rct == NULL) rct = &srcrect;
				return PDisplayImage::AvgLuminance(rct);
			}
			// Compute histogram over image or selected region
	bool		GetHistogram(double minmax[2],
					unsigned long *hist, int hlen,
					PHistoChan hc = PHCluminance,
					bool useSel = true) {
				const ImgRect * rct = NULL;
				if (useSel && SomeSelected() > 0)
						rct = &srcrect;
				return PDisplayImage::GetHistogram(minmax,
						hist, hlen, hc, rct);
			}
};

#endif	// ! _PDISPIMG_H_
