/*
 *  panimage.h
 *  panlib
 *
 *  Basic Pancine image class
 *
 *  Created by Greg Ward on 6/24/06.
 *  Copyright 2023 Anyhere Software. All rights reserved.
 *
 */

#ifndef _PANIMAGE_H_
#define _PANIMAGE_H_

#include "pimage.h"
#include "imgwriter.h"

/** @brief Wrapper class for Pancine image type
 *
 * Class to provide automatic allocation and deallocation,
 * reading and writing of Pancine images with basic accessors.
 * Simple operators are provided for floating-point images.
 * Other operations are left to C routines using ImgStruct.
 */

#define PAN_MAXWRITER	15		// limit on image writer interfaces

/// Image state flags:  uninitialized, color space only, unallocated, readable, r/w, clear
enum PImgCond { PICvoid=0, PICspace=1, PICfree=3, PICread=4, PICready=6, PICzero=8 };

/// Math flag settings: strict color typing, report warnings as errors
enum { PIMstrictcolor=1, PIMwarn2error=2 };

/// Image rotation & flip aliases for Reorient() call
enum { PIOcw90=IOrighttop, PIO180=IObotright, PIOccw90=IOleftbot,
	PIOhflip=IOtopright, PIOvflip=IObotleft,
	PIOcw90hflip=IOlefttop, PIOcw90vflip=IOrightbot };

/// Pancine image class
class PanImage {
private:
	ImgColorSpace *	sCSp;		// saved color space
	ImgInfo *	iinfo;		// image metadata
protected:
	static const ImgWriterInterface *
			imgWriter[PAN_MAXWRITER+1];
	ImgStruct	ims;		// image size, pointers, etc.
	short		ro;		// read-only (-1 not alterable)
			/// Zero out pointers
	void		Clear() {
				ims.xres = ims.yres = 0; ims.rowsize = 0;
				ims.csp = NULL;
				ims.mbase = NULL; ims.img = NULL;
				sCSp = NULL;
				iinfo = NULL;
				ro = 0;
				mathFlags = defMathFlags;
			}
			/// Protect color space against PfreeImage()
	bool		PreserveCS(const ImgColorSpace *csp = NULL);
			/// Check operands for valid operation
	bool		OperatorChecks(const PanImage &src, int lineno) const;
public:
			/// default math operation flags
	static short	defMathFlags;
			/// math operation flags for this object
	short		mathFlags;
			PanImage(ImgStruct *im = NULL) {
				Clear();
				if (im) Link(im);
			}
			PanImage(const PanImage &orig) {
				Clear();
				*this = orig;
			}
			PanImage(int xr, int yr, const ImgColorSpace &imgCS,
						PImgCond act = PICfree) {
				Clear();
				Init(xr, yr, imgCS, act);
			}
			PanImage(const ImgColorSpace &imgCS) {
				Clear();
				Init(imgCS);
			}
			PanImage(const char *imgFile, bool qt = true) {
				Clear();
				Load(imgFile, qt);
			}
			~PanImage() {
				PfreeImage(&ims);
				delete sCSp;
				delete iinfo;
			}
			/// Initialize image
	bool		Init(int xr, int yr, const ImgColorSpace &imgCS,
						PImgCond act = PICready);
			/// Initialize with color space only
	bool		Init(const ImgColorSpace &imgCS) {
				return Init(0, 0, imgCS, PICspace);
			}
			/// Check/change initialization
	bool		Init(PImgCond act = PICvoid);
			/// Determine image state (doesn't detect PICzero)
	PImgCond	Ready() const {
				if (!ims.csp) return PICvoid;
				if ((ims.xres <= 0) | (ims.yres <= 0)) return PICspace;
				if (!ims.img) return PICfree;
				return ro ? PICread : PICready;
			}
			/// Is image currently read-only?
	bool		ReadOnly() const {
				return (ims.img != NULL) & (ro != 0);
			}
			/// Set/clear read-only property
	bool		SetReadOnly(bool set2 = true) {
				if (!ims.img) return !set2;
				if (ro < 0 && (!ims.mbase ||
						ims.mbase->RetainCount() > 1))
					return set2;
				ro = set2; return true;
			}
			/// How many links to image?
	int		NLinks() const {
				if (!ims.img) return 0;
				if (!ims.mbase) return -1;
				return ims.mbase->RetainCount();
			}
			/// Query image width
	int		Width() const {
				return ims.xres;
			}
			/// Query image height
	int		Height() const {
				return ims.yres;
			}
			/// Query number of pixel components
	int		NComp() const {
				if (!ims.csp) return 0;
				return ImgPixelLen[ims.csp->format];
			}
			/// Image pixel size in bytes
	int		PSize() const {
				if (!ims.csp) return 0;
				return ImgPixelSize(ims.csp);
			}
			/// Check if image colorspace (& size) are compatible
	bool		Compat(const ImgStruct *im, int cfl = PICMequiv,
					bool matchSize = true) const {
				if (!im) return false;
				if (matchSize) {
					if (!ims.img | !im->img)
						return false;
					if ((ims.xres != im->xres) | (ims.yres != im->yres))
						return false;
				}
				return PmatchColorSpace(ims.csp, im->csp, cfl);
			}
	bool		Compat(const PanImage &other, int cfl = PICMequiv,
					bool matchSize = true) const {
				return Compat(&other.ims, cfl, matchSize);
			}
			/// Take ownership of existing image
	bool		Take(ImgStruct *im) {
				if (im == &ims) return false;
				Init();
				if (!im || !im->img) return false;
				ims = *im; im->img = NULL;
				if (ims.csp == MobjMem(ims.mbase)) im->csp = NULL;
				else PreserveCS();
				return true;
			}
			/// Take ownership, which may clear source entirely
	bool		Take(PanImage *srcp) {
				if (!srcp) return false;
				if (!Take(&srcp->ims)) return false;
				iinfo = srcp->iinfo; srcp->iinfo = NULL;
				ro = srcp->ro; mathFlags = srcp->mathFlags;
				return true;
			}
			/// Load image from file or object
			/// Use current size & color space if initialized
	bool		Load(ImgReader *ir, bool quiet = false);
	bool		Load(const char *imgFile, bool quiet = false);
	bool		Load(const PanImage &src);
			/// Register image writer
	static bool	AddImageWriter(const ImgWriterInterface *iwriter);
			/// Write image to file (suffix determines type)
	long		Write(const char *imgFile, int qual = -1) const;
			/// Get pointer to required image metadata
	const ImgInfo *	GetInfo(long flgs = 0) const {
				if (!iinfo || !iinfo->flags) return NULL;
				if ((iinfo->flags & flgs) != flgs) return NULL;
				return iinfo;
			}
			/// Copy image metadata, returning what's set
	long		GetInfo(ImgInfo *infp) const {
				if (!infp) return 0;
				if (!iinfo) *infp = defImgInfo;
				else if (infp != iinfo) *infp = *iinfo;
				return infp->flags;
			}
			/// Clear metadata
	void		ClearInfo(long flgs = ~0L) {
				if (!iinfo) return;
				if (!(iinfo->flags &= ~flgs)) {
					delete iinfo; iinfo = NULL;
					return;
				}
				if (flgs & IIFview) iinfo->view[0] = '\0';
				if (flgs & IIFparams) iinfo->params[0] = '\0';
				if (flgs & IIFcomments) iinfo->comments[0] = '\0';
			}
			/// Access image metadata to modify
	ImgInfo *	Info(long flg = 0) {
				if (!iinfo)
					*(iinfo = new ImgInfo) = defImgInfo;
				iinfo->flags |= flg;
				return iinfo;
			}
			/// Reorient image from "orienta" (or metadata if 0)
	bool		Reorient(int orienta = 0, bool force = false);
			/// Get pointer to read-only image struct
	const ImgStruct *
			GetImg() const {
				if (!ims.img) return NULL;
				return &ims;
			}
			/// Access image struct to modify
	ImgStruct *	Img() {
				if (!ims.img) ro = 0;
				else if (ro) return NULL;
				return &ims;
			}
			/// Get image color space
	const ImgColorSpace *
			GetCS() const {
				return ims.csp;
			}
	bool		GetCS(ImgColorSpace *cs) const {
				if (!cs | !ims.csp) return false;
				PcopyCS(cs, ims.csp);
				return true;
			}
			/// Get a single pixel value
	PixelVal	GetPixel(int x, int y) const {
				int		n = PSize();
				PixelVal	p;
				const uby8 *	bp;
				if (!n || (x < 0) | (x >= Width()) ||
						!(bp = GetRow(y)))
					return Pblack;
				p.csp = ims.csp;
				bp += (x+1)*n;
				while (n--) p.v.b[n] = *--bp;
				return p;
			}
			/// Assign a single pixel value (may be slow)
	bool		SetPixel(int x, int y, PixelVal p) {
				int	n = PSize();
				uby8 *	bp;
				if (!n || (x < 0) | (x >= Width()) ||
						!(bp = Row(y)))
					return false;
				if (p.csp != ims.csp)
					p = PconvertPixel(p, ims.csp);
				bp += (x+1)*n;
				while (n--) *--bp = p.v.b[n];
				return true;
			}
			/// Draw pixel-wide line between image points
	int		DrawLine(int x0, int y0, int x1, int y1, PixelVal p);
			/// Get pixel row pointer (no typing)
	const uby8 *	GetRow(int y) const {
				if (!ims.img) return NULL;
				if ((y < 0) | (y >= Height())) return NULL;
				return ProwPtr(&ims, y);
			}
			/// Array operator gets untyped row pointer w/o checks
	const uby8 *	operator[](int y) const {
				return ProwPtr(&ims, y);
			}
			/// Get pointer to untyped pixel
	const uby8 *	GetP(int x, int y) const {
				if ((x < 0) | (x >= Width())) return NULL;
				const uby8 *	pv = GetRow(y);
				if (!pv) return NULL;
				return pv + x*PSize();
			}
			/// Get 8-bit pixel row
	const uby8 *	GetRowByte(int y) const {
				if (!ims.csp) return NULL;
				if (ims.csp->dtype != IDTubyte) return NULL;
				return GetRow(y);
			}
			/// Get 16-bit pixel row
	const unsigned short *
			GetRowShort(int y) const {
				if (!ims.csp) return NULL;
				if (ims.csp->dtype != IDTushort) return NULL;
				return (const unsigned short *)GetRow(y);
			}
			/// Get floating-point pixel row
	const float *	GetRowFloat(int y) const {
				if (!ims.csp) return NULL;
				if (ims.csp->dtype != IDTfloat) return NULL;
				return (const float *)GetRow(y);
			}
			/// Pixel row pointer for writing (no typing)
	uby8 *		Row(int y) {
				if (!ims.img | ro) return NULL;
				if ((y < 0) | (y >= Height())) return NULL;
				return ProwPtr(&ims, y);
			}
			/// Array operator gets untyped row pointer w/o checks
	uby8 *		operator[](int y) {
				// if (ro) return NULL;	XXX not even this
				return ProwPtr(&ims, y);
			}
			/// Get writeable pointer to untyped pixel
	uby8 *		P(int x, int y) {
				if ((x < 0) | (x >= Width())) return NULL;
				uby8 *	pv = Row(y);
				if (!pv) return NULL;
				return pv + x*PSize();
			}
			/// 8-bit pixel row for writing
	uby8 *		RowByte(int y) {
				if (!ims.csp) return NULL;
				if (ims.csp->dtype != IDTubyte) return NULL;
				return Row(y);
			}
			/// 16-bit pixel row for writing
	unsigned short *
			RowShort(int y) {
				if (!ims.csp) return NULL;
				if (ims.csp->dtype != IDTushort) return NULL;
				return (unsigned short *)Row(y);
			}
			/// Floating-point pixel row for writing
	float *		RowFloat(int y) {
				if (!ims.csp) return NULL;
				if (ims.csp->dtype != IDTfloat) return NULL;
				return (float *)Row(y);
			}
			/// In situ color conversion (may delink image)
	bool		ConvertCS(const ImgColorSpace &imgCS,
						float sf = 1.f);
			/// Expert-only method to munge color space
	bool		AssertCS(const ImgColorSpace &imgCS) {
				return PreserveCS(&imgCS);
			}
			/// Link writeable (sub)image from ImgStruct
	bool		Link(ImgStruct *im, const ImgRect *r = NULL);
			/// Link (sub)image from PanImage object
	bool		Link(PanImage *srcp, const ImgRect *r = NULL);
			/// Link (sub)image from read-only ImgStruct
	bool		LinkRO(const ImgStruct *im, const ImgRect *r = NULL);
			/// Link (sub)image from const PanImage object
	bool		LinkRO(const PanImage &src, const ImgRect *r = NULL);
			/// Make private copy of our image (clears ReadOnly)
	bool		Delink() {
				if (!SetReadOnly(false)) {
					PanImage	imcpy(*this);
					return Take(&imcpy);
				}
				return PdelinkImage(&ims);
			}
			/// Crop image using given rectangle (or metadata)
	bool		Crop(const ImgRect *r = NULL) {
				if (!ims.img) return false;
				if (!r) {
					if (!GetInfo(IIFcrop)) return false;
					r = &iinfo->crop;
				}
				if (!PlegalRect(r, ims.xres, ims.yres))
					return false;
				return Link(&ims, r);
			}
			/// Crop image using given corners (more forgiving)
	bool		Crop(int x1, int y1, int x2, int y2);
			/// Access writeable image section
	PanImage	operator[](const ImgRect &rct) {
				PanImage	im;
				im.Link(this, &rct);
				return im;
			}
			/// Access read-only image section
	const PanImage	operator[](const ImgRect &rct) const {
				PanImage	im;
				im.LinkRO(*this, &rct);
				return im;
			}
			/// Compute mean pixel value (in the given CS)
	PixelVal	GetAverage(const ImgColorSpace *csp = NULL) const {
				return PcomputeAverage(&ims, csp);
			}
			/// Alpha blending (Load() semantics but no resize)
	bool		Blend(const PanImage &a0, const PanImage &alpha,
					const PanImage &a1);
			/// Alpha blending, current image as background
	bool		BlendIn(const PanImage &alpha, const PanImage &a1) {
				return PblendImages(Img(), GetImg(),
						alpha.GetImg(), a1.GetImg());
			}
			/// Operation blending with alpha mask
	bool		BlendOp(const PanImage &alpha, PimageOp *op,
						void *udp = NULL) {
				return PblendImageOpCB(Img(), alpha.GetImg(),
							op, udp);
			}
			/// Image copy operator (clears read-only if not self)
	PanImage &	operator=(const PanImage &src);
			/// Scalar operators
	PanImage &	operator=(const PixelVal &pval);
	PanImage &	operator=(const float fval) {
				return *this = PgraY(fval);
			}
	PanImage &	operator+=(const PixelVal &pval);
	PanImage &	operator+=(const float fval) {
				return *this += PgraY(fval);
			}
	PanImage &	operator-=(const PixelVal &pval);
	PanImage &	operator-=(const float fval) {
				return *this += -fval;
			}
	PanImage &	operator*=(const PixelVal &pval);
	PanImage &	operator*=(const float fval);
	PanImage &	operator/=(const PixelVal &pval);
	PanImage &	operator/=(const float fval);
			/// Clamp to some minimum
	PanImage &	operator|=(const PixelVal &pval);
	PanImage &	operator|=(const float vmin) {
				return *this |= PgraY(vmin);
			}
			/// Clamp to some maximum
	PanImage &	operator&=(const PixelVal &pval);
	PanImage &	operator&=(const float vmax) {
				return *this &= PgraY(vmax);
			}
			/// Dilate image using disk of given radius
	PanImage &	operator<<=(const float rad);
			/// Erode image using disk of given radius
	PanImage &	operator>>=(const float rad) {
				return *this <<= -rad;
			}
			/// Gaussian blur (or sharpen if rad < 0.5)
	PanImage &	operator^=(const float rad);
			/// Image (pixel-by-pixel) operators, same semantics
	PanImage &	operator+=(const PanImage &src);
	PanImage &	operator-=(const PanImage &src);
	PanImage &	operator*=(const PanImage &src);
	PanImage &	operator/=(const PanImage &src);
	PanImage &	operator|=(const PanImage &src);
	PanImage &	operator&=(const PanImage &src);
			/// Apply the given convolution kernel
	PanImage &	operator^=(const PanImage &kern);
};

/// Choose destination color space based on two operands
extern const ImgColorSpace *	PickCS(const PanImage &im1, const PanImage &im2,
					bool makeReal = false);

// Alpha blending (Load() semantics but no resize)
inline bool
PanImage::Blend(const PanImage &a0, const PanImage &alpha, const PanImage &a1)
{
	if (!ims.csp) Init(*PickCS(a0, a1));
	return PblendImages(Img(), a0.GetImg(), alpha.GetImg(), a1.GetImg());
}

/// Arithmetic operators hoping on return value optimization
inline PanImage
operator+(const PanImage &imLeft, const PanImage &imRight)
{
	PanImage	imResult(*PickCS(imLeft, imRight));
	if (imResult.Load(imLeft)) imResult += imRight;
	return imResult;
}

inline PanImage
operator-(const PanImage &imLeft, const PanImage &imRight)
{
	PanImage	imResult(*PickCS(imLeft, imRight));
	if (imResult.Load(imLeft)) imResult -= imRight;
	return imResult;
}

inline PanImage
operator*(const PanImage &imLeft, const PanImage &imRight)
{
	PanImage	imResult(*PickCS(imLeft, imRight, true));
	if (imResult.Load(imLeft)) imResult *= imRight;
	return imResult;
}

inline PanImage
operator/(const PanImage &imLeft, const PanImage &imRight)
{
	PanImage	imResult(*PickCS(imLeft, imRight, true));
	if (imResult.Load(imLeft)) imResult /= imRight;
	return imResult;
}

inline PanImage
operator|(const PanImage &imLeft, const PanImage &imRight)
{
	PanImage	imResult(*PickCS(imLeft, imRight));
	if (imResult.Load(imLeft)) imResult |= imRight;
	return imResult;
}

inline PanImage
operator&(const PanImage &imLeft, const PanImage &imRight)
{
	PanImage	imResult(*PickCS(imLeft, imRight));
	if (imResult.Load(imLeft)) imResult &= imRight;
	return imResult;
}

inline PanImage
operator^(PanImage imLeft, const PanImage &imRight)
{
	return imLeft ^= imRight;
}

inline PanImage
operator^(PanImage imLeft, float rad)
{
	return imLeft ^= rad;
}

/// Compare pixel values for equality
inline bool
operator==(const PixelVal &p1, const PixelVal &p2)
{
	return PequalVal(p1, p2);
}

inline bool
operator!=(const PixelVal &p1, const PixelVal &p2)
{
	return !PequalVal(p1, p2);
}

/// Compare pixel Y values for inequality
inline bool
operator>(const PixelVal &p1, const PixelVal &p2)
{
	return (PgetY(p1) > PgetY(p2));
}

inline bool
operator<(const PixelVal &p1, const PixelVal &p2)
{
	return (PgetY(p1) < PgetY(p2));
}

inline bool
operator>=(const PixelVal &p1, const PixelVal &p2)
{
	return (PgetY(p1) >= PgetY(p2));
}

inline bool
operator<=(const PixelVal &p1, const PixelVal &p2)
{
	return (PgetY(p1) <= PgetY(p2));
}

// Pixel arithmetic ignores log/gamma, so is not always same as image ops

/// Add two pixels
extern PixelVal		operator+(PixelVal p1, PixelVal p2);

/// Subtract two pixels
extern PixelVal		operator-(PixelVal p1, PixelVal p2);

/// Multiply two pixels
extern PixelVal		operator*(PixelVal p1, PixelVal p2);

/// Multiply pixel by scalar
extern PixelVal		operator*(PixelVal p, float f);

/// Pixel negation (returns Pblack if not float)
inline PixelVal
operator-(PixelVal p)
{
	if (!p.csp || p.csp->dtype != IDTfloat) return Pblack;
	int i = ImgPixelLen[p.csp->format];
	while (i--) p.v.f[i] = -p.v.f[i];
	return p;
}

inline PixelVal
operator*(float f, const PixelVal &p)
{
	return p * f;
}

inline PixelVal
operator*=(PixelVal &p, float f)
{
	return p = p * f;
}

inline PixelVal
operator+=(PixelVal &p1, const PixelVal &p2)
{
	if (!p1.csp) return p1 = p2;
	return p1 = p1 + PconvertPixel(p2, p1.csp);
}

inline PixelVal
operator-=(PixelVal &p1, const PixelVal &p2)
{
	if (!p1.csp) return p1 = -p2;
	return p1 = p1 - PconvertPixel(p2, p1.csp);
}

/// Divide pixel by scalar
inline PixelVal
operator/(const PixelVal &p, float f)
{
	if (f == 0) return Pblack;
	return p * 1.f/f;
}

inline PixelVal
operator/=(PixelVal &p1, float f)
{
	return p1 = p1 / f;
}

/// Load standard image writers
extern bool		PloadStandardWriters();

/// Callback for quantization increment (max. difference to lesser or greater)
typedef int	PquantMethod(int x, int y, int c, void *udp);

/// Callback using integer image to look up quantization value
PquantMethod		PQmethodImage;

/// Simple call-back to return constant value from integer's address
PquantMethod		PQmethodConstant;

/// Dequantize image (using optional quantization from call-back at each pixel)
extern bool		Pdequantize(PanImage *dstp, const PanImage &isrc,
					float sf=1.f,
					PquantMethod *qm=NULL, void *udp=NULL);

/// Estimate quantization map for image based on neighbor values
extern bool		Pquantization(PanImage *dstp, const PanImage &isrc);

/// Combined call to Pquantization() and Pdequantize()
inline bool
Pdecontour(PanImage *dstp, const PanImage &isrc, float sf=1.f)
{
	PanImage	quantIm;

	return Pquantization(&quantIm, isrc) &&
		Pdequantize(dstp, isrc, sf, PQmethodImage, quantIm.Img());
}

#endif	// ! _PANIMAGE_H_
