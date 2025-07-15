/*
 *  panimage.cpp
 *  panlib
 *
 *  Basic Pancine image class implementation
 *
 *  Created by Greg Ward on 6/23/06.
 *  Copyright 2006 Anyhere Software. All rights reserved.
 *
 */

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "system.h"
#include "dmessage.h"
#include "panimage.h"
#include "pstrings.h"
#include "abitmap.h"

const ImgWriterInterface *	PanImage::imgWriter[PAN_MAXWRITER+1] = {NULL};
short				PanImage::defMathFlags = PIMstrictcolor;

#define	MATHWARN		(mathFlags&PIMwarn2error ? DMCparameter : DMCwarning)

// Initialize image
bool
PanImage::Init(int xr, int yr, const ImgColorSpace &imgCS, PImgCond act)
{
	if (act == PICvoid) {
		DMESG(DMCparameter, "Illegal call to PanImage::Init(...PICvoid)");
		return false;
	}
	if (&imgCS == ims.csp) {
		Init(PICspace);
	} else {
		Init(PICvoid);
		if (act & PICspace)
			PreserveCS(&imgCS);
		else
			ims.csp = &imgCS;
	}
	if (act == PICspace)
		return true;
	if ((xr <= 0) | (yr <= 0))
		return false;
	ims.xres = xr; ims.yres = yr;
	if (act == PICfree) {
		ims.rowsize = xr * ImgPixelSize(ims.csp);
		return true;
	}
	return Init(act);	// allocates (and optionally clears) image
}

// Check/change initialization
bool
PanImage::Init(PImgCond act)
{
	if (act == PICread) {
		DMESG(DMCparameter, "Illegal call to PanImage::Init(PICread)");
		act = PICvoid;	// probably not what they wanted
	}
	if (!SetReadOnly(false) || act <= PICfree) {
		if (act == PICvoid)
			delete sCSp;
		else
			PreserveCS();
		PfreeImage(&ims);
		delete iinfo;
		if (act == PICvoid) {
			Clear();
			return true;
		}
		iinfo = NULL;
		ro = 0;
		if (act == PICspace)
			ims.xres = ims.yres = 0;
		if (act & PICspace)
			return (Ready() == act);
	}
	if (!ims.img && !PnewImage(&ims, .0))
		return false;
	delete sCSp; sCSp = NULL;
	if (act == PICzero) {
		if (mathFlags & PIMstrictcolor && ims.csp->logorig > 0)
			DMESG(MATHWARN, "Setting log image to zeroes");
		PclearImage(&ims, NULL);
	}
	return true;
}

// Protect color space against call to PfreeImage()
bool
PanImage::PreserveCS(const ImgColorSpace *csp)
{
	if (!csp)
		csp = ims.csp;
	else if (ims.img && (csp->dtype != ims.csp->dtype) |
			(ImgPixelLen[csp->format] != ImgPixelLen[ims.csp->format]))
		return false;
	if (!csp || csp->cstatic) {
		delete sCSp; sCSp = NULL;
		return ((ims.csp = csp) != NULL);
	}
	if (csp == sCSp)
		return true;
	if (!sCSp)
		sCSp = new ImgColorSpace;
	PcopyCS(sCSp, csp);
	ims.csp = sCSp;
	return true;
}

// Load image from open reader (use current size & CS if initialized)
bool
PanImage::Load(ImgReader *ir, bool quiet)
{
	if (!ir)
		return false;
	switch (Ready()) {
	case PICread:
		DMESGF(DMCparameter, "Cannot load '%s' into locked image",
				ir->file);
		return false;
	case PICspace:					// our color space
		ims.xres = ir->xres;
		ims.yres = ir->yres;
		break;
	case PICvoid:					// initialize
		if (!Init(ir->xres, ir->yres, ir->cs, PICfree))
			return false;
		break;
	default:
		break;
	}
	ro = 0;
	int	len = sprintf(dmessage_buf, "Reading %s '%s'",
			ir->encoding ? ir->encoding : PdescribeCS(&ir->cs, NULL),
			ir->file);
	if (ims.csp && !PmatchColorSpace(ims.csp, &ir->cs, PICMall)) {
		if (PmatchColorSpace(ims.csp, &ir->cs, PICMptype|PICMgamma)) {
			strcpy(dmessage_buf+len, " and transforming colors");
		} else {
			strcpy(dmessage_buf+len, " as "); len += 4;
			PdescribeCS(ims.csp, dmessage_buf+len);
		}
	}
	DMESG(DMCtrace, dmessage_buf);
	if (!PrenderImageR(&ims, ir, quiet))		// load the image
		return false;
							// assign info.
	if (IRgetInfo(ir, Info()) != IREnone || !iinfo->flags) {
		ClearInfo();
		return true;
	}
	if ((ims.xres != ir->xres) | (ims.yres != ir->yres)) {
		const float	xsf = float(ims.xres)/float(ir->xres);
		const float	ysf = float(ims.yres)/float(ir->yres);
		if (iinfo->flags & IIFhdensity)
			iinfo->hdensity *= xsf;
		if (iinfo->flags & IIFcrop) {
			iinfo->crop.xleft = int(xsf*iinfo->crop.xleft + .5f);
			iinfo->crop.xright = int(xsf*iinfo->crop.xright + .5f);
			iinfo->crop.ytop = int(ysf*iinfo->crop.ytop + .5f);
			iinfo->crop.ybottom = int(ysf*iinfo->crop.ybottom + .5f);
		}
	}
	return true;
}

// Load image from file (use current size & CS if initialized)
bool
PanImage::Load(const char *imgFile, bool quiet)
{
	if (!imgFile || !*imgFile)
		return false;

	ImgReader *	ir = PopenImageF(imgFile, quiet, NULL);

	bool	ok = Load(ir, quiet);

	IRclose(ir);

	return ok;
}

// Load image from another image (use current size & CS if initialized)
bool
PanImage::Load(const PanImage &src)
{
	if (!src.ims.img)
		return false;
	if (&src == this)
		return true;
	switch (Ready()) {
	case PICread:
		DMESG(DMCparameter, "Cannot load image into locked destination");
		return false;
	case PICspace:					// our color space
		ims.xres = src.ims.xres;
		ims.yres = src.ims.yres;
		break;
	case PICvoid:					// fresh destination
		Init(src.ims.xres, src.ims.yres, *src.ims.csp, PICfree);
		break;
	default:
		break;
	}
	ro = 0;
	ClearInfo();					// transfer image
	if (!PrenderImageI(&ims, &src.ims))
		return false;
	if (!src.GetInfo())
		return true;
	src.GetInfo(Info());				// assign info.
	if (ims.xres != src.ims.xres) {
		const float	sf = float(ims.xres)/float(src.ims.xres);
		if (iinfo->flags & IIFhdensity)
			iinfo->hdensity *= sf;
		if (iinfo->flags & IIFcrop) {
			iinfo->crop.xleft = int(sf*iinfo->crop.xleft + .5f);
			iinfo->crop.xright = int(sf*iinfo->crop.xright + .5f);
			iinfo->crop.ytop = int(sf*iinfo->crop.ytop + .5f);
			iinfo->crop.ybottom = int(sf*iinfo->crop.ybottom + .5f);
		}
	}
	return true;
}

// In situ color conversion (may delink image)
bool
PanImage::ConvertCS(const ImgColorSpace &imgCS, float sf)
{
	if (!ims.img) {
		DTEST((0.999f > sf) | (sf > 1.001f), DMCwarning,
				"Cannot apply multiplier to unallocated image");
		return PreserveCS(&imgCS);
	}
	if (ro) {
		DMESG(DMCparameter, "ConvertCS() called on locked image");
		return false;
	}
	if ((0.999f <= sf) & (sf <= 1.001f) &&		// safe to coerce CS?
			PmatchColorSpace(&imgCS, ims.csp, PICMequiv))
		return PreserveCS(&imgCS);

	if ((ims.csp->dtype == IDTfloat) & (ims.csp->logorig <= 0) &&
			PmatchColorSpace(&imgCS, ims.csp, PICMall)) {
		if (NLinks() > 1 && !Delink())		// careful of linked image(s)
			return false;
		*this *= sf;				// multiply (adjusts stonits)
		return true;
	}
	if (!PconvertColorSpace(&ims, &imgCS, sf))	// else delink/convert pixels
		return false;
	if (GetInfo(IIFstonits))
		iinfo->stonits /= sf;			// & adjust stonits
	return true;
}

// Register image writer
bool
PanImage::AddImageWriter(const ImgWriterInterface *iwriter)
{
	int	i;
	for (i = 0; imgWriter[i]; i++)
		if (imgWriter[i] == iwriter)
			return true;		// already got this one
	if (i >= PAN_MAXWRITER) {
		DMESG(DMCparameter, "Too many image writers");
		return false;
	}
	imgWriter[i] = iwriter;
	imgWriter[i+1] = NULL;
	return true;
}

// Write image to file (suffix determines type)
long
PanImage::Write(const char *imgFile, int qual) const
{
	if (!imgFile || !*imgFile)
		return 0;
	if (!ims.img)
		return 0;
	const char *	sfx = PgetSuffix(imgFile);
	if (!sfx) {
		DMESGF(DMCparameter, "Missing suffix for image file '%s'", imgFile);
		return 0;
	}
	int	i;
	for (i = 0; imgWriter[i]; i++)
		if (!strcasecmp(sfx, imgWriter[i]->suffix))
			break;
	if (!imgWriter[i]) {
		DMESGF(DMCparameter, "Unknown image suffix '%s'", sfx);
		return 0;
	}
	ImgWriteBuf	iwb;
	if (!PsetWriteBuf(&iwb, &ims))
		return 0;
	GetInfo(&iwb.info);
	if ((0 <= qual) & (qual <= 100)) {
		iwb.info.quality = qual;
		iwb.info.flags |= IIFquality;
	} else
		qual = (iwb.info.flags & IIFquality) ? iwb.info.quality : -1;
	const char *	itype = (*imgWriter[i]->SupportedCS)(iwb.csp, qual);
	if (!itype) {
		sprintf(dmessage_buf, "%s not supported for %s output",
				PdescribeCS(iwb.csp,NULL), imgWriter[i]->suffix);
		DMESG(DMCparameter, dmessage_buf);
		return 0;
	}
	if (iwb.info.flags & IIFview &&	// check pixel aspect ratio
			!FindImgInfoParam(&iwb.info, "SAMP360")) {
		const char *	cp = strstr(iwb.info.view, "-vt");
		const int	vtype = cp ? cp[3] : 'v';
		double		hsiz=0, vsiz=0;
		cp = strstr(iwb.info.view, "-vh ");
		if (cp) hsiz = atof(cp+4);
		cp = strstr(iwb.info.view, "-vv ");
		if (cp) vsiz = atof(cp+4);
		if ((hsiz > .1) & (vsiz > .1)) {
			if (vtype == 'v')
				iwb.pixAspect = iwb.yres * tan(M_PI/360.*hsiz) /
					 	(iwb.xres * tan(M_PI/360.*vsiz));
			else if (vtype == 'l')
				iwb.pixAspect = iwb.yres*hsiz / (iwb.xres*vsiz);

			if ((.995 <= iwb.pixAspect) & (iwb.pixAspect <= 1.005))
				iwb.pixAspect = 1;	// let it go if < 0.5%
		}
	}
	sprintf(dmessage_buf, "Writing %s '%s'", itype, imgFile);
	DMESG(DMCtrace, dmessage_buf);
	long	nbwritten = (*imgWriter[i]->WriteImage)(imgFile, &iwb);
	if (nbwritten <= 0)
		DMESGF(DMCresource, "Error writing image file '%s'", imgFile);
	return nbwritten;
}

// Reorient image from "orienta" (or metadata if 0)
bool
PanImage::Reorient(int orienta, bool force)
{
	if (!ims.img)
		return false;
					// figure out what to do
	if (orienta <= 0) {
		if (!GetInfo(IIForientation))
			return true;
		orienta = iinfo->orientation;
	}
	if (orienta == IOtopleft)	// already oriented properly?
		return true;
	const int	psiz = PSize();
	if (ro | (ims.rowsize != ims.xres*psiz) || NLinks() > 1) {
		if (!force && (ro || NLinks() > 1)) {
			DMESGF(DMCparameter, "Reorient called on %s image",
					ro ? "locked" : "linked");
			return false;
		} else if (!Delink())	// force || coalesce private data
			return false;
	}
	bool		hflip=false, vflip=false;
	bool		xyswap=false;
	switch (orienta) {
	case IOtopright:
		hflip = true;
		break;
	case IObotright:
		hflip = vflip = true;
		break;
	case IObotleft:
		if (!force) {		// we're really cheating here...
			ims.img += (ssize_t)(ims.yres-1)*ims.rowsize;
			ims.rowsize = -ims.rowsize;
			if (GetInfo(IIForientation))
				iinfo->orientation = IOtopleft;
			return true;
		}
		vflip = true;
		break;
	case IOlefttop:
		xyswap = true;
		break;
	case IOrighttop:
		xyswap = vflip = true;
		break;
	case IOrightbot:
		xyswap = hflip = vflip = true;
		break;
	case IOleftbot:
		xyswap = hflip = true;
		break;
	default:
		DMESGF(DMCparameter, "Bad orientation value %d", orienta);
		return false;
	}
#define MapToOrig(xo, yo, xc, yc) \
		if (xyswap) { xo = yc; yo = xc; } \
		else { xo = xc; yo = yc; } \
		if (hflip) xo = ims.xres-1 - xo; \
		if (vflip) yo = ims.yres-1 - yo;
	/*
	 * The following code copies pixels in natural loops,
	 * but requires that the source and destination images
	 * occupy exactly the same memory footprint.
	 */
	uby8		ptmp[MaxPixelLen*MaxDataSize];
	int		xs, ys, xd, yd;
	const uby8 *	sptr;
	ImgStruct	idst = ims;
	if (xyswap) { idst.xres = ims.yres; idst.yres = ims.xres; }
	idst.rowsize = idst.xres*psiz;
	ABitMap2	pmark(idst.xres, idst.yres);
	for (xd = yd = 0; pmark.Find(&xd, &yd, false); xd++) {
					// start copy loop
		memcpy(ptmp, PpixP(&idst,xd,yd,psiz), psiz);
		pmark.Set(xd, yd);
		int	xd1 = xd, yd1 = yd;
		int	xd0, yd0;
		for ( ; ; ) {		// run full circle
			xd0 = xd1; yd0 = yd1;
			MapToOrig(xs, ys, xd0, yd0);
			sptr = PpixP(&ims,xs,ys,psiz);
			if (!PpixPos(&xd1, &yd1, &idst, sptr))
				DMESG(DMCassert, "PpixPos failed!");
			if (!pmark.TestAndSet(xd1, yd1))
				break;
			memcpy(PpixP(&idst,xd0,yd0,psiz), sptr, psiz);
		}
					// complete loop
		memcpy(PpixP(&idst,xd0,yd0,psiz), ptmp, psiz);
	}
#undef MapToOrig
	ims = idst;
	if (!GetInfo())
		return true;
					// fix metadata
	iinfo->orientation = IOtopleft;
	iinfo->flags |= IIForientation;
	if (xyswap && GetInfo(IIFhvangle))
		iinfo->hvangle = 360./M_PI*atan( (double)ims.yres/ims.xres *
						tan(M_PI/360.*iinfo->hvangle) );
	return true;
}

// Determine if rectangle covers entire image
static inline bool
rectIsFull(const ImgRect *r, const ImgStruct &im)
{
	if (!im.img)
		return false;
	if (!r)				// convention for entire image
		return true;
	if ((r->xleft > 0) | (r->ytop > 0))
		return false;
	return ((r->xright == im.xres) & (r->ybottom == im.yres));
}

// Link (sub)image from writeable ImgStruct
bool
PanImage::Link(ImgStruct *im, const ImgRect *r)
{
	if (im == &ims) {	// subimage of self
		if (!ims.img)
			return false;
		if (rectIsFull(r, ims))
			return true;
		ClearInfo(IIFcrop|IIFhvangle|IIFview);
		return PlinkSubimage(&ims, &ims, r);
	}
	Init();
	if (!im || !im->img)
		return false;
	return PlinkSubimage(&ims, im, r);
}

// Link (sub)image from PanImage object
bool
PanImage::Link(PanImage *srcp, const ImgRect *r)
{
	if (!srcp) {
		Init();
		return false;
	}
	if (rectIsFull(r, srcp->ims))
		r = NULL;
				// unsets strict read-only if not linked
	srcp->SetReadOnly(srcp->ReadOnly());
	if (!Link(&srcp->ims, r))
		return false;
	if (srcp->GetInfo()) {
		srcp->GetInfo(Info());
		if (r) ClearInfo(IIFcrop|IIFhvangle|IIFview);
	}
	ro = srcp->ro;		// copy read-only status
	mathFlags = srcp->mathFlags;
	return true;
}

// Link (sub)image from read-only ImgStruct
bool
PanImage::LinkRO(const ImgStruct *im, const ImgRect *r)
{
	if (im == &ims) {	// link to self(!)
		SetReadOnly();
		return Link(&ims, r);
	}
	Init();
	if (!im || !im->img)
		return false;
	ro = -1;		// strict read-only
	return PlinkSubimage(&ims, im, r);
}

// Link (sub)image from read-only PanImage object
bool
PanImage::LinkRO(const PanImage &src, const ImgRect *r)
{
	if (rectIsFull(r, src.ims))
		r = NULL;
	if (!LinkRO(&src.ims, r))
		return false;
	if (src.GetInfo()) {
		src.GetInfo(Info());
		if (r) ClearInfo(IIFcrop|IIFhvangle|IIFview);
	}
	mathFlags = src.mathFlags;
	return true;
}

// Crop image using given rectangle corners (more forgiving)
bool
PanImage::Crop(int x1, int y1, int x2, int y2)
{
	ImgRect	r;

	if (!ims.img)
		return false;

	PsetRect(&r, x1, y1, x2, y2);

	if (r.xleft < 0) r.xleft = 0;
	if (r.ytop < 0) r.ytop = 0;
	if (r.xright > ims.xres) r.xright = ims.xres;
	if (r.ybottom > ims.yres) r.ybottom = ims.yres;
	return Crop(&r);
}

// Image copy operator (clears read-only if not self)
PanImage &
PanImage::operator=(const PanImage &src)
{
	if (&src == this)
		return *this;
	if (src.ims.img)
		Init(src.ims.xres, src.ims.yres, *src.ims.csp, PICready);
	else if (src.ims.csp)
		Init(*src.ims.csp);
	else
		Init();
	if (!ims.img)
		return *this;
	PmapImage(&ims, &src.ims, 1.f);
	if (src.GetInfo())
		src.GetInfo(Info());
	mathFlags = src.mathFlags;
	return *this;
}

static inline int
iround(double x)
{
	return int(x + .5) - int(x < -.5);
}

// Draw pixel-wide line between image points (no antialiasing)
int
PanImage::DrawLine(int x0, int y0, int x1, int y1, PixelVal p)
{
	if (Ready() != PICready)
		return 0;

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
	p = PconvertPixel(p, ims.csp);
	while (n--)
		np += SetPixel(iround(x0 + n*xstep), iround(y0 + n*ystep), p);

	return np;
}

// Scalar assignment operator
PanImage &
PanImage::operator=(const PixelVal &pval)
{
	switch (Ready()) {
	case PICvoid:
	case PICspace:
		DMESG(DMCparameter, "Cannot assign value to uninitialized image");
		break;
	case PICread:
		DMESG(DMCparameter, "Cannot assign value to locked image");
		break;
	default:
		PsetImage(&ims, pval);
		break;
	}
	return *this;
}

// Scalar addition operator
PanImage &
PanImage::operator+=(const PixelVal &pval)
{
	if (!pval.csp)
		return *this;
	switch (Ready()) {
	case PICfree:
		return *this = pval;
	case PICvoid:
	case PICspace:
		DMESG(DMCparameter, "Cannot add to uninitialized image");
		return *this;
	case PICread:
		DMESG(DMCparameter, "Cannot operate on locked image");
		return *this;
	default:
		break;
	}
	if (ims.csp->dtype != IDTfloat) {
		DMESG(DMCparameter, "Scalar addition requires real image");
		return *this;
	}
	if (mathFlags & PIMstrictcolor &&
			!PmatchColorSpace(ims.csp, &ICS_Y, PICMgamma)) {
		DMESG(DMCparameter,
			"Unsupported scalar addition on non-linear image");
		return *this;
	}
	const float	eps = (pval.csp->dtype == IDTfloat) ? 1e-6f : 0.003f;
	const PixelVal	padd = PconvertPixel(pval, ims.csp);
	const int	plen = ImgPixelLen[ims.csp->format];
	int		y, n, i;
	for (i = plen; i--; )
		if ((-eps > padd.v.f[i]) | (padd.v.f[i] > eps))
			break;
	if (i < 0)
		return *this;
	for (y = ims.yres; y--; ) {
		float *	fp = (float *)ProwPtr(&ims, y);
		for (n = ims.xres; n--; )
			for (i = 0; i < plen; i++)
				*fp++ += padd.v.f[i];
	}
	return *this;
}

// Scalar subtraction operator
PanImage &
PanImage::operator-=(const PixelVal &pval)
{
	if (!pval.csp)
		return *this;
	if (!ims.csp || ims.csp->dtype != IDTfloat) {
		DMESG(DMCparameter, "Scalar subtraction requires real image");
		return *this;
	}
	PixelVal	psub = PconvertPixel(pval, ims.csp);
	int		i = ImgPixelLen[ims.csp->format];
	while (i--)
		psub.v.f[i] = -psub.v.f[i];
	return *this += psub;
}

// Scalar multiplication operator (color)
PanImage &
PanImage::operator*=(const PixelVal &pval)
{
	switch (Ready()) {
	case PICfree:
		return *this = pval;
	case PICvoid:
	case PICspace:
		DMESG(DMCparameter, "Cannot multiply uninitialized image");
		return *this;
	case PICread:
		DMESG(DMCparameter, "Cannot operate on locked image");
		return *this;
	default:
		break;
	}
	if (!pval.csp) {
		Init(PICzero);
		return *this;
	}
	if (ims.csp->dtype != IDTfloat || (mathFlags & PIMstrictcolor &&
						ims.csp->logorig > 0)) {
		DMESG(DMCparameter, "Scalar multiplication requires real image");
		return *this;
	}
	const float	eps = (pval.csp->dtype == IDTfloat) ? 1e-6f : 0.003f;
	const PixelVal	pmul = PconvertPixel(pval, ims.csp);
	const int	plen = ImgPixelLen[ims.csp->format];
	int		y, n, i;
	for (i = plen; i--; )
		if ((1.f-eps > pmul.v.f[i]) | (pmul.v.f[i] > 1.f+eps))
			break;
	if (i < 0)
		return *this;
	for (y = ims.yres; y--; ) {
		float *	fp = (float *)ProwPtr(&ims, y);
		for (n = ims.xres; n--; )
			for (i = 0; i < plen; i++)
				*fp++ *= pmul.v.f[i];
	}
	if (GetInfo(IIFstonits))		// adjust stonits
		iinfo->stonits /= PgetY(pmul);
	return *this;
}

// Scalar multiplication operator (grayscale)
PanImage &
PanImage::operator*=(const float fval)
{
	switch (Ready()) {
	case PICfree:
		return *this = fval;
	case PICvoid:
	case PICspace:
		DMESG(DMCparameter, "Cannot multiply uninitialized image");
		return *this;
	case PICread:
		DMESG(DMCparameter, "Cannot operate on locked image");
		return *this;
	default:
		break;
	}
	if (fval == 0) {
		Init(PICzero);
		return *this;
	}
	if ((0.999999f <= fval) & (fval <= 1.000001f))
		return *this;
	if (ims.csp->dtype != IDTfloat || (mathFlags & PIMstrictcolor &&
						ims.csp->logorig > 0)) {
		ConvertCS(*ims.csp, fval);	// adjusts stonits, also
		return *this;
	}
	const float	fvadj = (mathFlags & PIMstrictcolor) ?
					(float)pow(fval, 1./ims.csp->gamma) :
					fval ;
	for (int y = ims.yres; y--; ) {
		float *	fp = (float *)ProwPtr(&ims, y);
		for (int n = ims.xres*ImgPixelLen[ims.csp->format]; n--; )
			*fp++ *= fvadj;
	}
	if (GetInfo(IIFstonits))
		iinfo->stonits /= fval;		// adjust stonits
	return *this;
}

// Scalar division operator (color)
PanImage &
PanImage::operator/=(const PixelVal &pval)
{
	ImgColorSpace	destCS;
	PixelVal	pdiv;
	int		i;
	if (!pval.csp)
		goto divByZero;
	PcopyCS(&destCS, pval.csp);
	PrealCS(&destCS);
	pdiv = PconvertPixel(pval, &destCS);
	if (!pval.csp)
		goto divByZero;
	i = ImgPixelLen[destCS.format];
	while (i--) {
		if (pdiv.v.f[i] == 0)
			goto divByZero;
		pdiv.v.f[i] = 1.f/pdiv.v.f[i];
	}
	return *this *= pdiv;
divByZero:
	DMESG(MATHWARN, "Image divide by zero component");
	Init(PICzero);
	return *this;
}

// Scalar division operator (grayscale)
PanImage &
PanImage::operator/=(const float fval)
{
	if (fval == 0) {
		if (Ready() != PICready)
			return *this;
		DMESG(MATHWARN, "Image divide by zero scalar");
		Init(PICzero);
		return *this;
	}
	return *this *= 1.f/fval;
}

// Clamp to some minimum
PanImage &
PanImage::operator|=(const PixelVal &pval)
{
	switch (Ready()) {
	case PICfree:
		return *this = pval;
	case PICvoid:
	case PICspace:
		DMESG(DMCparameter, "Scalar operation on uninitialized image");
		return *this;
	case PICread:
		DMESG(DMCparameter, "Cannot operate on locked image");
		return *this;
	default:
		break;
	}
	const PixelVal	pmin = PconvertPixel(pval, ims.csp);
	const int	plen = ImgPixelLen[ims.csp->format];
	int		y, n, i;
	for (y = ims.yres; y--; )
		switch (ims.csp->dtype) {
		case IDTfloat: {
			float *			fp = (float *)ProwPtr(&ims, y);
			for (n = ims.xres; n--; fp += plen)
				for (i = plen; i--; )
					if (fp[i] < pmin.v.f[i])
						fp[i] = pmin.v.f[i];
			} break;
		case IDTubyte: {
			uby8 *			bp = ProwPtr(&ims, y);
			for (n = ims.xres; n--; bp += plen)
				for (i = plen; i--; )
					if (bp[i] < pmin.v.b[i])
						bp[i] = pmin.v.b[i];
			} break;
		case IDTushort: {
			unsigned short *	sp = (unsigned short *)ProwPtr(&ims, y);
			for (n = ims.xres; n--; sp += plen)
				for (i = plen; i--; )
					if (sp[i] < pmin.v.s[i])
						sp[i] = pmin.v.s[i];
			} break;
		default:
			DMESG(DMCparameter,
				"Unsupported color space for image operation");
			return *this;
		}
	return *this;
}

// Clamp to some maximum
PanImage &
PanImage::operator&=(const PixelVal &pval)
{
	switch (Ready()) {
	case PICfree:
		return *this = pval;
	case PICvoid:
	case PICspace:
		DMESG(DMCparameter, "Scalar operation on uninitialized image");
		return *this;
	case PICread:
		DMESG(DMCparameter, "Cannot operate on locked image");
		return *this;
	default:
		break;
	}
	const PixelVal	pmax = PconvertPixel(pval, ims.csp);
	const int	plen = ImgPixelLen[ims.csp->format];
	int		y, n, i;
	for (y = ims.yres; y--; )
		switch (ims.csp->dtype) {
		case IDTfloat: {
			float *			fp = (float *)ProwPtr(&ims, y);
			for (n = ims.xres; n--; fp += plen)
				for (i = plen; i--; )
					if (fp[i] > pmax.v.f[i])
						fp[i] = pmax.v.f[i];
			} break;
		case IDTubyte: {
			uby8 *			bp = ProwPtr(&ims, y);
			for (n = ims.xres; n--; bp += plen)
				for (i = plen; i--; )
					if (bp[i] > pmax.v.b[i])
						bp[i] = pmax.v.b[i];
			} break;
		case IDTushort: {
			unsigned short *	sp = (unsigned short *)ProwPtr(&ims, y);
			for (n = ims.xres; n--; sp += plen)
				for (i = plen; i--; )
					if (sp[i] > pmax.v.s[i])
						sp[i] = pmax.v.s[i];
			} break;
		default:
			DMESG(DMCparameter,
				"Unsupported color space for image operation");
			return *this;
		}
	return *this;
}

// Dilate image using disk of the given radius
PanImage &
PanImage::operator<<=(const float rad)
{
	switch (Ready()) {
	case PICvoid:
	case PICspace:
	case PICfree:
		DMESGF(DMCparameter, "%s on empty image",
				(rad>0) ? "Dilation" : "Erosion");
		return *this;
	case PICread:
		DMESG(DMCparameter, "Cannot operate on locked image");
		return *this;
	default:
		break;
	}
	PdilateImage(&ims, &ims, rad);
	return *this;
}

// Gaussian blur (or sharpen if rad < 0.4)
PanImage &
PanImage::operator^=(const float rad)
{
	switch (Ready()) {
	case PICvoid:
	case PICspace:
	case PICfree:
		DMESG(DMCparameter, "Blur/sharpen on empty image");
		return *this;
	case PICread:
		DMESG(DMCparameter, "Cannot operate on locked image");
		return *this;
	default:
		break;
	}
	if (rad >= 0.45f) {		// blur
		PblurImage(&ims, &ims, rad);
		return *this;
	}
	if (rad > 0.39f)		// don't bother
		return *this;
					// else sharpen
	float		usmData[3][3] = {
				-.75, -1, -.75,
				  -1,  7,   -1,
				-.75, -1, -.75
			};
	PanImage	usmKern(3, 3, ICS_Y, PICfree);
	usmKern.Img()->img = (uby8 *)usmData;
				// boost effect based on value
	usmKern *= 10.f*(.5f - rad)*(.5f - rad);
	if (rad > 0)		// rad <= 0 gives edge detect
		usmData[1][1] += 1;

	return *this ^= usmKern;
}

// Sum in another image
PanImage &
PanImage::operator+=(const PanImage &src)
{
	if (src.Ready() < PICread) {
		DMESG(DMCparameter, "Addition of uninitialized image");
		return *this;
	}
	switch (Ready()) {
	case PICfree:
		ims.xres = src.ims.xres;
		ims.yres = src.ims.yres;
		// fall through
	case PICvoid:
	case PICspace:
		Load(src);
		return *this;
	case PICread:
		DMESG(DMCparameter, "Cannot operate on locked image");
		return *this;
	default:
		break;
	}
	if ((ims.xres != src.ims.xres) | (ims.yres != src.ims.yres)) {
		DMESG(DMCparameter, "Image addition resolution mismatch");
		return *this;
	}
	if ((ims.csp->dtype != IDTfloat) | (src.ims.csp->dtype != IDTfloat) ||
			((mathFlags & PIMstrictcolor) ?
				!PmatchColorSpace(ims.csp, src.ims.csp, PICMall) :
				ImgPixelLen[ims.csp->format] !=
					ImgPixelLen[src.ims.csp->format])) {
		PsumImage(&ims, &src.ims, 1.f);
		return *this;
	}
	for (int y = ims.yres; y--; ) {
		const float *	sp = (const float *)src[y];
		float *		dp = (float *)ProwPtr(&ims, y);
		for (int n = ims.xres*ImgPixelLen[ims.csp->format]; n--; )
			*dp++ += *sp++;
	}
	return *this;
}

// Subtract another image
PanImage &
PanImage::operator-=(const PanImage &src)
{
	if (src.Ready() < PICread) {
		DMESG(DMCparameter, "Subtraction of uninitialized image");
		return *this;
	}
	switch (Ready()) {
	case PICvoid:
	case PICspace:
	case PICfree:
		DMESG(DMCparameter, "Subtraction from empty image");
		return *this;
	case PICread:
		DMESG(DMCparameter, "Cannot operate on locked image");
		return *this;
	default:
		break;
	}
	if ((ims.xres != src.ims.xres) | (ims.yres != src.ims.yres)) {
		DMESG(DMCparameter, "Image subtraction resolution mismatch");
		return *this;
	}
	if (src.ims.img == ims.img) {
		Init(PICzero);
		return *this;
	}
	if ((ims.csp->dtype != IDTfloat) | (src.ims.csp->dtype != IDTfloat) ||
			((mathFlags & PIMstrictcolor) ?
				!PmatchColorSpace(ims.csp, src.ims.csp, PICMall) :
				ImgPixelLen[ims.csp->format] !=
					ImgPixelLen[src.ims.csp->format])) {
		PsumImage(&ims, &src.ims, -1.f);
		return *this;
	}
	for (int y = ims.yres; y--; ) {
		const float *	sp = (const float *)src[y];
		float *		dp = (float *)ProwPtr(&ims, y);
		for (int n = ims.xres*ImgPixelLen[ims.csp->format]; n--; )
			*dp++ -= *sp++;
	}
	return *this;
}

// Check image source for operator
bool
PanImage::OperatorChecks(const PanImage &src, int lineno) const
{
	if (src.Ready() < PICread) {
		dmessage(DMCparameter, "Uninitialized image given to operator",
				__FILE__, lineno);
		return false;
	}
	if ((ims.xres != src.ims.xres) | (ims.yres != src.ims.yres)) {
		dmessage(DMCparameter, "Image operand resolution mismatch",
				__FILE__, lineno);
		return false;
	}
	if (ims.csp->dtype != IDTfloat) {
		dmessage(DMCparameter, "Operator requires real destination",
				__FILE__, lineno);
		return false;
	}
	return true;
}

// Callback for PanImage::operator*=(const PanImage &)
static int
mulSL(const uby8 *slp, int len, int y, void *udp)
{
	const float *	ps = (const float *)slp;
	ImgStruct *	ib = (ImgStruct *)udp;
	int		n = len * ImgPixelLen[ib->csp->format];
	float *		pd = (float *)ProwPtr(ib, y);

	while (n--)
		*pd++ *= *ps++;
	return 0;
}

// Callback for PanImage::operator*=(const PanImage &) into log image
static int
addSL(const uby8 *slp, int len, int y, void *udp)
{
	const float *	ps = (const float *)slp;
	ImgStruct *	ib = (ImgStruct *)udp;
	int		n = len * ImgPixelLen[ib->csp->format];
	float *		pd = (float *)ProwPtr(ib, y);

	while (n--)
		*pd++ += *ps++;
	return 0;
}

// Multiply in another image
PanImage &
PanImage::operator*=(const PanImage &src)
{
	PscanlineMethod *	cvtSL = mulSL;
	ImgColorSpace		myCS;

	myCS.cstatic = -1;

	switch (Ready()) {
	case PICvoid:			// need to initialize
		PcopyCS(&myCS, src.ims.csp);
	    // fall through
	case PICspace:
	case PICfree:
		if (src.Ready() < PICread)
			return *this;
		if (myCS.cstatic < 0)
			PcopyCS(&myCS, ims.csp);
		PrealCS(&myCS);
		Init(myCS);
		Load(src);
		return *this;
	case PICread:
		DMESG(DMCparameter, "Cannot operate on locked image");
		return *this;
	default:
		break;
	}
	if (!OperatorChecks(src, __LINE__))
		return *this;
	if (mathFlags & PIMstrictcolor) {
		PcopyCS(&myCS, ims.csp); // maintain strict color
		if (myCS.logorig > 0) {
			myCS.logorig = 1;
			cvtSL = addSL;
		}
	} else {			// or just convert pixel format
		PcopyCS(&myCS, src.ims.csp);
		myCS.dtype = ims.csp->dtype;
		myCS.format = ims.csp->format;
	}
	PconvertImageCB(&src.ims, &myCS, 1.f, cvtSL, &ims);
	ClearInfo(IIFstonits);		// invalidates calibration
	return *this;
}

// Callback for PanImage::operator/=(const PanImage &)
static int
divSL(const uby8 *slp, int len, int y, void *udp)
{
	const float *	ps = (const float *)slp;
	ImgStruct *	ib = (ImgStruct *)udp;
	int		n = len * ImgPixelLen[ib->csp->format];
	float *		pd = (float *)ProwPtr(ib, y);
	int		nz = 0;

	while (n--)
		if (*ps == 0) {
			*pd++ = 0; ps++;
			++nz;
		} else
			*pd++ /= *ps++;
	return nz;
}

// Callback for PanImage::operator/=(const PanImage &) into log image
static int
subSL(const uby8 *slp, int len, int y, void *udp)
{
	const float *	ps = (const float *)slp;
	ImgStruct *	ib = (ImgStruct *)udp;
	int		n = len * ImgPixelLen[ib->csp->format];
	float *		pd = (float *)ProwPtr(ib, y);

	while (n--)
		*pd++ -= *ps++;
	return 0;
}

// Divide in another image
PanImage &
PanImage::operator/=(const PanImage &src)
{
	PscanlineMethod *	cvtSL = divSL;
	ImgColorSpace		myCS;

	myCS.cstatic = -1;

	switch (Ready()) {
	case PICvoid:			// need to initialize
		PcopyCS(&myCS, src.ims.csp);
	    // fall through
	case PICspace:
	case PICfree:			// compute reciprocal
		if (src.Ready() < PICread)
			return *this;
		if (myCS.cstatic < 0)
			PcopyCS(&myCS, ims.csp);
		PrealCS(&myCS);
		Init(src.ims.xres, src.ims.yres, myCS, PICfree);
		*this = 1.f;
		if (src.GetInfo())
			src.GetInfo(Info());
		break;
	case PICread:
		DMESG(DMCparameter, "Cannot operate on locked image");
		return *this;
	default:			// check for divide-by-self
		if (src.ims.img == ims.img && (src.ims.xres == ims.xres) &
						(src.ims.yres == ims.yres))
			return *this = 1.f;
		break;
	}
	if (!OperatorChecks(src, __LINE__))
		return *this;
	if (mathFlags & PIMstrictcolor) {
		PcopyCS(&myCS, ims.csp); // maintain strict color
		if (myCS.logorig > 0) {
			myCS.logorig = 1;
			cvtSL = subSL;
		}			
	} else {			// or just convert pixel format
		PcopyCS(&myCS, src.ims.csp);
		myCS.dtype = ims.csp->dtype;
		myCS.format = ims.csp->format;
	}
	int	nzeroes = PconvertImageCB(&src.ims, &myCS, 1.f, cvtSL, &ims);
	if (nzeroes > 0)
		DMESGF(MATHWARN, "Divide-by-zero in %.3f%% of image",
				100.f*nzeroes /
			(ImgPixelLen[ims.csp->format]*ims.xres*ims.yres));
	ClearInfo(IIFstonits);		// invalidates calibration
	return *this;
}

// Callback for float destination of operator|=(const PanImage &)
static int
unionSLf(const uby8 *slp, int len, int y, void *udp)
{
	const float *	ps = (const float *)slp;
	ImgStruct *	ib = (ImgStruct *)udp;
	int		n = len * ImgPixelLen[ib->csp->format];
	float *		pd = (float *)ProwPtr(ib, y);

	while (n--) {
		if (*pd < *ps)
			*pd = *ps;
		++pd; ++ps;
	}
	return 0;
}

// Callback for unsigned byte destination of operator|=(const PanImage &)
static int
unionSLb(const uby8 *ps, int len, int y, void *udp)
{
	ImgStruct *	ib = (ImgStruct *)udp;
	int		n = len * ImgPixelLen[ib->csp->format];
	uby8 *		pd = ProwPtr(ib, y);

	while (n--) {
		if (*pd < *ps)
			*pd = *ps;
		++pd; ++ps;
	}
	return 0;
}

// Callback for unsigned short destination of operator|=(const PanImage &)
static int
unionSLs(const uby8 *slp, int len, int y, void *udp)
{
	const unsigned short *	ps = (const unsigned short *)slp;
	ImgStruct *		ib = (ImgStruct *)udp;
	int			n = len * ImgPixelLen[ib->csp->format];
	unsigned short *	pd = (unsigned short *)ProwPtr(ib, y);

	while (n--) {
		if (*pd < *ps)
			*pd = *ps;
		++pd; ++ps;
	}
	return 0;
}

// Callback for unsupported data type
static int
badType(const uby8 *slp, int len, int y, void *udp)
{
	DMESG(DMCparameter, "Unsupported data type");
	return -1;
}

// Maximum union operator
PanImage &
PanImage::operator|=(const PanImage &src)
{
	static PscanlineMethod	*spa[] = {
		unionSLb, unionSLs, badType, unionSLf
	};

	if (src.Ready() < PICread) {
		DMESG(DMCparameter, "Union of uninitialized image");
		return *this;
	}
	switch (Ready()) {
	case PICfree:
		ims.xres = src.ims.xres;
		ims.yres = src.ims.yres;
		// fall through
	case PICvoid:
	case PICspace:
		Load(src);
		return *this;
	case PICread:
		DMESG(DMCparameter, "Cannot operate on locked image");
		return *this;
	default:
		break;
	}
	if ((ims.xres != src.ims.xres) | (ims.yres != src.ims.yres)) {
		DMESG(DMCparameter, "Image union resolution mismatch");
		return *this;
	}
	if (src.ims.img == ims.img)
		return *this;

	PconvertImageCB(&src.ims, ims.csp, 1.f, spa[ims.csp->dtype], &ims);
	return *this;
}

// Callback for float destination of operator&=(const PanImage &)
static int
interSLf(const uby8 *slp, int len, int y, void *udp)
{
	const float *	ps = (const float *)slp;
	ImgStruct *	ib = (ImgStruct *)udp;
	int		n = len * ImgPixelLen[ib->csp->format];
	float *		pd = (float *)ProwPtr(ib, y);

	while (n--) {
		if (*pd > *ps)
			*pd = *ps;
		++pd; ++ps;
	}
	return 0;
}

// Callback for unsigned byte destination of operator&=(const PanImage &)
static int
interSLb(const uby8 *ps, int len, int y, void *udp)
{
	ImgStruct *	ib = (ImgStruct *)udp;
	int		n = len * ImgPixelLen[ib->csp->format];
	uby8 *		pd = ProwPtr(ib, y);

	while (n--) {
		if (*pd > *ps)
			*pd = *ps;
		++pd; ++ps;
	}
	return 0;
}

// Callback for unsigned short destination of operator&=(const PanImage &)
static int
interSLs(const uby8 *slp, int len, int y, void *udp)
{
	const unsigned short *	ps = (const unsigned short *)slp;
	ImgStruct *		ib = (ImgStruct *)udp;
	int			n = len * ImgPixelLen[ib->csp->format];
	unsigned short *	pd = (unsigned short *)ProwPtr(ib, y);

	while (n--) {
		if (*pd > *ps)
			*pd = *ps;
		++pd; ++ps;
	}
	return 0;
}

// Minimum intersection operator
PanImage &
PanImage::operator&=(const PanImage &src)
{
	static PscanlineMethod	*spa[] = {
		interSLb, interSLs, badType, interSLf
	};

	if (src.Ready() < PICread) {
		DMESG(DMCparameter, "Intersection of uninitialized image");
		return *this;
	}
	switch (Ready()) {
	case PICfree:
		ims.xres = src.ims.xres;
		ims.yres = src.ims.yres;
		// fall through
	case PICvoid:
	case PICspace:
		Load(src);
		return *this;
	case PICread:
		DMESG(DMCparameter, "Cannot operate on locked image");
		return *this;
	default:
		break;
	}
	if ((ims.xres != src.ims.xres) | (ims.yres != src.ims.yres)) {
		DMESG(DMCparameter, "Image intersection resolution mismatch");
		return *this;
	}
	if (src.ims.img == ims.img)
		return *this;

	PconvertImageCB(&src.ims, ims.csp, 1.f, spa[ims.csp->dtype], &ims);
	return *this;
}

// Apply the given convolution kernel to our image
PanImage &
PanImage::operator^=(const PanImage &kern)
{
	switch (Ready()) {
	case PICvoid:
	case PICspace:
	case PICfree:
		DMESG(DMCparameter, "Convolution on empty image");
		return *this;
	case PICread:
		DMESG(DMCparameter, "Cannot operate on locked image");
		return *this;
	default:
		break;
	}
	if (kern.Ready() < PICread) {
		DMESG(DMCparameter, "Empty convolution kernel");
		return *this;
	}
	PanImage	iTmp;		// may need temporary float image
	ImgColorSpace	fltCS;
	PcopyCS(&fltCS, ims.csp);
	PrealCS(&fltCS);
	iTmp.ims.csp = &fltCS;
	if (PmatchColorSpace(ims.csp, &fltCS, PICMdtype|PICMgamma)) {
		iTmp.ims = ims;
		ims.img = NULL;
	} else if (!PmapImage(&iTmp.ims, &ims, 1.f))
		return *this;
					// PconvolveImage() does the rest
	if (!PconvolveImage(&ims, &iTmp.ims, &kern.ims))
		DMESG(DMCdata, "PconvolveImage() failed!");
					// fix calibration if any
	if (GetInfo(IIFstonits)) {
		float	meanY = kern.GetAverage(&ICS_Y).v.f[0];
		if (meanY <= 1e-6f)
			ClearInfo(IIFstonits);
		else
			iinfo->stonits /= meanY *
				float(kern.Width()*kern.Height());
	}
	return *this;
}

// Choose destination color space based on two operands
const ImgColorSpace *
PickCS(const PanImage &im1, const PanImage &im2, bool makeReal)
{
	static ImgColorSpace	myCS;
	const int		nc1 = im1.NComp();
	const int		nc2 = im2.NComp();
	int			ds1, ds2;

	if (!nc1 & !nc2)
		return &ICS_Y;		// doesn't matter
	if (makeReal) {
		if (nc1 == nc2) {
			if ((im1.GetCS()->dtype == IDTfloat) |
					(im2.GetCS()->dtype != IDTfloat))
				PcopyCS(&myCS, im1.GetCS());
			else
				PcopyCS(&myCS, im2.GetCS());
		} else
			PcopyCS(&myCS, (nc1>nc2) ? im1.GetCS() : im2.GetCS());
		PrealCS(&myCS);
		return &myCS;
	}
	if (!nc2) return im1.GetCS();
	if (!nc1) return im2.GetCS();
	ds1 = ImgDataSize[im1.GetCS()->dtype];
	ds2 = ImgDataSize[im2.GetCS()->dtype];
	if (ds1 >= ds2) {
		if (nc1 >= nc2) return im1.GetCS();
		if (ds1 == ds2) return im2.GetCS();
		PcopyCS(&myCS, im2.GetCS());
		myCS.dtype = im1.GetCS()->dtype;
	} else {
		if (nc2 >= nc1) return im2.GetCS();
		PcopyCS(&myCS, im1.GetCS());
		myCS.dtype = im2.GetCS()->dtype;
	}
	return &myCS;
}

// Prepare two pixels by choosing and converting to wider color space
static inline bool
PromoteOperands(PixelVal &p1, PixelVal &p2)
{
	if (!p1.csp | !p2.csp)
		return false;
	if (p1.csp == p2.csp)
		return true;
	int	ncdiff = ImgPixelLen[p2.csp->format] - ImgPixelLen[p1.csp->format];
	if (ncdiff > 0 || !ncdiff &
			(ImgDataSize[p2.csp->dtype] > ImgDataSize[p1.csp->dtype]))
		p1 = PconvertPixel(p1, p2.csp);
	else
		p2 = PconvertPixel(p2, p1.csp);
	return true;
}

// Add two pixels
PixelVal
operator+(PixelVal p1, PixelVal p2)
{
	if (!p2.csp)
		return p1;
	if (!p1.csp)
		return p2;
	if (!PromoteOperands(p1, p2))
		return Pblack;

	int	i = ImgPixelLen[p1.csp->format];

	switch (p1.csp->dtype) {
	case IDTfloat:
		while (i--) p1.v.f[i] += p2.v.f[i];
		break;
	case IDTubyte:
		while (i--)
			if (p1.v.b[i] + p2.v.b[i] >= 0xff) p1.v.b[i] = 0xff;
			else p1.v.b[i] += p2.v.b[i];
		break;
	case IDTushort:
		while (i--)
			if (p1.v.s[i] + p2.v.s[i] >= 0xffff) p1.v.s[i] = 0xffff;
			else p1.v.s[i] += p2.v.s[i];
		break;
	default:
		return Pblack;
	}
	return p1;
}

// Subtract two pixels
PixelVal
operator-(PixelVal p1, PixelVal p2)
{
	if (!p2.csp)
		return p1;
	if (!p1.csp)
		return -p2;
	if (!PromoteOperands(p1, p2))
		return Pblack;

	int	i = ImgPixelLen[p1.csp->format];

	switch (p1.csp->dtype) {
	case IDTfloat:
		while (i--) p1.v.f[i] -= p2.v.f[i];
		break;
	case IDTubyte:
		while (i--)
			if (p1.v.b[i] <= p2.v.b[i]) p1.v.b[i] = 0;
			else p1.v.b[i] -= p2.v.b[i];
		break;
	case IDTushort:
		while (i--)
			if (p1.v.s[i] <= p2.v.s[i]) p1.v.s[i] = 0;
			else p1.v.s[i] -= p2.v.s[i];
		break;
	default:
		return Pblack;
	}
	return p1;
}

// Multiply two pixels
PixelVal
operator*(PixelVal p1, PixelVal p2)
{
	if (!PromoteOperands(p1, p2))
		return Pblack;

	int	i = ImgPixelLen[p1.csp->format];

	switch (p1.csp->dtype) {
	case IDTfloat:
		while (i--) p1.v.f[i] *= p2.v.f[i];
		break;
	case IDTubyte:
		while (i--)
			p1.v.b[i] = ((p1.v.b[i] + 1) * p2.v.b[i]) >> 8;
		break;
	case IDTushort:
		while (i--)
			p1.v.s[i] = (unsigned(p1.v.s[i] + 1) * p2.v.s[i]) >> 16;
		break;
	default:
		return Pblack;
	}
	return p1;
}

// Multiply pixel by scalar
PixelVal
operator*(PixelVal p, float f)
{
	if (f == 1)
		return p;
	if (!p.csp | (f == 0))
		return Pblack;
	if ((p.csp->dtype != IDTfloat) & (f < 0))
		return Pblack;

	int	i = ImgPixelLen[p.csp->format];

	switch (p.csp->dtype) {
	case IDTfloat:
		while (i--) p.v.f[i] *= f;
		break;
	case IDTubyte:
		while (i--)
			if (p.v.b[i]*f >= 0xff) p.v.b[i] = 0xff;
			else p.v.b[i] = int(p.v.b[i]*f + .5f);
		break;
	case IDTushort:
		while (i--)
			if (p.v.s[i]*f >= 0xffff) p.v.s[i] = 0xffff;
			else p.v.s[i] = int(p.v.s[i]*f + .5f);
		break;
	default:
		return Pblack;
	}
	return p;
}
