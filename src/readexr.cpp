/*
 *  readexr.cpp
 *  
 *  Open and read ILM EXR high dynamic-range RGBA format (16-bit floats).
 *
 *  Created by gward on Fri Aug 23 2002.
 *  Copyright (c) 2002 Anyhere Software. All rights reserved.
 *
 */

#ifdef freebsd
#define PLATFORM_DARWIN_PPC
#endif

#include "rtio.h"
#include <stdlib.h>
#include <math.h>
#include "imgreader.h"
#include <ImfRgbaFile.h>
#include <ImfStandardAttributes.h>

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

// Tone-mapping data conversion routines
#include "color.h"
#include "tonemap.h"

#ifndef EXR_MAXHISTO
#define EXR_MAXHISTO	(512L*512L)	/* maximum pixels to histogram */
#endif

using namespace Imf;
using namespace std;

// Floating point comparison of primary chromaticities
#define	FEQ(a,b)	((a) < (b)+1e-5 && (b) < (a)+1e-5)
#define	PRIMEQ(p1,p2)	(FEQ((p1)[0][0],(p2)[0][0])&&FEQ((p1)[0][1],(p2)[0][1])\
			&&FEQ((p1)[1][0],(p2)[1][0])&&FEQ((p1)[1][1],(p2)[1][1])\
			&&FEQ((p1)[2][0],(p2)[2][0])&&FEQ((p1)[2][1],(p2)[2][1])\
			&&FEQ((p1)[3][0],(p2)[3][0])&&FEQ((p1)[3][1],(p2)[3][1]))

// Conversion buffer for tone-mapping operations
struct MyTMbuf {
	int			bsiz;		// number of pixels allocated
	bool			grayscale;	// allocated for grayscale
	TMbright *		lbuf;		// luminance buffer for TM
	float *			fbuf;		// float pixel data buffer
				MyTMbuf() {
					lbuf=NULL; fbuf=NULL;
					grayscale = false; bsiz = 0;
				}
				~MyTMbuf() {
					Free();
				}
	bool			IncreaseSize(int npix, bool gs = false) {
					if ((npix <= bsiz) & (gs >= grayscale))
						return true;
					Free();
					lbuf = new TMbright [npix];
					if ((grayscale = gs))
						fbuf = new float [npix];
					else
						fbuf = new float [3*npix];
					bsiz = npix;
					return true;
				}
	void			Free() {
					delete [] lbuf; lbuf = NULL;
					delete [] fbuf; fbuf = NULL;
					bsiz = 0;
				}
};

// Derive our input class from ImgReader struct
struct ExrReader : ImgReader {
private:
	char			encdesc[64];	// encoding description
	RgbaInputFile		exr;		// OpenEXR library object
	Rgba *			ibuf;		// EXR input buffer
	int			stripy;		// top y of loaded strip
	int			striplen;	// scanlines per compressed strip
	RGBPRIMS		monprims;	// target display primaries
	double			LdDyn;		// dynamic range of display
	double			LdMax;		// maximum display luminance
	bool			humanVis;	// match human visibility?
	TMstruct *		ts;		// tone-mapping structure
	MyTMbuf			tb;		// tone-mapping buffers
	bool			FillRGB96(ImgReadBuf *rb);
	bool			FillY32(ImgReadBuf *rb);
	bool			FillA32(ImgReadBuf *rb);
	bool			FillRGB24(ImgReadBuf *rb);
	bool			FillY8(ImgReadBuf *rb);
	bool			(ExrReader::*fillF)(ImgReadBuf *rb);
	bool			ReadStrip();
	bool			FloatSpan(float *fp, int x, int y, int len,
						int rate = 1);
	bool			FillFloat(float *fbuf, ImgReadBuf *rb);
	bool			SetConversion(const ImgColorSpace *ocs);
public:
				ExrReader(const char fname[]);
				~ExrReader() {
					if (ts != NULL) tmDone(ts);
					delete [] ibuf;
				}
	bool			GetInfo(ImgInfo *info);
	void			ToneMapping(double dyn, double max, int hv);
	bool			ReadRec(ImgReadBuf *rb);
};

// Open EXR image and return reader object
static ImgReader *
ERopen(const char *fname)
{
	if (fname == NULL)
		return NULL;
	ExrReader *	exrr;
	try {
		exrr = new ExrReader(fname);
	}
	catch (const exception &e) {
#if 1
		cerr << e.what() << '\n';
#endif
		return NULL;
	}
	return exrr;
}

// Close EXR reader and free object
static void
ERclose(ImgReader *ir)
{
	delete (ExrReader *)ir;
}

// Get EXR information
static ImgReadErr
ERgetInfo(ImgReader *ir, ImgInfo *info)
{
	ir->errCode = IREnone;
	if (!(*(ExrReader *)ir).GetInfo(info) && !ir->errCode)
		ir->errCode = IREunknown;
	return ir->errCode;
}

// Set EXR tone-mapping
static void
ERtoneMapping(ImgReader *ir, double dyn, double max, int hv)
{
	(*(ExrReader *)ir).ToneMapping(dyn, max, hv);
}

// Read EXR rectangle
static ImgReadErr
ERreadRec(ImgReader *ir, ImgReadBuf *rb)
{
	ir->errCode = IREnone;
	if (!(*(ExrReader *)ir).ReadRec(rb) && ir->errCode == IREnone) {
		strcpy(ir->errMsg, "EXR reader routine failure");
		ir->errCode = IREunknown;
	}
	return ir->errCode;
}

// EXR reader interface
extern const ImgReaderInterface IRInterfaceEXR;
const ImgReaderInterface IRInterfaceEXR = {
	"EXR.ilm",
	&ERopen, NULL, &ERgetInfo, NULL,
	&ERtoneMapping, &ERreadRec, &ERclose
};

// Open a new EXR file
ExrReader::ExrReader(const char fname[]) : exr(fname)
{
	ri = &IRInterfaceEXR;
	ibuf = NULL; ts = NULL;
	LdDyn = 100.; LdMax = 100.;
	humanVis = false;
	errCode = IREnone;
	strcpy(file, fname);
	xres = exr.dataWindow().max.x - exr.dataWindow().min.x + 1;
	yres = exr.dataWindow().max.y - exr.dataWindow().min.y + 1;
	pixAspect = exr.pixelAspectRatio();
	encoding = encdesc;			// get compression
	striplen = 1;
	switch (exr.compression()) {
	case NO_COMPRESSION:
		strcpy(encdesc, "uncompressed");
		break;
	case RLE_COMPRESSION:
		strcpy(encdesc, "RLE");
		break;
	case ZIP_COMPRESSION:
		striplen = 16;
	    // fall through
	case ZIPS_COMPRESSION:
		strcpy(encdesc, "ZIP");
		break;
	case PIZ_COMPRESSION:
		strcpy(encdesc, "PIZ");
		break;
	default:
		strcpy(encdesc, "compressed");
		break;
	}
	PcopyCS(&cs, &ICS_RGB709);		// get color encoding
	if (hasChromaticities(exr.header())) {
		const Chromaticities &	chrom = chromaticities(exr.header());
		cs.chroma[0][0] = chrom.red[0];
		cs.chroma[0][1] = chrom.red[1];
		cs.chroma[1][0] = chrom.green[0];
		cs.chroma[1][1] = chrom.green[1];
		cs.chroma[2][0] = chrom.blue[0];
		cs.chroma[2][1] = chrom.blue[1];
		cs.chroma[3][0] = chrom.white[0];
		cs.chroma[3][1] = chrom.white[1];
	}
	fillF = &ExrReader::FillRGB96;
	switch (exr.channels()) {
	case WRITE_RGB:
		strcat(encdesc, " RGB");
		break;
	case WRITE_RGBA:
		strcat(encdesc, " RGBA");
		break;
	case WRITE_A:
		strcat(encdesc, " Alpha");
		cs.format = IPFa;
		fillF = &ExrReader::FillA32;
		break;
	case WRITE_R:
		strcat(encdesc, " Red");
		cs.format = IPFy;
		fillF = &ExrReader::FillY32;
		break;
	case WRITE_G:
		strcat(encdesc, " Green");
		cs.chroma[0][0] = cs.chroma[1][0];
		cs.chroma[0][1] = cs.chroma[1][1];
		cs.format = IPFy;
		fillF = &ExrReader::FillY32;
		break;
	case WRITE_Y:
		strcat(encdesc, " Luminance");
		cs.format = IPFy;
		fillF = &ExrReader::FillY32;
		break;
	case WRITE_B:
		strcat(encdesc, " Blue");
		cs.chroma[0][0] = cs.chroma[2][0];
		cs.chroma[0][1] = cs.chroma[2][1];
		cs.format = IPFy;
		fillF = &ExrReader::FillY32;
		break;
	default:
		strcat(encdesc, " ???");
		sprintf(errMsg, "Unsupported input color space (%d)", exr.channels());
		errCode = IREunsupported;
		return;
	}
	fr.xleft = 0; fr.xright = xres;		// initialize scan position
	switch (exr.lineOrder()) {
	case INCREASING_Y:			// scans proceed downwards
		fr.ytop = 0;
		break;
	case DECREASING_Y:			// scans proceed upwards
		fr.ytop = yres - striplen;
		break;
	default:
		sprintf(errMsg, "Unsupported scanline order (%d)", exr.lineOrder());
		errCode = IREunsupported;
		return;
	}
	fr.ybottom = fr.ytop + striplen;
	nr = fr;
	ibuf = new Rgba [xres*striplen];	// allocate our scanline buffer
	if (striplen == 1)			// assign once if we can
		exr.setFrameBuffer(ibuf, 1, 0);
	stripy = -1000;
	frame = 0; nframes = 1;
	frameType = IRFnone; frameRate = 0;
}

// Set new tone-mapping parameters
void
ExrReader::ToneMapping(double dyn, double max, int hv)
{
	LdDyn = dyn;			// store new parameters
	LdMax = max;
	humanVis = (bool)hv;
	if (ts != NULL) {		// reset tone-mapping
		tmDone(ts); ts = NULL;
	}
}

// Get header information
bool
ExrReader::GetInfo(ImgInfo *info)
{
	if (info == NULL)
		return false;
	*info = defImgInfo;
	if (exr.displayWindow() != exr.dataWindow()) {
		info->crop.xleft = exr.displayWindow().min.x -
					exr.dataWindow().min.x;
		info->crop.xright = exr.displayWindow().max.x -
					exr.dataWindow().min.x + 1;
		info->crop.ytop = exr.displayWindow().min.y -
					exr.dataWindow().min.y;
		info->crop.ybottom = exr.displayWindow().max.y -
					exr.dataWindow().min.y + 1;
		info->flags |= IIFcrop;
	}
	const Header::ConstIterator &	hend = exr.header().end();
	Header::ConstIterator		ha;
	if (hasCapDate(exr.header())) {
		strncpy(info->capdate, capDate(exr.header()).c_str(), 20);
		if (info->capdate[18] && !info->capdate[19])
			info->flags |= IIFcapdate;
	}
	if (hasLongitude(exr.header()) && hasLatitude(exr.header())) {
		info->latlong[0] = latitude(exr.header());
		info->latlong[1] = longitude(exr.header());
		info->flags |= IIFlatlong;
	}
	if (hasWhiteLuminance(exr.header())) {
		info->stonits = whiteLuminance(exr.header());
		if (info->stonits > 0)
			info->flags |= IIFstonits;
	}
	if (hasXDensity(exr.header())) {
		info->hdensity = xDensity(exr.header());
		if (info->hdensity > 0)
			info->flags |= IIFhdensity;
	}
	if ((ha = exr.header().find("screenWindowWidth")) != hend &&
			ha.attribute().typeName() == FloatAttribute::staticTypeName()) {
		info->hvangle = static_cast<const FloatAttribute&>
					(ha.attribute()).value();
		if (info->hvangle > 0) {
			info->hvangle = 2.*180./M_PI * atan(0.5*info->hvangle);
			info->flags |= IIFhvangle;
		}
	}
	if (hasExpTime(exr.header())) {
		info->exptime = expTime(exr.header());
		if (info->exptime > 0)
			info->flags |= IIFexptime;
	}
	if (hasAperture(exr.header())) {
		info->aperture = aperture(exr.header());
		if (info->aperture > 0)
			info->flags |= IIFaperture;
	}
	if (hasIsoSpeed(exr.header())) {
		info->asa = isoSpeed(exr.header());
		if (info->asa > 0)
			info->flags |= IIFasa;
	}
	if (hasFocus(exr.header())) {
		info->focus = focus(exr.header());
		if (info->focus > 0)
			info->flags |= IIFfocus;
	}
	if ((ha = exr.header().find("view")) != hend &&
			ha.attribute().typeName() == StringAttribute::staticTypeName()) {
		strlcpy(info->view,
				static_cast<const StringAttribute&>
					(ha.attribute()).value().c_str(),
					sizeof(info->view));
		info->flags |= IIFview;
	}
	if (hasOwner(exr.header())) {
		strlcpy(info->owner, owner(exr.header()).c_str(),
				sizeof(info->owner));
		info->flags |= IIFowner;
	}
	if (hasComments(exr.header())) {
		strlcpy(info->comments, comments(exr.header()).c_str(),
				sizeof(info->comments));
		info->flags |= IIFcomments;
	}
	return true;
}

// Read the next strip of scanlines
bool
ExrReader::ReadStrip()
{
	if (nr.ytop >= nr.ybottom)
		return false;
	if (striplen > 1)			// really ugly...
		exr.setFrameBuffer(ibuf-(exr.dataWindow().min.y+nr.ytop)*xres, 1, xres);
	try {
		exr.readPixels(exr.dataWindow().min.y+nr.ytop,
				exr.dataWindow().min.y+nr.ybottom-1);
	}
	catch (const exception &e) {
		strcpy(errMsg, e.what());
		errCode = IREread;
		stripy = -1;
		return false;
	}
	stripy = nr.ytop;			// advance to next strip
	switch (exr.lineOrder()) {
	case INCREASING_Y:
		nr.ytop += striplen;
		nr.ybottom += striplen;
		if (nr.ybottom > yres)
			nr.ybottom = yres;
		break;
	case DECREASING_Y:
		nr.ytop -= striplen;
		nr.ybottom -= striplen;
		if (nr.ytop < 0)
			nr.ytop = 0;
		break;
	default:
		return false;
	}
	return true;
}

// Load the specified pixels into a floating-point buffer
bool
ExrReader::FloatSpan(float *fp, int x, int y, int len, int rate)
{
	if ((x < 0) | (x >= xres) | (len <= 0))
		return false;
						// load relevant strip
	if ((y < stripy) | (y >= stripy+striplen)) {
		if ((y < nr.ytop) | (y >= nr.ybottom))
			switch (exr.lineOrder()) {
			case INCREASING_Y:
				nr.ytop = (y / striplen) * striplen;
				nr.ybottom = nr.ytop + striplen;
				if (nr.ybottom > yres)
					nr.ybottom = yres;
				break;
			case DECREASING_Y:
				nr.ybottom = yres - ((yres-1-y) / striplen) * striplen;
				nr.ytop = nr.ybottom - striplen;
				break;
			default:
				return false;
			}
		if (!ReadStrip())
			return false;
	}
	x += (y - stripy)*xres;				// adj. buffer position
	switch (exr.channels()) {
	case WRITE_RGBA:
		if (cs.format == IPFa) {		// getting alpha layer
			for ( ; len--; x += rate)
				*fp++ = ibuf[x].a;
			break;
		}
		for ( ; len--; x += rate) {
			float	anorm = ibuf[x].a;
			if (anorm == 0)
				anorm = 1.f;
			else				// undo premult
				anorm = 1.f / anorm;
			*fp++ = anorm * ibuf[x].r;
			*fp++ = anorm * ibuf[x].g;
			*fp++ = anorm * ibuf[x].b;
		}
		break;
	case WRITE_RGB:
		for ( ; len--; x += rate) {
			*fp++ = ibuf[x].r;
			*fp++ = ibuf[x].g;
			*fp++ = ibuf[x].b;
		}
		break;
	case WRITE_R:
		for ( ; len--; x += rate)
			*fp++ = ibuf[x].r;
		break;
	case WRITE_G:
	case WRITE_Y:
		for ( ; len--; x += rate)
			*fp++ = ibuf[x].g;
		break;
	case WRITE_B:
		for ( ; len--; x += rate)
			*fp++ = ibuf[x].b;
		break;
	case WRITE_A:
		for ( ; len--; x += rate)
			*fp++ = ibuf[x].a;
	default:
		return false;
	}
	return true;
}

// Read strip(s) to fill the given float buffer
bool
ExrReader::FillFloat(float *fbuf, ImgReadBuf *rb)
{
	const int	slen = (rb->r.xright - rb->r.xleft)/rb->subsample;
	const int	plen = ImgPixelLen[cs.format];
	int		ypos, ystep, yend;
	if (fr.ytop) {			// DECREASING_Y
		ypos = rb->r.ybottom-1;
		ystep = -rb->subsample;
		yend = rb->r.ytop-1;
	} else {			// INCREASING_Y
		ypos = rb->r.ytop;
		ystep = rb->subsample;
		yend = rb->r.ybottom;
	}
	for ( ; ypos != yend; ypos += ystep)
		if (!FloatSpan(fbuf + (ypos-rb->r.ytop)/rb->subsample*slen*plen,
				rb->r.xleft, ypos, slen, rb->subsample))
			return false;
	return true;
}

// Read strip(s) to fill the given RGB float buffer
bool
ExrReader::FillRGB96(ImgReadBuf *rb)
{
	if ((rb->cs.dtype != IDTfloat) | (rb->cs.format != IPFrgb))
		return false;
	return FillFloat((float *)rb->buf, rb);
}

// Read strip(s) to fill the given Gray float buffer
bool
ExrReader::FillY32(ImgReadBuf *rb)
{
	if ((rb->cs.dtype != IDTfloat) | (rb->cs.format != IPFy))
		return false;
	return FillFloat((float *)rb->buf, rb);
}

// Read strip(s) to fill the given Alpha float buffer
bool
ExrReader::FillA32(ImgReadBuf *rb)
{
	if ((rb->cs.dtype != IDTfloat) | (rb->cs.format != IPFa))
		return false;
	return FillFloat((float *)rb->buf, rb);
}

// Read strip(s) to fill the given RGB byte buffer
bool
ExrReader::FillRGB24(ImgReadBuf *rb)
{
	if ((rb->cs.dtype != IDTubyte) | (rb->cs.format != IPFrgb))
		return false;
	const int	npix = ((rb->r.xright - rb->r.xleft)/rb->subsample) *
				((rb->r.ybottom - rb->r.ytop)/rb->subsample);
	if (!tb.IncreaseSize(npix, false))
		return false;
	if (!FillFloat(tb.fbuf, rb))
		return false;
	if (tmCvColors(ts, tb.lbuf, rb->buf, (COLOR *)tb.fbuf, npix) != TM_E_OK)
		goto tmerr;
	if (tmMapPixels(ts, rb->buf, tb.lbuf, rb->buf, npix) != TM_E_OK)
		goto tmerr;
	return true;
tmerr:
	strcpy(errMsg, tmErrorMessage[ts->lastError]);
	errCode = IREunknown;
	return false;
}

// Read strip(s) to fill the given Gray byte buffer
bool
ExrReader::FillY8(ImgReadBuf *rb)
{
	if ((rb->cs.dtype != IDTubyte) | (rb->cs.format != IPFy))
		return false;
	const int	npix = ((rb->r.xright - rb->r.xleft)/rb->subsample) *
				((rb->r.ybottom - rb->r.ytop)/rb->subsample);
	if (!tb.IncreaseSize(npix, true))
		return false;
	if (!FillFloat(tb.fbuf, rb))
		return false;
	if (tmCvGrays(ts, tb.lbuf, tb.fbuf, npix) != TM_E_OK)
		goto tmerr;
	if (tmMapPixels(ts, rb->buf, tb.lbuf, TM_NOCHROM, npix) != TM_E_OK)
		goto tmerr;
	return true;
tmerr:
	strcpy(errMsg, tmErrorMessage[ts->lastError]);
	errCode = IREunknown;
	return false;
}

// Set up tone-mapping conversion if needed
bool
ExrReader::SetConversion(const ImgColorSpace *ocs)
{
	if (PmatchColorSpace(ocs, &cs, PICMall))
		return true;		// no conversion needed
					// else see what we can do...
	if (ocs->dtype != IDTubyte)
		return false;
	if (cs.dtype != IDTfloat)	// don't handle anything else
		return false;
	if (ocs->logorig > 0)		// can't handle log output
		return false;
	if (!((ocs->format == IPFrgb) & (cs.format == IPFrgb)) &&
			!((ocs->format == IPFy) & (cs.format == IPFy)))
		return false;
					// check/set tone-mapping
	if (ts != NULL && PRIMEQ(monprims, ocs->chroma))
		return (tb.lbuf != NULL);
					// set output space
	memcpy(monprims, ocs->chroma, sizeof(RGBPRIMS));
	int	tmfl = TM_F_NOSTDERR;
	tmfl |= (humanVis ? TM_F_HUMAN : TM_F_CAMERA);
	if (cs.format == IPFy)
		tmfl |= TM_F_BW;
	if (ts != NULL)			// redo tone-mapping
		tmDone(ts);
	ts = tmInit(tmfl, monprims, ocs->gamma);
	if (ts == NULL)
		return false;
	Header::ConstIterator	ha;	// set input space
	double			stonits;
	int			ypos, ystep, yend;
	if ((ha = exr.header().find("sampToNits")) != exr.header().end() &&
			ha.attribute().typeName() == FloatAttribute::staticTypeName())
		stonits = double(static_cast<const FloatAttribute&>
					(ha.attribute()).value());
	else
		stonits = 100.;		// reasonable assumption?
	if (tmSetSpace(ts, (RGBPRIMP)cs.chroma, stonits) != TM_E_OK)
		goto tmerr;
	nr = fr;			// compute histogram
	if (!tb.IncreaseSize(xres, cs.format==IPFy))
		goto tmerr;
	ystep = 1;
	if (xres*yres > EXR_MAXHISTO)
		ystep = xres*yres/EXR_MAXHISTO;
	if (fr.ytop) {			// DECREASING_Y
		ypos = yres/ystep*ystep - 1;
		yend = -1;
		ystep = -ystep;
	} else {			// INCREASING_Y
		ypos = 0;
		yend = yres/ystep*ystep;
	}
	while (ypos != yend) {		// read the image
		int	rval;
		if (!FloatSpan(tb.fbuf, 0, ypos, xres))
			goto tmerr;
		if (tmfl & TM_F_BW)
			rval = tmCvGrays(ts, tb.lbuf, tb.fbuf, xres);
		else
			rval = tmCvColors(ts, tb.lbuf, TM_NOCHROM, (COLOR *)tb.fbuf, xres);
		if (rval == TM_E_OK)
			rval = tmAddHisto(ts, tb.lbuf, xres, 1);
		if (rval != TM_E_OK)
			goto tmerr;
		ypos += ystep;
	}
	if (tmComputeMapping(ts, ocs->gamma, LdDyn, LdMax) == TM_E_OK) {
		switch (exr.channels()) {
		case WRITE_RGB:
		case WRITE_RGBA:
			fillF = &ExrReader::FillRGB24;
			break;
		case WRITE_R:
		case WRITE_G:
		case WRITE_B:
		case WRITE_Y:
			fillF = &ExrReader::FillY8;
			break;
		default:
			goto tmerr;	// how did we get here??
		}
		return true;		// ready for action
	}
tmerr:
	tb.Free();			// avoids repeat business
	fillF = (cs.format==IPFy ? &ExrReader::FillY32 : &ExrReader::FillRGB96);
	return false;
}

// Read a rectangle from EXR image file
bool
ExrReader::ReadRec(ImgReadBuf *rb)
{
	if (rb == NULL || rb->subsample <= 0)
		return false;
	if (!PlegalRect(&rb->r, xres, yres)) {
		strcpy(errMsg, "bogus requested image rectangle");
		errCode = IREunknown;
		return false;
	}
					// can we convert to requested CS?
	if (!SetConversion(&rb->cs)) {
		rb->cs = cs;		// make caller work
		rb->buf = NULL;
	}
	if (striplen > rb->r.ybottom - rb->r.ytop)
		rb->buf = NULL;
	if (rb->buf == NULL) {		// allocate convenient buffer
		rb->r.xleft = 0;
		rb->r.xright = xres - (xres % rb->subsample);
		rb->r.ybottom = rb->r.ytop + (striplen > rb->subsample ? striplen : rb->subsample);
		int	buflen = ImgReadBufLen(rb);
		if (buflen <= 0) {
			strcpy(errMsg, "internal buffer size error");
			errCode = IREunknown;
			return false;
		}
		rb->buf = (uby8 *)malloc(buflen);
		if (rb->buf == NULL) {
			strcpy(errMsg, "cannot allocate new read buffer");
			errCode = IREmemory;
			return false;
		}
	}
					// fill routine does actual work
	return (this->*fillF)(rb);
}
