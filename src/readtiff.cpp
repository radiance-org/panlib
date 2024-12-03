/*
 *  readtiff.cpp
 *  panlib
 *
 *  Class to read TIFF files in Pancine.
 *  Depends on Sam Leffler's TIFF library.
 *
 *  Created by gward on Fri Jun 01 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#include <iostream>
using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "system.h"
#include "imgreader.h"
#include "tiffio.h"
#include "tiffmsg.h"
#include "tiffin.h"
#include "exif.h"

// Tone-mapping data conversion routines
#include "color.h"
#include "tonemap.h"
#include "tmaptiff.h"

// Floating point comparison of primary chromaticities
#define	FEQ(a,b)	((a) < (b)+1e-5 && (b) < (a)+1e-5)
#define	PRIMEQ(p1,p2)	(FEQ((p1)[0][0],(p2)[0][0])&&FEQ((p1)[0][1],(p2)[0][1])\
			&&FEQ((p1)[1][0],(p2)[1][0])&&FEQ((p1)[1][1],(p2)[1][1])\
			&&FEQ((p1)[2][0],(p2)[2][0])&&FEQ((p1)[2][1],(p2)[2][1])\
			&&FEQ((p1)[3][0],(p2)[3][0])&&FEQ((p1)[3][1],(p2)[3][1]))

// Conversion buffer for tone-mapping operations
struct MyTMbuf {
	int			buflen;		// number of pixels allocated
	size_t			pixsiz;		// input bytes/pixel
	TMbright *		lbuf;		// luminance buffer for TM
	uby8 *			rawbuf;		// raw input data buffer
				MyTMbuf() {
					lbuf=NULL; rawbuf=NULL;
					buflen=0; pixsiz=0;
				}
				~MyTMbuf() {
					Free();
				}
	bool			IncreaseSize(int npix, size_t bpp) {
					if ((npix <= buflen &
							npix*bpp <= buflen*pixsiz))
						return true;
					Free();
					lbuf = new TMbright [npix];
					rawbuf = new uby8 [npix*bpp];
					buflen = npix;
					pixsiz = bpp;
					return true;
				}
	void			Free() {
					delete [] lbuf; lbuf = NULL;
					delete [] rawbuf; rawbuf = NULL;
					buflen = 0; pixsiz = 0;
				}
};

// Derive TIFF reader class from ImgReader base struct
struct TIFFReader : ImgReader {
private:
	char			encdesc[64];	// encoding description
	TIFF *			tif;		// TIFF library pointer
	uint16			samp_pixel;	// TIFF samples per pixel
	uint16			bits_samp;	// TIFF bits per sample
	uint16			comp;		// TIFF compression scheme
	uint32			sxr, syr;	// TIFF strip/tile size
	TIFFin			tin;		// Exif input object
	tdata_t			ibuf;		// strip/tile input buffer
	RGBPRIMS		monprims;	// target display primaries
	double			LdDyn;		// dynamic range of display
	double			LdMax;		// maximum display luminance
	bool			humanVis;	// match human visibility?
	TMstruct *		ts;		// tone-mapping structure
	MyTMbuf			tb;		// tone-mapping buffers
	float			stonits;	// sample-to-nits conversion
	void			SetPosition(int x, int y);
	bool			ReadStrip();
	bool			ReadStripSeparate();
	bool			ReadTile();
	bool			ReadTileSeparate();
	bool			(TIFFReader::*readF)();
	bool			FillCopy(ImgReadBuf *rb);
	bool			FillRGB24fromRGB96(ImgReadBuf *rb);
	bool			FillY8fromGry32(ImgReadBuf *rb);
	bool			FillRGB24fromRGB48(ImgReadBuf *rb);
	bool			FillY8fromGry16(ImgReadBuf *rb);
	bool			FillRGB24fromLuv24(ImgReadBuf *rb);
	bool			FillRGB24fromLuv32(ImgReadBuf *rb);
	bool			FillY8fromL16(ImgReadBuf *rb);
	bool			(TIFFReader::*fillF)(ImgReadBuf *rb);
	bool			Init();
	float			GetS2Nits();
	bool			SetLogConversion(const ImgColorSpace *ocs);
	bool			SetShortConversion(const ImgColorSpace *ocs);
	bool			SetConversion(const ImgColorSpace *ocs);
public:
				TIFFReader(const char *fname = NULL);
				~TIFFReader() {
					Close();
				}
	bool			Open(const char *fname);
	bool			SeekFrame(int offs, ImgSeekMode sm);
	bool			IsOpen() const {
					return (tif != NULL);
				}
	bool			GetInfo(ImgInfo *info);
	bool			GetThumbnail(ImgStruct *thm);
	void			ToneMapping(double dyn, double max, int hv);
	bool			ReadRec(ImgReadBuf *rb);
	void			Close();
};

// Open TIFF image and return reader object
static ImgReader *
TRopen(const char *fname)
{
	if (fname == NULL)
		return NULL;
	TIFFReader *	tifr = new TIFFReader(fname);
	if (!tifr->IsOpen()) {
		delete tifr;
		return NULL;		// not a TIFF
	}
	return tifr;
}

// Open next TIFF layer
static ImgReadErr
TRseekFrame(ImgReader *ir, int offs, ImgSeekMode sm)
{
	ir->errCode = IREnone;
	(*(TIFFReader *)ir).SeekFrame(offs, sm);
	return ir->errCode;
}

// Close TIFF reader and free object
static void
TRclose(ImgReader *ir)
{
	delete (TIFFReader *)ir;
}

// Get TIFF information
static ImgReadErr
TRgetInfo(ImgReader *ir, ImgInfo *info)
{
	ir->errCode = IREnone;
	if (!(*(TIFFReader *)ir).GetInfo(info) && !ir->errCode)
		ir->errCode = IREunknown;
	return ir->errCode;
}

// Get TIFF thumbnail
static ImgReadErr
TRgetThumbnail(ImgReader *ir, ImgStruct *thm)
{
	ir->errCode = IREnone;
	if (!(*(TIFFReader *)ir).GetThumbnail(thm)) {
		strcpy(ir->errMsg, "no thumbnail");
		ir->errCode = IREunsupported;
	}
	return ir->errCode;
}

// Set TIFF tone-mapping
static void
TRtoneMapping(ImgReader *ir, double dyn, double max, int hv)
{
	(*(TIFFReader *)ir).ToneMapping(dyn, max, hv);
}

// Read TIFF rectangle
static ImgReadErr
TRreadRec(ImgReader *ir, ImgReadBuf *rb)
{
	ir->errCode = IREnone;
	if (!(*(TIFFReader *)ir).ReadRec(rb) && ir->errCode == IREnone) {
		strcpy(ir->errMsg, "TIFF reader routine failure");
		ir->errCode = IREunknown;
	}
	return ir->errCode;
}

// TIFF reader interface
extern "C" {
extern const ImgReaderInterface IRInterfaceTIFF;
const ImgReaderInterface IRInterfaceTIFF = {
	"TIFF.tif",
	&TRopen, &TRseekFrame, &TRgetInfo, &TRgetThumbnail,
	&TRtoneMapping, &TRreadRec, &TRclose
};
}

// Constructor for TIFF reader
TIFFReader::TIFFReader(const char *fname)
{
	ri = &IRInterfaceTIFF;
	tif = NULL; ibuf = NULL;
	LdDyn = 100.; LdMax = 100.;
	humanVis = false;
	ts = NULL;
	errCode = IREnone;
	errMsg[0] = '\0';
	stonits = -1.f;
	Open(fname);
}

// Close TIFF input and free resources
void
TIFFReader::Close()
{
	tin.Close();
	if (tif != NULL) {
		TIFFClose(tif); tif = NULL;
	}
	if (ibuf != NULL) {
		_TIFFfree(ibuf); ibuf = NULL;
	}
	if (ts != NULL) {
		tmDone(ts); ts = NULL;
	}
	tb.Free();
}

// Open TIFF input file, returning true if all is well
bool
TIFFReader::Open(const char *fname)
{
	Close();			// close old file first
	if (fname == NULL)
		return false;
	// replace TIFF stderr reporting with buffer writers
 	TIFFSetErrorHandler(BufferTiffError);
	TIFFSetWarningHandler(BufferTiffWarning);
	TiffMessageBuffer[0] = '\0';
	errCode = IREnone;
					// call TIFF library open
	tif = TIFFOpen(fname, "r");
	if (tif == NULL)
		return false;		// cannot open or wrong file type
	strcpy(file, fname);		// save file name
	frame = 0; nframes = 0;
	frameType = IRFlayer; frameRate = 0;
	return Init();			// finish initialization
}

// Advance TIFF directory and initialize corresponding image
bool
TIFFReader::SeekFrame(int offs, ImgSeekMode sm)
{
	if (!IsOpen()) {
		frame = -1;
		return false;
	}
	int	fn = -1;
	switch (sm) {
	case IRSabs:
		fn = offs; break;
	case IRSrel:
	case IRSadv:
		fn = frame+offs; break;
	}
	if (fn == frame) {
		SetPosition(0, 0);
		return true;
	}
	if (!TIFFSetDirectory(tif, fn)) {
		if (fn == frame+1)
			nframes = fn;
		errCode = IREtruncated;
		strcpy(errMsg, "no more layers");
		frame = -1;
		return false;
	}
	return Init();
}

// Initialize image corresponding to current TIFF directory
bool
TIFFReader::Init()
{
	uint32	iwidth, iheight;	// get image dimensions
	if (!TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &iwidth) ||
			!TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &iheight) ||
			(iwidth == 0 | iheight == 0)) {
		strcpy(errMsg, "missing TIFF image width/length");
		errCode = IREformat;
		return false;		// this is hopeless
	}
	xres = iwidth; yres = iheight;	// orientation reported in GetInfo()
					// get pixel aspect ratio
	float	hdens, vdens;
	if (TIFFGetField(tif, TIFFTAG_XRESOLUTION, &hdens) &&
			TIFFGetField(tif, TIFFTAG_YRESOLUTION, &vdens))
		pixAspect = vdens / hdens;
	else
		pixAspect = 1.f;
					// initialize color space
	PcopyCS(&cs, &ICS_sRGB);
	encoding = NULL;
	uint16	phot, pconf;
	TIFFGetFieldDefaulted(tif, TIFFTAG_PHOTOMETRIC, &phot);
	TIFFGetFieldDefaulted(tif, TIFFTAG_PLANARCONFIG, &pconf);
	TIFFGetFieldDefaulted(tif, TIFFTAG_COMPRESSION, &comp);
	switch (phot) {
	case PHOTOMETRIC_RGB:
		break;
	case PHOTOMETRIC_YCBCR:
		if (comp == COMPRESSION_JPEG && pconf == PLANARCONFIG_CONTIG) {
			TIFFSetField(tif, TIFFTAG_JPEGCOLORMODE,
				JPEGCOLORMODE_RGB);
		} else {
			strcpy(errMsg, "cannot convert YCbCr color space");
			errCode = IREunsupported;
			return false;
		}
		break;
	case PHOTOMETRIC_MINISBLACK:
		PcopyCS(&cs, &ICS_Y8);
		break;
	case PHOTOMETRIC_LOGL:
		if (comp != COMPRESSION_SGILOG) {
			strcpy(errMsg, "LogL format must be SGILOG compressed");
			errCode = IREunsupported;
			return false;
		}
		TIFFSetField(tif, TIFFTAG_SGILOGDATAFMT, SGILOGDATAFMT_FLOAT);
		PcopyCS(&cs, &ICS_Y);
		break;
	case PHOTOMETRIC_LOGLUV:
		if ((comp != COMPRESSION_SGILOG & comp != COMPRESSION_SGILOG24)) {
			strcpy(errMsg, "LogLuv format must be SGILOG compressed");
			errCode = IREunsupported;
			return false;
		}
		if (pconf != PLANARCONFIG_CONTIG) {
			strcpy(errMsg, "LogLuv format must be contiguous");
			errCode = IREunsupported;
			return false;
		}
		TIFFSetField(tif, TIFFTAG_SGILOGDATAFMT, SGILOGDATAFMT_FLOAT);
		PcopyCS(&cs, &ICS_XYZ);
		break;
	default:
		sprintf(errMsg, "unsupported photometric interpretation (%d)",
				phot);
		errCode = IREunsupported;
		return false;
	}
	uint16	samp_fmt;
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT, &samp_fmt);
	TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &bits_samp);
	switch (samp_fmt) {
	case SAMPLEFORMAT_UINT:
		if (bits_samp == 8) {
			cs.dtype = IDTubyte;
		} else if (bits_samp == 16) {
			const char	*softver;
			cs.dtype = IDTushort;
			if (TIFFGetField(tif, TIFFTAG_SOFTWARE, &softver) &&
					strstr(softver, "dcraw") != NULL)
				cs.gamma = 1.0;	// hack for dcraw linear
			else
				cs.gamma = 2.2;	// XXX this is a BIG assumption!
		} else {
			sprintf(errMsg, "unsupported bits/sample (%d)",
					bits_samp);
			errCode = IREunsupported;
			return false;
		}
		break;
	case SAMPLEFORMAT_IEEEFP:
		if (bits_samp != 8*sizeof(float)) {
			sprintf(errMsg, "unsupported float size (%d bits)",
					bits_samp);
			errCode = IREunsupported;
			return false;
		}
		PrealCS(&cs);
		break;
	default:
		sprintf(errMsg, "unsupported data format (%d)", samp_fmt);
		errCode = IREunsupported;
		return false;
	}
	TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &samp_pixel);
	if (samp_pixel > ImgPixelLen[cs.format] && pconf == PLANARCONFIG_SEPARATE)
		samp_pixel = ImgPixelLen[cs.format];
	else if (samp_pixel != ImgPixelLen[cs.format]) {
		strcpy(errMsg, "inconsistent samples/pixel (alpha?)");
		errCode = IREunsupported;
		return false;
	}
	encoding = encdesc;			// set encoding description
	switch (comp) {
	case COMPRESSION_NONE:
		if (samp_fmt == SAMPLEFORMAT_IEEEFP)
			strcpy(encdesc, "real");
		else
			sprintf(encdesc, "%d-bit", samp_pixel*bits_samp);
		break;
	case COMPRESSION_LZW:
		strcpy(encdesc, "LZW");
		break;
	case COMPRESSION_JPEG:
		strcpy(encdesc, "DCT");
		break;
	case COMPRESSION_PACKBITS:
	case COMPRESSION_THUNDERSCAN:
	case COMPRESSION_IT8LW:
		strcpy(encdesc, "RLE");
		break;
	case COMPRESSION_PIXARFILM:
		strcpy(encdesc, "LZW-Log");
		break;
	case COMPRESSION_PIXARLOG:
		strcpy(encdesc, "Zip-Log");
		break;
	case COMPRESSION_DEFLATE:
	case COMPRESSION_ADOBE_DEFLATE:
		strcpy(encdesc, "Deflated");
		break;
	case COMPRESSION_SGILOG:
		strcpy(encdesc, "RLE-Log");
		break;
	case COMPRESSION_SGILOG24:
		strcpy(encdesc, "24-bit-Log");
		break;
	default:
		strcpy(encdesc, "compressed");
		break;
	}
	switch (phot) {
	case PHOTOMETRIC_RGB:
		strcat(encdesc, " RGB");
		break;
	case PHOTOMETRIC_YCBCR:
		strcat(encdesc, " YCbCr");
		break;
	case PHOTOMETRIC_MINISBLACK:
		strcat(encdesc, " grayscale");
		break;
	case PHOTOMETRIC_LOGL:
		strcat(encdesc, " luminance");
		break;
	case PHOTOMETRIC_LOGLUV:
		strcat(encdesc, " Luv");
		break;
	}
	float *	fa;	// XXX known bug in libtiff -- only positive chroma
	if (TIFFGetField(tif, TIFFTAG_PRIMARYCHROMATICITIES, &fa)) {
		cs.chroma[0][0] = fa[0]; cs.chroma[0][1] = fa[1];
		cs.chroma[1][0] = fa[2]; cs.chroma[1][1] = fa[3];
		cs.chroma[2][0] = fa[4]; cs.chroma[2][1] = fa[5];
	}
	if (TIFFGetField(tif, TIFFTAG_WHITEPOINT, &fa)) {
		cs.chroma[3][0] = fa[0]; cs.chroma[3][1] = fa[1];
	}
	// XXX should get gamma from transfer function, but dunno how...
	fillF = &TIFFReader::FillCopy;
					// initialize scanning
	tsize_t	iblen;
	if (TIFFGetField(tif, TIFFTAG_TILEWIDTH, &sxr) &&
			TIFFGetField(tif, TIFFTAG_TILELENGTH, &syr)) {
		if ((int)sxr > xres) {
			strcpy(errMsg, "tiles wider than image");
			errCode = IREunsupported;
			return false;
		}
		if ((int)syr > yres) {
			strcpy(errMsg, "tiles taller than image");
			errCode = IREunsupported;
			return false;
		}
		iblen = TIFFTileSize(tif);
		if (pconf == PLANARCONFIG_SEPARATE && samp_pixel > 1) {
			iblen *= samp_pixel+1;
			readF = &TIFFReader::ReadTileSeparate;
		} else
			readF = &TIFFReader::ReadTile;
	} else {
		TIFFGetFieldDefaulted(tif, TIFFTAG_ROWSPERSTRIP, &syr);
		if ((int)syr > yres)
			syr = yres;
		sxr = xres;
		iblen = TIFFStripSize(tif);
		if (pconf == PLANARCONFIG_SEPARATE && samp_pixel > 1) {
			iblen *= samp_pixel+1;
			readF = &TIFFReader::ReadStripSeparate;
		} else
			readF = &TIFFReader::ReadStrip;
	}
	if (ibuf == NULL)
		ibuf = _TIFFmalloc(iblen);
	else
		ibuf = _TIFFrealloc(ibuf, iblen);
	if (ibuf == NULL) {
		strcpy(errMsg, "cannot allocate new input buffer");
		errCode = IREmemory;
		return false;
	}
	fr.xleft = 0; fr.xright = sxr;
	fr.ytop = 0; fr.ybottom = syr;
	nr = fr;
	return true;			// ready to roll
}

// Reposition input and prepare to read next strip/tile
void
TIFFReader::SetPosition(int x, int y)
{
	nr.xleft = x - (x % sxr);
	nr.xright = nr.xleft + sxr;
	nr.ytop = y - (y % syr);
	nr.ybottom = nr.ytop + syr;
}

// Fill read buffer without sample translation
bool
TIFFReader::FillCopy(ImgReadBuf *rb)
{
	const int	bytes_pixel = samp_pixel*(bits_samp >> 3);
					// position input
	SetPosition(rb->r.xleft, rb->r.ytop);
	while (nr.ytop < rb->r.ybottom) {
		ImgRect	ibr = nr;	// read next rectangle
		int	y;
		if (!(this->*readF)())
			return false;
					// copy samples
		if (rb->subsample > 1) {
			y = ibr.ytop + rb->subsample - 1;
			y -= y % rb->subsample;
			if (y < rb->r.ytop) y = rb->r.ytop;
			for ( ; y < ibr.ybottom; y += rb->subsample) {
				if (y >= rb->r.ybottom)
					break;
				int	x = ibr.xleft + rb->subsample - 1;
				x -= x % rb->subsample;
				if (x < rb->r.xleft) x = rb->r.xleft;
				uby8 *	pd = rb->buf + bytes_pixel *
					( (y - rb->r.ytop)/rb->subsample *
				((rb->r.xright - rb->r.xleft)/rb->subsample)
					+ (x - rb->r.xleft)/rb->subsample );
				uby8 *	ps = (uby8 *)ibuf + bytes_pixel *
					((y - ibr.ytop)*sxr + (x - ibr.xleft));
				for ( ; x < ibr.xright; x += rb->subsample) {
					if (x >= rb->r.xright)
						break;
					memcpy(pd, ps, bytes_pixel);
					pd += bytes_pixel;
					ps += rb->subsample * bytes_pixel;
				}
			}
		} else {		// no subsampling -- faster
			y = ibr.ytop;
			if (y < rb->r.ytop) y = rb->r.ytop;
			for ( ; y < ibr.ybottom; y++) {
				if (y >= rb->r.ybottom)
					break;
				int	x = ibr.xleft;
				if (x < rb->r.xleft) x = rb->r.xleft;
				uby8 *	pd = rb->buf + bytes_pixel *
					( (y - rb->r.ytop) *
						(rb->r.xright - rb->r.xleft)
					+ (x - rb->r.xleft) );
				uby8 *	ps = (uby8 *)ibuf + bytes_pixel *
					((y - ibr.ytop)*sxr + (x - ibr.xleft));
				int	nbytes = bytes_pixel*(ibr.xright - x);
				if (rb->r.xright < ibr.xright)
					nbytes = bytes_pixel*(rb->r.xright - x);
				memcpy(pd, ps, nbytes);
			}
		}
		if ((int)sxr < xres)	// skip unneeded tiles
			if (nr.xleft > rb->r.xright)
				SetPosition(rb->r.xleft, nr.ybottom);
			else if (nr.xright <= rb->r.xleft)
				SetPosition(rb->r.xleft, nr.ytop);
	}
	return true;
}

// Fill rectangle, tone-mapping from RGB float to RGB24
bool
TIFFReader::FillRGB24fromRGB96(ImgReadBuf *rb)
{
	int	npixels = (rb->r.xright - rb->r.xleft)/rb->subsample *
			(rb->r.ybottom - rb->r.ytop)/rb->subsample;
	if (!tb.IncreaseSize(npixels, 3*sizeof(float))) {
		strcpy(errMsg, "cannot allocate tone-mapping buffers");
		errCode = IREmemory;
		return false;
	}
	ImgReadBuf	tempb;		// set up copy buffer
	tempb.r = rb->r;		// use their rectangle
	tempb.subsample = rb->subsample;
	tempb.cs = cs;			// pretend our color space
	tempb.buf = tb.rawbuf;
	if (!FillCopy(&tempb))		// fill raw buffer
		return false;
					// convert pixels
	int	tmr;
	tmr = tmCvColors(ts, tb.lbuf, rb->buf, (COLOR *)tb.rawbuf, npixels);
	if (tmr != TM_E_OK)
		goto tmerr;
					// apply tone-mapping
	tmr = tmMapPixels(ts, rb->buf, tb.lbuf, rb->buf, npixels);
	if (tmr == TM_E_OK)
		return true;
tmerr:
	strcpy(errMsg, tmErrorMessage[tmr]);
	errCode = IREunknown;
	return false;
}

// Fill rectangle, tone-mapping from Gray float to Y8
bool
TIFFReader::FillY8fromGry32(ImgReadBuf *rb)
{
	int	npixels = (rb->r.xright - rb->r.xleft)/rb->subsample *
			(rb->r.ybottom - rb->r.ytop)/rb->subsample;
	if (!tb.IncreaseSize(npixels, sizeof(float))) {
		strcpy(errMsg, "cannot allocate tone-mapping buffers");
		errCode = IREmemory;
		return false;
	}
	ImgReadBuf	tempb;		// set up copy buffer
	tempb.r = rb->r;		// use their rectangle
	tempb.subsample = rb->subsample;
	tempb.cs = cs;			// pretend our color space
	tempb.buf = tb.rawbuf;
	if (!FillCopy(&tempb))		// fill raw buffer
		return false;
					// convert pixels
	int	tmr;
	tmr = tmCvGrays(ts, tb.lbuf, (float *)tb.rawbuf, npixels);
	if (tmr != TM_E_OK)
		goto tmerr;
					// apply tone-mapping
	tmr = tmMapPixels(ts, (uby8 *)rb->buf, tb.lbuf, TM_NOCHROM, npixels);
	if (tmr == TM_E_OK)
		return true;
tmerr:
	strcpy(errMsg, tmErrorMessage[tmr]);
	errCode = IREunknown;
	return false;
}

// Fill rectangle, tone-mapping from RGB unsigned short to RGB24
bool
TIFFReader::FillRGB24fromRGB48(ImgReadBuf *rb)
{
	int	npixels = (rb->r.xright - rb->r.xleft)/rb->subsample *
			(rb->r.ybottom - rb->r.ytop)/rb->subsample;
	if (!tb.IncreaseSize(npixels, 3*sizeof(uint16))) {
		strcpy(errMsg, "cannot allocate tone-mapping buffers");
		errCode = IREmemory;
		return false;
	}
	ImgReadBuf	tempb;		// set up copy buffer
	tempb.r = rb->r;		// use their rectangle
	tempb.subsample = rb->subsample;
	tempb.cs = cs;			// pretend our color space
	tempb.buf = tb.rawbuf;
	if (!FillCopy(&tempb))		// fill raw buffer
		return false;
					// convert pixels
	int	tmr;
	tmr = tmCvRGB48(ts, tb.lbuf, rb->buf, (uint16 (*)[3])tb.rawbuf,
				npixels, cs.gamma);
	if (tmr != TM_E_OK)
		goto tmerr;
					// apply tone-mapping
	tmr = tmMapPixels(ts, rb->buf, tb.lbuf, rb->buf, npixels);
	if (tmr == TM_E_OK)
		return true;
tmerr:
	strcpy(errMsg, tmErrorMessage[tmr]);
	errCode = IREunknown;
	return false;
}

// Fill rectangle, tone-mapping from Gray unsigned short to Y8
bool
TIFFReader::FillY8fromGry16(ImgReadBuf *rb)
{
	int	npixels = (rb->r.xright - rb->r.xleft)/rb->subsample *
			(rb->r.ybottom - rb->r.ytop)/rb->subsample;
	if (!tb.IncreaseSize(npixels, sizeof(uint16))) {
		strcpy(errMsg, "cannot allocate tone-mapping buffers");
		errCode = IREmemory;
		return false;
	}
	ImgReadBuf	tempb;		// set up copy buffer
	tempb.r = rb->r;		// use their rectangle
	tempb.subsample = rb->subsample;
	tempb.cs = cs;			// pretend our color space
	tempb.buf = tb.rawbuf;
	if (!FillCopy(&tempb))		// fill raw buffer
		return false;
					// convert pixels
	int	tmr;
	tmr = tmCvGray16(ts, tb.lbuf, (uint16 *)tb.rawbuf, npixels, cs.gamma);
	if (tmr != TM_E_OK)
		goto tmerr;
					// apply tone-mapping
	tmr = tmMapPixels(ts, (uby8 *)rb->buf, tb.lbuf, TM_NOCHROM, npixels);
	if (tmr == TM_E_OK)
		return true;
tmerr:
	strcpy(errMsg, tmErrorMessage[tmr]);
	errCode = IREunknown;
	return false;
}

// Fill rectangle, tone-mapping from LogLuv24 to RGB24
bool
TIFFReader::FillRGB24fromLuv24(ImgReadBuf *rb)
{
	int	npixels = (rb->r.xright - rb->r.xleft)/rb->subsample *
			(rb->r.ybottom - rb->r.ytop)/rb->subsample;
	if (!tb.IncreaseSize(npixels, sizeof(uint32))) {
		strcpy(errMsg, "cannot allocate tone-mapping buffers");
		errCode = IREmemory;
		return false;
	}
	ImgReadBuf	tempb;		// set up copy buffer
	tempb.r = rb->r;		// use their rectangle
	tempb.subsample = rb->subsample;
	tempb.cs = cs;			// pretend our color space
	tempb.buf = tb.rawbuf;
	if (!FillCopy(&tempb))		// fill raw buffer
		return false;
					// convert pixels
	int	tmr;
	tmr = tmCvLuv24(ts, tb.lbuf, rb->buf, (uint32 *)tb.rawbuf, npixels);
	if (tmr != TM_E_OK)
		goto tmerr;
					// apply tone-mapping
	tmr = tmMapPixels(ts, rb->buf, tb.lbuf, rb->buf, npixels);
	if (tmr == TM_E_OK)
		return true;
tmerr:
	strcpy(errMsg, tmErrorMessage[tmr]);
	errCode = IREunknown;
	return false;
}

// Fill rectangle, tone-mapping from LogLuv32 to RGB24
bool
TIFFReader::FillRGB24fromLuv32(ImgReadBuf *rb)
{
	int	npixels = (rb->r.xright - rb->r.xleft)/rb->subsample *
			(rb->r.ybottom - rb->r.ytop)/rb->subsample;
	if (!tb.IncreaseSize(npixels, sizeof(uint32))) {
		strcpy(errMsg, "cannot allocate tone-mapping buffers");
		errCode = IREmemory;
		return false;
	}
	ImgReadBuf	tempb;		// set up copy buffer
	tempb.r = rb->r;		// use their rectangle
	tempb.subsample = rb->subsample;
	tempb.cs = cs;			// pretend our color space
	tempb.buf = tb.rawbuf;
	if (!FillCopy(&tempb))		// fill raw buffer
		return false;
					// convert pixels
	int	tmr;
	tmr = tmCvLuv32(ts, tb.lbuf, rb->buf, (uint32 *)tb.rawbuf, npixels);
	if (tmr != TM_E_OK)
		goto tmerr;
					// apply tone-mapping
	tmr = tmMapPixels(ts, rb->buf, tb.lbuf, rb->buf, npixels);
	if (tmr == TM_E_OK)
		return true;
tmerr:
	strcpy(errMsg, tmErrorMessage[tmr]);
	errCode = IREunknown;
	return false;
}

// Fill rectangle, tone-mapping from LogL16 to Y8
bool
TIFFReader::FillY8fromL16(ImgReadBuf *rb)
{
	int	npixels = (rb->r.xright - rb->r.xleft)/rb->subsample *
			(rb->r.ybottom - rb->r.ytop)/rb->subsample;
	if (!tb.IncreaseSize(npixels, sizeof(uint16))) {
		strcpy(errMsg, "cannot allocate tone-mapping buffers");
		errCode = IREmemory;
		return false;
	}
	ImgReadBuf	tempb;		// set up copy buffer
	tempb.r = rb->r;		// use their rectangle
	tempb.subsample = rb->subsample;
	tempb.cs = cs;			// pretend our color space
	tempb.buf = tb.rawbuf;
	if (!FillCopy(&tempb))		// fill raw buffer
		return false;
					// convert pixels
	int	tmr;
	tmr = tmCvL16(ts, tb.lbuf, (uint16 *)tb.rawbuf, npixels);
	if (tmr != TM_E_OK)
		goto tmerr;
					// apply tone-mapping
	tmr = tmMapPixels(ts, (uby8 *)rb->buf, tb.lbuf, TM_NOCHROM, npixels);
	if (tmr == TM_E_OK)
		return true;
tmerr:
	strcpy(errMsg, tmErrorMessage[tmr]);
	errCode = IREunknown;
	return false;
}

// Read next strip of interleaved pixel samples into preallocated buffer
bool
TIFFReader::ReadStrip()
{
	if (nr.ytop >= yres)
		return false;		// at EOD

	if (TIFFReadEncodedStrip(tif, (tstrip_t)(nr.ytop/syr),
						ibuf, (tsize_t)(-1)) <= 0) {
		strcpy(errMsg, TiffMessageBuffer);
		errCode = IREread;
		return false;
	}
	nr.ytop = nr.ybottom;		// advance to next strip
	if ((nr.ybottom += syr) > yres)
		nr.ybottom = yres;
	return true;
}

// Read next strip of separated pixel samples into preallocated buffer
bool
TIFFReader::ReadStripSeparate()
{
	if (nr.ytop >= yres)
		return false;		// atEOD
	/*
	 *  Read each strip into extra space at end of input buffer,
	 *  then reswizzle them into front of buffer.
	 */
	const tstrip_t	strip_samp = (yres + (syr-1)) / syr;
	const int	bytes_samp = bits_samp >> 3;
	tdata_t		dp = (tdata_t)((uby8 *)ibuf +
					sxr*syr*samp_pixel*bytes_samp);
	tstrip_t	strip = nr.ytop / syr;
	int		s;		// read each sample's strip
	for (s = 0; s < samp_pixel; s++) {
		if (TIFFReadEncodedStrip(tif, strip, dp, (tsize_t)(-1)) <= 0) {
			strcpy(errMsg, TiffMessageBuffer);
			errCode = IREread;
			return false;
		}
		uby8 *	pd = (uby8 *)ibuf + s*bytes_samp;
		uby8 *	ps = (uby8 *)dp;
		int		n = sxr*(nr.ybottom - nr.ytop);
		if (bytes_samp == 1)	// interleave samples
			while (n--) {
				*pd = *ps++;
				pd += samp_pixel;
			}
		else
			while (n--) {
				memcpy(pd, ps, bytes_samp);
				pd += samp_pixel*bytes_samp;
				ps += bytes_samp;
			}
		strip += strip_samp;	// do next sample plane
	}
	nr.ytop = nr.ybottom;		// advance to next strip
	if ((nr.ybottom += syr) > yres)
		nr.ybottom = yres;
	return true;
}

// Read next tile of interleaved pixel samples into preallocated buffer
bool
TIFFReader::ReadTile()
{
	if (nr.ytop >= yres)
		return false;		// at EOD

	if (TIFFReadTile(tif, ibuf, nr.xleft, nr.ytop, 0, 0) <= 0) {
		strcpy(errMsg, TiffMessageBuffer);
		errCode = IREread;
		return false;
	}
	nr.xleft = nr.xright;		// advance to next column
	if ((nr.xright += sxr) > xres)
		nr.xright = xres;
	if (nr.xleft >= xres) {		// beginning of next row
		nr.xleft = 0; nr.xright = sxr;
		nr.ytop = nr.ybottom;
		if ((nr.ybottom += syr) > yres)
			nr.ybottom = yres;
	}
	return true;
}

// Read next set of separate pixel tiles into preallocated buffer
bool
TIFFReader::ReadTileSeparate()
{
	if (nr.ytop >= yres)
		return false;		// at EOD
	/*
	 *  Read each tile into extra space at end of input buffer,
	 *  then reswizzle them into front of buffer.
	 */
	const int	bytes_samp = bits_samp >> 3;
	tdata_t		dp = (tdata_t)((uby8 *)ibuf +
					sxr*syr*samp_pixel*bytes_samp);
	int		s;		// read each sample's tile
	for (s = 0; s < samp_pixel; s++) {
		if (TIFFReadTile(tif, dp, nr.xleft, nr.ytop, 0, s) <= 0) {
			strcpy(errMsg, TiffMessageBuffer);
			errCode = IREread;
			return false;
		}
		uby8 *	pd = (uby8 *)ibuf + s*bytes_samp;
		uby8 *	ps = (uby8 *)dp;
		int		n = (nr.xright - nr.xleft) *
					(nr.ybottom - nr.ytop);
		if (bytes_samp == 1)	// interleave samples
			while (n--) {
				*pd = *ps++;
				pd += samp_pixel;
			}
		else
			while (n--) {
				memcpy(pd, ps, bytes_samp);
				pd += samp_pixel*bytes_samp;
				ps += bytes_samp;
			}
	}
	nr.xleft = nr.xright;		// advance to next column
	if ((nr.xright += sxr) > xres)
		nr.xright = xres;
	if (nr.xleft >= xres) {		// beginning of next row
		nr.xleft = 0; nr.xright = sxr;
		nr.ytop = nr.ybottom;
		if ((nr.ybottom += syr) > yres)
			nr.ybottom = yres;
	}
	return true;
}

// Get Sample-to-nits conversion factor
float
TIFFReader::GetS2Nits()
{
	if (stonits <= 0) {
		ImgInfo	iinf;
		iinf.flags = 0;
		if (tin.Ready() || tin.OpenFile(file))
			GetExifInfo(&iinf, &tin);
		if (iinf.flags & IIFstonits)
			stonits = iinf.stonits;
		else
			stonits = 100.f;	// XXX assumption
	}
	return stonits;
}

// Set new tone-mapping parameters
void
TIFFReader::ToneMapping(double dyn, double max, int hv)
{
	LdDyn = dyn;			// store new parameters
	LdMax = max;
	humanVis = (bool)hv;
	if (ts != NULL) {		// reset tone-mapping
		tmDone(ts); ts = NULL;
	}
}

// Set up color conversion and tone-mapping from LogLuv
bool
TIFFReader::SetLogConversion(const ImgColorSpace *ocs)
{
	if (!((ocs->format == IPFrgb) & (cs.format == IPFxyz)) &&
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
	int	npixels;		// set input space
	if (tmSetSpace(ts, TM_XYZPRIM, GetS2Nits()) != TM_E_OK)
		goto tmerr;
	if (cs.format == IPFxyz) {	// set raw input format & data size
		TIFFSetField(tif, TIFFTAG_SGILOGDATAFMT, SGILOGDATAFMT_RAW);
		bits_samp = 32;
		samp_pixel = 1;
	} else /* cs.format == IPFy */ {
		TIFFSetField(tif, TIFFTAG_SGILOGDATAFMT, SGILOGDATAFMT_16BIT);
		bits_samp = 16;
	}
				// don't reallocate ibuf -- it's big enough
	SetPosition(0, 0);		// compute histogram
	npixels = sxr*syr;
	if (!tb.IncreaseSize(npixels, cs.format==IPFy ?
				sizeof(uint16) : sizeof(uint32)))
		goto tmerr;
	while (npixels > 0) {		// read the whole image
		if (!(this->*readF)())
			goto tmerr;
		int	rv = TM_E_TMINVAL;
		if (cs.format == IPFy)
			rv = tmCvL16(ts, tb.lbuf, (uint16 *)ibuf, npixels);
		else if (comp == COMPRESSION_SGILOG24)
			rv = tmCvLuv24(ts, tb.lbuf, TM_NOCHROM, (uint32 *)ibuf, npixels);
		else if (comp == COMPRESSION_SGILOG)
			rv = tmCvLuv32(ts, tb.lbuf, TM_NOCHROM, (uint32 *)ibuf, npixels);
		if (rv == TM_E_OK)
			rv = tmAddHisto(ts, tb.lbuf, npixels, 1);
		if (rv != TM_E_OK)
			goto tmerr;
		npixels = (nr.xright - nr.xleft)*(nr.ybottom - nr.ytop);
	}
	if (tmComputeMapping(ts, ocs->gamma, LdDyn, LdMax) == TM_E_OK) {
		if (cs.format == IPFy)
			fillF = &TIFFReader::FillY8fromL16;
		else if (comp == COMPRESSION_SGILOG24)
			fillF = &TIFFReader::FillRGB24fromLuv24;
		else /* comp == COMPRESSION_SGILOG */
			fillF = &TIFFReader::FillRGB24fromLuv32;
		return true;		// ready for action
	}
tmerr:
	tb.Free();			// avoids repeat business
	fillF = &TIFFReader::FillCopy;
	return false;
}

// Set up color conversion and tone-mapping from 16-bit samples
bool
TIFFReader::SetShortConversion(const ImgColorSpace *ocs)
{
	if (!(ocs->format == IPFrgb & cs.format == IPFrgb) &&
			!(ocs->format == IPFy & cs.format == IPFy))
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
	int	npixels;		// set input space
	if (tmSetSpace(ts, (RGBPRIMP)cs.chroma, GetS2Nits()) != TM_E_OK)
		goto tmerr;
	SetPosition(0, 0);		// compute histogram
	npixels = sxr*syr;
	if (!tb.IncreaseSize(npixels, sizeof(uint16)*(cs.format==IPFy?1:3)))
		goto tmerr;
	while (npixels > 0) {		// read the whole image
		if (!(this->*readF)())
			goto tmerr;
		int	rv = TM_E_TMINVAL;
		if (cs.format == IPFy)
			rv = tmCvGray16(ts, tb.lbuf, (uint16 *)ibuf, npixels, cs.gamma);
		else /* cs.format == IPFrgb */
			rv = tmCvRGB48(ts, tb.lbuf, TM_NOCHROM, (uint16 (*)[3])ibuf,
					npixels, cs.gamma);
		if (rv == TM_E_OK)
			rv = tmAddHisto(ts, tb.lbuf, npixels, 1);
		if (rv != TM_E_OK)
			goto tmerr;
		npixels = (nr.xright - nr.xleft)*(nr.ybottom - nr.ytop);
	}
	if (tmComputeMapping(ts, ocs->gamma, LdDyn, LdMax) == TM_E_OK) {
		if (cs.format == IPFy)
			fillF = &TIFFReader::FillY8fromGry16;
		else /* cs.format == IPFrgb */
			fillF = &TIFFReader::FillRGB24fromRGB48;
		return true;		// ready for action
	}
tmerr:
	tb.Free();			// avoids repeat business
	fillF = &TIFFReader::FillCopy;
	return false;
}

// Set up color conversion and tone-mapping if necessary
bool
TIFFReader::SetConversion(const ImgColorSpace *ocs)
{
	if (PmatchColorSpace(ocs, &cs, PICMall))
		return true;		// no conversion needed
					// else see what we can do...
	if (ocs->dtype != IDTubyte)
		return false;
	if (ocs->logorig > 0)		// don't handle direct log output
		return false;
	if ((comp == COMPRESSION_SGILOG | comp == COMPRESSION_SGILOG24))
		return SetLogConversion(ocs);
	if (cs.dtype == IDTushort)
		return SetShortConversion(ocs);
	if (cs.dtype != IDTfloat)
		return false;
	if (!(ocs->format == IPFrgb & cs.format == IPFrgb) &&
			!(ocs->format == IPFy & cs.format == IPFy))
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
	int	npixels;		// set input space
	if (tmSetSpace(ts, (RGBPRIMP)cs.chroma, GetS2Nits()) != TM_E_OK)
		goto tmerr;
	SetPosition(0, 0);		// compute histogram
	npixels = sxr*syr;
	if (!tb.IncreaseSize(npixels, sizeof(float)*(cs.format==IPFy?1:3)))
		goto tmerr;
	while (npixels > 0) {		// read the whole image
		if (!(this->*readF)())
			goto tmerr;
		int	rv = TM_E_TMINVAL;
		if (cs.format == IPFy)
			rv = tmCvGrays(ts, tb.lbuf, (float *)ibuf, npixels);
		else /* cs.format == IPFrgb */
			rv = tmCvColors(ts, tb.lbuf, TM_NOCHROM, (COLOR *)ibuf, npixels);
		if (rv == TM_E_OK)
			rv = tmAddHisto(ts, tb.lbuf, npixels, 1);
		if (rv != TM_E_OK)
			goto tmerr;
		npixels = (nr.xright - nr.xleft)*(nr.ybottom - nr.ytop);
	}
	if (tmComputeMapping(ts, ocs->gamma, LdDyn, LdMax) == TM_E_OK) {
		if (cs.format == IPFy)
			fillF = &TIFFReader::FillY8fromGry32;
		else /* cs.format == IPFrgb */
			fillF = &TIFFReader::FillRGB24fromRGB96;
		return true;		// ready for action
	}
tmerr:
	tb.Free();			// avoids repeat business
	fillF = &TIFFReader::FillCopy;
	return false;
}

// Read a rectangle from open TIFF
bool
TIFFReader::ReadRec(ImgReadBuf *rb)
{
	if (!IsOpen() || rb == NULL || rb->subsample <= 0)
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
	if (rb->buf == NULL) {		// allocate convenient buffer for us
		SetPosition(rb->r.xleft, rb->r.ytop);
		rb->r = nr;		// minimal buffer
		ImgFixSampling(rb);
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

// Extract Exif information from TIFF
bool
TIFFReader::GetInfo(ImgInfo *info)
{
	if (info == NULL)
		return false;
	*info = defImgInfo;		// clear info
	if (!IsOpen())
		return false;
	if (!tin.Ready() && !tin.OpenFile(file))
		return false;
	GetExifInfo(info, &tin);	// ignore return value
	return true;
}

// Extract thumbnail from TIFF
bool
TIFFReader::GetThumbnail(ImgStruct *thm)
{
	if (!IsOpen())
		return false;
	if (!tin.Ready() && !tin.OpenFile(file))
		return false;
	return GetExifThumbnail(thm, &tin);
}
