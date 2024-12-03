/*
 *  readjpeg.cpp
 *  panlib
 *
 *  Class to read (HDR) JPEG files in Pancine.
 *  Depends on Tome Lane's JPEG library with Sunnybrook Tech's extensions.
 *
 *  Created by gward on Thu May 24 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
using namespace std;
#include <string.h>
#include "system.h"
#include "imgreader.h"
#include "jstreamsrc.h"
#include "jpeghdr.h"
#include "tiffin.h"
#include "exif.h"

#ifndef SF_JPEG_IDCT
#define SF_JPEG_IDCT	(JPEG_LIB_VERSION==62)	// use custom DCT decoder?
#endif

#define MAXSAVELEN		0xfffd		// maximum saved marker length

enum JPState {JPSclosed, JPSopen, JPSstart, JPSscan, JPSerror};

// Derive JPEG reader class from base ImgReader struct
struct JPEGReader : ImgReader {
private:
	JPState			state;		// reader state
	istream *		istr;		// JPEG input stream
	jpeghdr_decompress_struct
				jhinf;		// HDR JPEG info struct
	jpeg_error_mgr		jerr;		// error handler
	JSAMPLE *		ldrscan;	// LDR scanline buffer
	JHSAMPLE *		hdrscan;	// HDR scanline buffer
	bool			HasHDR() {
					return (cs.dtype == IDTfloat);
				}
	bool			ReadingHDR() {
					return (hdrscan != NULL);
				}
	bool			DetermineCS();
	bool			JPreset();
	bool			JPstart(bool getHDR);
	bool			JPnextScanline(JSAMPLE *sl = NULL,
						JHSAMPLE *hsl = NULL);
	bool			JPsetRec(ImgReadBuf *rb);
	void			JumpSet(jumper_struct *js = NULL) {
					jhinf.cinfo.client_data = (void *)js;
				}
	void			DiagnoseError();
public:
				JPEGReader(const char *fname = NULL);
				~JPEGReader() {
					Close();
				}
	bool			Open(const char *fname);
	bool			IsOpen() const {
					return state != JPSclosed;
				}
	bool			GetInfo(ImgInfo *info);
	bool			GetThumbnail(ImgStruct *thm);
	bool			ReadRec(ImgReadBuf *rb);
	void			Close();
};


// Open JPEG image and return reader object
static ImgReader *
JRopen(const char *fname)
{
	if (fname == NULL)
		return NULL;
	JPEGReader *	jpr = new JPEGReader(fname);
	if (!jpr->IsOpen()) {
		delete jpr;
		return NULL;				// not a JPEG
	}
	return jpr;
}

// Close JPEG reader and free object
static void
JRclose(ImgReader *ir)
{
	delete (JPEGReader *)ir;
}

// Get JPEG information
static ImgReadErr
JRgetInfo(ImgReader *ir, ImgInfo *info)
{
	ir->errCode = IREnone;
	if (!(*(JPEGReader *)ir).GetInfo(info) && !ir->errCode)
		ir->errCode = IREunknown;
	return ir->errCode;
}

// Get JPEG thumbnail
static ImgReadErr
JRgetThumbnail(ImgReader *ir, ImgStruct *thm)
{
	ir->errCode = IREnone;
	if (!(*(JPEGReader *)ir).GetThumbnail(thm)) {
		strcpy(ir->errMsg, "no thumbnail");
		ir->errCode = IREunsupported;
	}
	return ir->errCode;
}

// Read JPEG rectangle
static ImgReadErr
JRreadRec(ImgReader *ir, ImgReadBuf *rb)
{
	ir->errCode = IREnone;
	if (!(*(JPEGReader *)ir).ReadRec(rb) && ir->errCode == IREnone) {
		strcpy(ir->errMsg, "JPEG reader routine failure");
		ir->errCode = IREunknown;
	}
	return ir->errCode;
}

// JPEG reader interface
extern "C" {
extern const ImgReaderInterface IRInterfaceJPEG;
const ImgReaderInterface IRInterfaceJPEG = {
	"JPEG.jpg.jfif.jff.jiff",
	&JRopen, NULL, &JRgetInfo, &JRgetThumbnail,
	NULL, &JRreadRec, &JRclose
};
}

// Constructor for JPEG reader
JPEGReader::JPEGReader(const char *fname)
{
	ri = &IRInterfaceJPEG;
	state = JPSclosed;
	istr = NULL;
	ldrscan = NULL;
	hdrscan = NULL;
	errCode = IREnone;
	errMsg[0] = '\0';
	Open(fname);
}

// Close our JPEG stream and free reader resources
void
JPEGReader::Close()
{
	if (ldrscan != NULL) {
		free(ldrscan);
		ldrscan = NULL;
	}
	if (hdrscan != NULL) {
		free(hdrscan);
		hdrscan = NULL;
	}
	if (state == JPSclosed)
		return;
	jpeghdr_destroy_decompress(&jhinf);
	delete istr; istr = NULL;
	state = JPSclosed;
}

// Open a JPEG file and read header, returning true on success
// Sets errCode to IREnone on basic errors (file won't open or not JPEG)
bool
JPEGReader::Open(const char *fname)
{
	if (state != JPSclosed)			// close previous file
		Close();
	if (fname == NULL)
		return false;
	errCode = IREnone;
	ifstream *	ifstr = new ifstream(fname, ios::in | ios::binary);
	if (!ifstr->is_open()) {
		delete ifstr;
		return false;			// cannot even open it!
	}
	istr = ifstr;
	strcpy(file, fname);			// save file name
	encoding = NULL;
	jhinf.cinfo.err = jpeg_std_error(&jerr);	// prepare JPEG struct
	jerr.output_message = jpeg_error_output;
	jerr.emit_message = jpeg_emit_message;
	jerr.error_exit = jpeg_error_jump;
	jumper_struct	myjumper;
	JumpSet(&myjumper);
	if (setjmp(myjumper.b)) {		// jpeg_error_jump to here
		DiagnoseError();
errorJ:
		if (jerr.msg_code == JERR_NO_SOI) {
			Close();
			errCode = IREnone;	// not even a JPEG
		}
		return false;
	}
	jpeghdr_create_decompress(&jhinf);
	state = JPSopen;			// set initial state
	JPreset();
	JumpSet();
	if (state != JPSstart)
		goto errorJ;			// bad header
	frame = 0; nframes = 1;
	frameType = IRFnone; frameRate = 0;
	return true;
}

// Read header and prepare for decompression (aborting any read in progress)
bool
JPEGReader::JPreset()
{
	if (state == JPSclosed)
		return false;
	if (state != JPSopen) {			// reset input
		jpeghdr_abort_decompress(&jhinf);
		istr->clear();
		istr->seekg(0);
		state = JPSopen;
	}
	// Exif thumbnail is in APP1 marker, along with camera parameters
	jpeg_save_markers(&jhinf.cinfo, JPEG_APP0+1, MAXSAVELEN);
	// JFIF thumbnail is in APP0 marker, along with camera parameters
	jpeg_save_markers(&jhinf.cinfo, JPEG_APP0+0, MAXSAVELEN);
	// Attach input stream as JPEG source
	jpeg_stream_src(&jhinf.cinfo, istr);
	// Check header for HDR input
	int	hrv = jpeghdr_read_header(&jhinf);
	if ((hrv != JPEG_HEADER_OK) & (hrv != JPEG_HEADER_HDR)) {
		DiagnoseError();
		return false;
	}
	state = JPSstart;			// indicates after header
						// set reader parameters
	xres = jhinf.cinfo.image_width;
	yres = jhinf.cinfo.image_height;
	if (jhinf.cinfo.density_unit > 0 &&
			(jhinf.cinfo.X_density > 0) &
			(jhinf.cinfo.Y_density > 0))
		pixAspect = float(jhinf.cinfo.Y_density) /
					float(jhinf.cinfo.X_density);
	else
		pixAspect = 1.f;
						// make color space consistent
	if (hrv == JPEG_HEADER_HDR) {
		PcopyCS(&cs, &ICS_RGB709);
		encoding = "HDR YCbCr";
	} else {
		if (!DetermineCS())		// determine Exif color space
			PcopyCS(&cs, &ICS_sRGB);
		switch (jhinf.cinfo.jpeg_color_space) {
		case JCS_YCbCr:
			encoding = "DCT YCbCr";
			break;
		case JCS_RGB:
			encoding = "DCT RGB";
			break;
		case JCS_GRAYSCALE:
			PcopyCS(&cs, &ICS_Y8);
			encoding = "DCT grayscale";
			break;
		case JCS_CMYK:
			cs.format = IPFcmyk;
			encoding = "DCT CMYK";
			break;
		default:			// in trouble with the rest
			sprintf(errMsg, "cannot convert JPEG color space (%d)",
					jhinf.cinfo.jpeg_color_space);
			errCode = IREunsupported;
			state = JPSerror;
			return false;
		}
	}
	fr.xleft = 0; fr.xright = xres;		// x range never changes
	fr.ytop = 0; fr.ybottom = 1;		// set y to top scanline
	nr = fr;				// first scanline is next
	return true;
}

// Start decompression cycle after header has been read
bool
JPEGReader::JPstart(bool getHDR)
{
	if (state != JPSstart)			// must be in start state!
		return false;
	jhinf.cinfo.dct_method = JDCT_FLOAT;
	if (getHDR && !HasHDR())
		return false;
						// our code needs scale_num=1
	jhinf.cinfo.scale_num = 1;
#if (JPEG_LIB_VERSION >= 90)
	jhinf.cinfo.scale_denom = 1;		// seems to be broken
#endif
	if (getHDR) {
		jpeghdr_start_decompress(&jhinf);
	} else {
		jpeg_start_decompress(&jhinf.cinfo);
#if SF_JPEG_IDCT
		jpeg_sf_idct(&jhinf.cinfo);	// set up custom decompression
#endif
	}
	nr.ytop = 0;				// reset scanline to top
	nr.ybottom = jhinf.cinfo.scale_denom;
	if (ldrscan != NULL) {			// (re)allocate scanline buffers
		free(ldrscan);
		ldrscan = NULL;
	}
	if (hdrscan != NULL) {
		free(hdrscan);
		hdrscan = NULL;
	}
	if (getHDR) {
		hdrscan = (JHSAMPLE *)malloc(jhinf.cinfo.output_width*
						3*sizeof(JHSAMPLE));
		if (hdrscan == NULL)
			goto memerr;
	} else {
		ldrscan = (JSAMPLE *)malloc(jhinf.cinfo.output_width*
			jhinf.cinfo.out_color_components*sizeof(JSAMPLE));
		if (ldrscan == NULL)
			goto memerr;
	}
	state = JPSscan;				// we're scanning, now
	return true;
memerr:
	strcpy(errMsg, "cannot allocate scanline buffer");
	errCode = IREmemory;
	return false;
}

// Read next scanline from JPEG input (into sl if not NULL)
bool
JPEGReader::JPnextScanline(JSAMPLE *sl, JHSAMPLE *hsl)
{
	if (state != JPSscan || nr.ytop >= yres)
		return false;
	if (ReadingHDR()) {
		if (hsl == NULL)
			hsl = hdrscan;
		if (!jpeghdr_read_scanline(&jhinf, hsl))
			return false;
		for (int x = (sl==NULL) ? 0 : xres; x-- > 0; )
			jpeghdr_scan_rgb24(&jhinf, sl+3*x, x);
	} else {
		if (sl == NULL)
			sl = ldrscan;
		if (!jpeg_read_scanlines(&jhinf.cinfo, &sl, 1))
			return false;
	}
	nr.ytop = nr.ybottom;
	nr.ybottom += jhinf.cinfo.scale_denom;
	return true;
}

/* Fix sampling region for image read buffer */
bool
JPEGReader::JPsetRec(ImgReadBuf *rb)
{
	if (rb->subsample <= 0)			// check/fix rectangle
		return false;
	if (!PlegalRect(&rb->r, xres, yres))
		return false;
	ImgFixSampling(rb);
	if (rb->r.ytop < nr.ytop)		// must abort current cycle
		if (!JPreset())
			return false;
	bool	getHDR = ((rb->cs.dtype == IDTfloat) && HasHDR());
	if (state == JPSstart) {		// new decompression cycle
		if (rb->subsample < 2)
			jhinf.cinfo.scale_denom = 1;
		else if (rb->subsample < 4)
			jhinf.cinfo.scale_denom = 2;
		else if (rb->subsample < 8)
			jhinf.cinfo.scale_denom = 4;
		else
			jhinf.cinfo.scale_denom = 8;
		switch (rb->cs.format) {
		case IPFrgb:
			if (!getHDR)
				jhinf.cinfo.out_color_space = JCS_RGB;
			break;
		case IPFy:
			getHDR = false;		// HDR Y not supported (yet)
			jhinf.cinfo.out_color_space = JCS_GRAYSCALE;
			break;
		case IPFcmy:
		case IPFcmyk:
			getHDR = false;		// HDR CMY really not supported
			if (jhinf.cinfo.jpeg_color_space == JCS_CMYK ||
					jhinf.cinfo.jpeg_color_space == JCS_YCCK)
				jhinf.cinfo.out_color_space = JCS_CMYK;
			break;
		default:			// forget the rest
			jhinf.cinfo.out_color_space = JCS_RGB;
			break;
		}
		if (!JPstart(getHDR))		// start 'er up
			return false;
	}
	while (nr.ytop < rb->r.ytop)		// advance to position
		if (!JPnextScanline())
			return false;
	return true;				// ready to read
}

// Read a rectangle from our JPEG
bool
JPEGReader::ReadRec(ImgReadBuf *rb)
{
	if ((state != JPSstart) & (state != JPSscan) || rb == NULL)
		return false;
	jumper_struct	myjumper;
	JumpSet(&myjumper);
	if (setjmp(myjumper.b)) {
		DiagnoseError();
		return false;
	}
	if (!JPsetRec(rb))
		return false;

	if (!PmatchColorSpace(&rb->cs, &cs, PICMformat|PICMchroma|PICMwhite) ||
			rb->cs.logorig > .0f || 
			(rb->cs.dtype == IDTfloat && !ReadingHDR())) {
		rb->cs = cs;			// make caller work
		if (!ReadingHDR())
			rb->cs.dtype = IDTubyte;
		rb->buf = NULL;
	}
	if (rb->subsample != (int)jhinf.cinfo.scale_denom) {
		rb->subsample = jhinf.cinfo.scale_denom;
		rb->buf = NULL;
	}
	if (!rb->buf) {				// need new buffer
		rb->r.xleft = 0;		// use scanline-sized buffer
		rb->r.xright = xres;
		rb->r.ybottom = rb->r.ytop + rb->subsample;
		int	buflen = ImgReadBufLen(rb);
		if (buflen <= 0) {
			strcpy(errMsg, "internal buffer size error");
			errCode = IREunknown;
			return false;
		}
						// do not free client memory!
		rb->buf = (uby8 *)malloc(buflen);
		if (rb->buf == NULL) {
			strcpy(errMsg, "cannot allocate new read buffer");
			errCode = IREmemory;
			return false;
		}
	}
						// read scanline(s)
	const int	nssamps = (rb->r.xright - rb->r.xleft)/rb->subsample *
					ImgPixelLen[rb->cs.format];
	const int	ssiz = nssamps * ImgDataSize[rb->cs.dtype];
	const int	nscn = (rb->r.ybottom - rb->r.ytop)/rb->subsample;
	int		nsrd = 0;
	JSAMPLE *	lscan = NULL;
	JHSAMPLE *	hscan = NULL;
	if ((rb->r.xleft == 0) & (rb->r.xright == xres)) {
		if (ReadingHDR()) hscan = (JHSAMPLE *)rb->buf;
		else lscan = (JSAMPLE *)rb->buf;
		for (nsrd = 0; nsrd < nscn; nsrd++) {
			if (!JPnextScanline(lscan, hscan))
				break;
			if (ReadingHDR()) hscan += nssamps;
			else lscan += nssamps;
		}
	} else {				// read partial scanlines
		void *	slp;
		if (ReadingHDR()) {
			hscan = hdrscan;
			slp = (void *)(hscan + rb->r.xleft/rb->subsample *
						ImgPixelSize(&rb->cs));
		} else {
			lscan = ldrscan;
			slp = (void *)(lscan + rb->r.xleft/rb->subsample *
						ImgPixelSize(&rb->cs));
		}
		for (nsrd = 0; nsrd < nscn; nsrd++) {
			if (!JPnextScanline(lscan, hscan))
				break;
			memcpy(rb->buf + nsrd*ssiz, slp, ssiz);
		}
	}
	JumpSet();
	return (nsrd == nscn);			// did we get them all?
}

// Diagnose a JPEG error and transfer to errMsg and errCode
void
JPEGReader::DiagnoseError()
{
	JumpSet();
	state = JPSerror;
						// copy jpeglib message
	strcpy(errMsg, jpeg_error_buffer);
						// examine code (jerror.h)
	switch (jerr.msg_code) {
					// unsupported features
	case JERR_SOF_UNSUPPORTED:
	case JWRN_ADOBE_XFORM:
		errCode = IREunsupported;
		break;
					// format errors
	case JERR_UNKNOWN_MARKER:
	case JERR_BAD_COMPONENT_ID:
	case JERR_BAD_DCT_COEF:
	case JERR_BAD_DCTSIZE:
	case JERR_BAD_HUFF_TABLE:
	case JERR_BAD_IN_COLORSPACE:
	case JERR_BAD_J_COLORSPACE:
	case JERR_BAD_LENGTH:
	case JERR_BAD_PRECISION:
	case JERR_COMPONENT_COUNT:
	case JERR_EMPTY_IMAGE:
	case JERR_EOI_EXPECTED:
	case JERR_SOF_DUPLICATE:
	case JERR_SOF_NO_SOS:
	case JERR_SOI_DUPLICATE:
	case JERR_NO_HUFF_TABLE:
	case JERR_NO_QUANT_TABLE:
	case JERR_NO_SOI:
	case JWRN_BOGUS_PROGRESSION:
	case JWRN_EXTRANEOUS_DATA:
	case JWRN_HIT_MARKER:
	case JWRN_HUFF_BAD_CODE:
	case JWRN_MUST_RESYNC:
	case JWRN_NOT_SEQUENTIAL:
	case JWRN_TOO_MUCH_DATA:
		errCode = IREformat;
		break;
					// Memory errors
	case JERR_OUT_OF_MEMORY:
		errCode = IREmemory;
		break;
					// Read errors
	case JERR_FILE_READ:
	case JERR_TFILE_CREATE:
	case JERR_TFILE_READ:
	case JERR_TFILE_SEEK:
	case JERR_TFILE_WRITE:
	case JERR_XMS_READ:
		errCode = IREread;
		break;
					// Truncation errors
	case JWRN_JPEG_EOF:
	case JERR_MISSING_DATA:
	case JERR_INPUT_EMPTY:
	case JERR_INPUT_EOF:
	case JERR_NO_IMAGE:
	case JERR_TOO_LITTLE_DATA:
		errCode = IREtruncated;
		break;
	default:			// Code error or ??
		errCode = IREunknown;
		break;
	}
}

// Check/load Exif color space
bool
JPEGReader::DetermineCS()
{
	if ((state != JPSstart) & (state != JPSscan))
		return false;
	jpeg_saved_marker_ptr	mp;
	for (mp = jhinf.cinfo.marker_list; mp != NULL; mp = mp->next)
		if (mp->marker == JPEG_APP0+1)		// Exif header
			break;
	if (mp == NULL || mp->data_length < 16)
		return false;
	if (memcmp(mp->data, "Exif\0\0", 6))
		return false;
	TIFFin	tin((const char *)mp->data+6, mp->data_length-6);
	return GetExifCS(&cs, &tin);
}

// Parse JPEG APP1 marker (Exif header)
static bool
getExif(const JOCTET FAR *data, unsigned int len, ImgInfo *info)
{
	if (len < 16)
		return false;		// too short even for header
	if (memcmp(data, "Exif\0\0", 6))
		return false;
	TIFFin	tin((const char *)data+6, len-6);
	return GetExifInfo(info, &tin);	// get info from Exif
}

// Parse JPEG APP0 marker (JFIF header)
static bool
getJFIF(const JOCTET FAR *data, unsigned int len, ImgInfo *info)
{
	if (len < 5)
		return false;			// too short even for header
	if (memcmp(data, "JFIF\0", 5))
		return false;
	// Don't know how to use this information right now...
	return true;
}

// Estimate JPEG quality setting from quantization table contents
static int
JRgetQuality(jpeg_decompress_struct *cinf)
{						// default table from IJG library
	static const UINT16 std_luminance_quant_tbl[DCTSIZE2] = {
		16,  11,  10,  16,  24,  40,  51,  61,
		12,  12,  14,  19,  26,  58,  60,  55,
		14,  13,  16,  24,  40,  57,  69,  56,
		14,  17,  22,  29,  51,  87,  80,  62,
		18,  22,  37,  56,  68, 109, 103,  77,
		24,  35,  55,  64,  81, 104, 113,  92,
		49,  64,  78,  87, 103, 121, 120, 101,
		72,  92,  95,  98, 112, 100, 103,  99
	};
	int	i;
	if (cinf->comp_info == NULL)
		return -1;			// no component info
	for (i = cinf->num_components; i-- > 0; )	// find Y channel
		if (cinf->comp_info[i].component_id == 1)
			break;
	if (i < 0)
		return -1;			// no Y channel
							// get Y quant. table
	int		tabno = cinf->comp_info[i].quant_tbl_no;
	JQUANT_TBL *	tbl;
	if ((tabno < 0) | (tabno >= NUM_QUANT_TBLS) ||
			(tbl = cinf->quant_tbl_ptrs[tabno]) == NULL)
		return -1;			// no quantization table
	float	scalesum = 0;				// average scalefactor
	int	nincluded = 0;
	int	row, col;
	for (row = 4; row--; )				// just LF components
		for (col = 4; col--; ) {
			i = 8*row + col;
			if (tbl->quantval[i] >= 255)
				continue;		// in case of baseline
			if (tbl->quantval[i] > 1)
				scalesum += (float)tbl->quantval[i] /
					    (float)std_luminance_quant_tbl[i];
			nincluded++;
		}
	if (!nincluded)
		return 0;			// QT max'ed out => 0 setting
	int	quality = (int)(100.f*scalesum/(float)nincluded + .5f);
	if (quality > 100)
		quality = 5000 / quality;			// setting < 50
	else
		quality = (200 - quality) / 2;			// setting >= 50
	return quality;				// return computed IJG setting
}

// Get information from Exif header (if one)
bool
JPEGReader::GetInfo(ImgInfo *info)
{
	if ((state != JPSstart) & (state != JPSscan))
		return false;
	if (info == NULL)
		return false;
	*info = defImgInfo;
	switch (jhinf.cinfo.density_unit) {
	case 1:						// pixels/inch
		info->hdensity = jhinf.cinfo.X_density;
		info->flags |= IIFhdensity;
		break;
	case 2:						// pixels/cm
		info->hdensity = jhinf.cinfo.X_density * 2.54f;
		info->flags |= IIFhdensity;
		break;
	}
	jpeg_saved_marker_ptr	mp;
	for (mp = jhinf.cinfo.marker_list; mp != NULL; mp = mp->next)
		switch (mp->marker) {
		case JPEG_APP0+0:			// JFIF header
			getJFIF(mp->data, mp->data_length, info);
			break;
		case JPEG_APP0+1:			// Exif header
			getExif(mp->data, mp->data_length, info);
			break;
		}
	info->quality = JRgetQuality(&jhinf.cinfo);
	info->flags |= IIFquality;
	if (HasHDR() && jhinf.samp2nits != JH_SAMP2NITS_UNK) {
		info->stonits = jhinf.samp2nits;
		info->flags |= IIFstonits;
	}
	return true;
}

// Get thumbnail from Exif if we can
bool
JPEGReader::GetThumbnail(ImgStruct *thm)
{
	if ((state != JPSstart) & (state != JPSscan))
		return false;
	jpeg_saved_marker_ptr	mp;
	for (mp = jhinf.cinfo.marker_list; mp != NULL; mp = mp->next)
		if (mp->marker == JPEG_APP0+1)
			break;
	if (mp == NULL || mp->data_length < 16)
		return false;				// no Exif marker
	TIFFin	tin((const char *)mp->data+6, mp->data_length-6);
	return GetExifThumbnail(thm, &tin);
}
