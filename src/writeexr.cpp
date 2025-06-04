/*
 *  writeexr.cpp
 *  panlib
 *
 *  Routines for writing out EXR image files.
 *
 *  Created by gward on Aug 26 2002.
 *  Copyright (c) 2002 Anyhere Software. All rights reserved.
 *
 */

#ifdef freebsd
#define PLATFORM_DARWIN_PPC
#endif

#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include "system.h"
#include "imgio.h"
#include "imgwriter.h"
#include "color.h"
#include <ImfRgbaFile.h>
#include <ImfStandardAttributes.h>

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

using namespace	Imf;
using namespace std;

static COLORMAT	xyz2srgbWB;		// standard conversion matrix

// Fill header attributes and return channels to write or -1 on error
static int
ExrFinishHeader(Header *hdr, const ImgWriteBuf *wb)
{
	if (wb->csp->dtype != IDTfloat)
		return -1;			// don't handle anything else
	if (wb->info.flags & IIFcrop) {
		hdr->displayWindow().min.x = wb->info.crop.xleft;
		hdr->displayWindow().max.x = wb->info.crop.xright - 1;
		hdr->displayWindow().min.y = wb->info.crop.ytop;
		hdr->displayWindow().max.y = wb->info.crop.ybottom - 1;
	}
	if (wb->info.flags & IIFcapdate)
		addCapDate(*hdr, string(wb->info.capdate));
		// hdr->insert("capDate", string(wb->info.capdate));
	if (wb->info.flags & IIFlatlong) {
		addLatitude(*hdr, wb->info.latlong[0]);
		// hdr->insert("latitude", wb->info.latlong[0]);
		addLongitude(*hdr, wb->info.latlong[1]);
		// hdr->insert("longitude", wb->info.latlong[1]);
	}
	if (wb->info.flags & IIFstonits)
		addWhiteLuminance(*hdr, wb->info.stonits);
		// hdr->insert("whiteLuminance", wb->info.stonits);
	if (wb->info.flags & IIFhdensity)
		addXDensity(*hdr, wb->info.hdensity);
		// hdr->insert("xDensity", wb->info.hdensity);
	if (wb->info.flags & IIFexptime)
		addExpTime(*hdr, wb->info.exptime);
		// hdr->insert("expTime", wb->info.exptime);
	if (wb->info.flags & IIFaperture)
		addAperture(*hdr, wb->info.aperture);
		// hdr->insert("aperture", wb->info.aperture);
	if (wb->info.flags & IIFasa)
		addIsoSpeed(*hdr, wb->info.asa);
		// hdr->insert("isoSpeed", wb->info.asa);
	if (wb->info.flags & IIFfocus)
		addFocus(*hdr, wb->info.focus);
		// hdr->insert("focus", wb->info.focus);
	if (wb->info.flags & IIFview)
		hdr->insert("view", StringAttribute(wb->info.view));
	if (wb->info.flags & IIFowner)
		addOwner(*hdr, string(wb->info.owner));
		// hdr->insert("owner", string(wb->info.owner));
	if (wb->info.flags & IIFcomments)
		addComments(*hdr, string(wb->info.comments));
		// hdr->insert("comments", string(wb->info.comments));
	if (wb->csp->format == IPFrgb) {
		Chromaticities  chrom;
		chrom.red[0] = wb->csp->chroma[0][0];
		chrom.red[1] = wb->csp->chroma[0][1];
		chrom.green[0] = wb->csp->chroma[1][0];
		chrom.green[1] = wb->csp->chroma[1][1];
		chrom.blue[0] = wb->csp->chroma[2][0];
		chrom.blue[1] = wb->csp->chroma[2][1];
		chrom.white[0] = wb->csp->chroma[3][0];
		chrom.white[1] = wb->csp->chroma[3][1];
		addChromaticities(*hdr, chrom);
		// hdr->insert("chromaticities", chrom);
		return WRITE_RGB;
	}
	if (wb->csp->format == IPFxyz) {
		compxyz2rgbWBmat(xyz2srgbWB, (const RGBPRIMP)ICS_RGB709.chroma);
		return WRITE_RGB;
	}
	if (wb->csp->format == IPFy)
		return WRITE_Y;
//		return WRITE_G;	// should be WRITE_Y but fails for some reason
	return -1;
}

// Write out the next scanline
static void
ExrWriteScan(RgbaOutputFile *fo, const ImgWriteBuf *wb, const int y)
{
	Rgba		shortScan[4096];
	Rgba *		rgbaScan = shortScan;
	if (wb->xres > 4096)
		rgbaScan = new Rgba [wb->xres];
	const float *	pp = (const float *)(wb->img + (ssize_t)y*wb->rowsize);
	int		x;
	switch (wb->csp->format) {
	case IPFrgb:
		for (x = 0; x < wb->xres; x++) {
			rgbaScan[x].r = *pp++;
			rgbaScan[x].g = *pp++;
			rgbaScan[x].b = *pp++;
			rgbaScan[x].a = 1.f;
		}
		break;
	case IPFxyz:
		for (x = 0; x < wb->xres; x++) {
			COLOR	clr;
			colortrans(clr, xyz2srgbWB, const_cast<float *>(pp));
			rgbaScan[x].r = clr[RED];
			rgbaScan[x].g = clr[GRN];
			rgbaScan[x].b = clr[BLU];
			rgbaScan[x].a = 1.f;
			pp += 3;
		}
		break;
	case IPFy:
		for (x = 0; x < wb->xres; x++) {
			rgbaScan[x].r =
			rgbaScan[x].g = 
			rgbaScan[x].b = *pp++;
			rgbaScan[x].a = 1.f;
		}
		break;
	default:
		break;
	}
	fo->setFrameBuffer(rgbaScan, 1, 0);
	fo->writePixels();
	if (rgbaScan != shortScan)
		delete [] rgbaScan;
}

// Check if the given color space is supported by our writer
static const char *
ExrSupportedCS(const ImgColorSpace *csp, int)
{
	if (csp->logorig > 0)
		return NULL;
	if (csp->dtype != IDTfloat)
		return NULL;			// float is all we handle
	if ((csp->format == IPFxyz) | (csp->format == IPFrgb))
		return "EXR PIZ RGB";
	if (csp->format == IPFy)
		return "EXR PIZ Luminance";
	return NULL;
}	

// Write out a EXR image
static long
ExrWriteImage(const char *fname, const ImgWriteBuf *wb)
{
					// check arguments
	if ((fname == NULL) | (wb == NULL) || !*fname)
		return 0;
					// create EXR header
	float			screenWidth = 1;
	if (wb->info.flags & IIFhvangle)
		screenWidth = 2.*tan(M_PI/180./2.*wb->info.hvangle);
	Header			myHead(wb->xres, wb->yres, wb->pixAspect,
					Imath::V2f(0,0), screenWidth,
					INCREASING_Y,
					PIZ_COMPRESSION);
	int			outChan = ExrFinishHeader(&myHead, wb);
	RgbaOutputFile *	fileOut = NULL;
	if (outChan < 0)
		return 0;
	try {
		fileOut = new RgbaOutputFile(fname, myHead, (RgbaChannels)outChan);
		for (int y = 0; y < wb->yres; y++)
			ExrWriteScan(fileOut, wb, y);
	}
	catch (const exception &e) {
#ifndef NDEBUG
		cerr << e.what() << '\n';
#endif
		delete fileOut;
		return 0;
	}
	delete fileOut;			// does this generate exceptions?
	int		fd = open(fname, O_RDONLY);
	if (fd < 0)
		return 0;
	SET_FD_BINARY(fd);
	off_t		flen = lseek(fd, 0L, SEEK_END);
	close(fd);
	return (long)flen;
}

extern "C" {
extern const ImgWriterInterface	IWInterfaceEXR;
const ImgWriterInterface	IWInterfaceEXR = {
	"exr", &ExrSupportedCS, &ExrWriteImage
};
}
