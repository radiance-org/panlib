/*
 *  exifthumb.cpp
 *  panlib
 *
 *  Extraction routines for Exif headers
 *
 *  Created by gward on 6/9/2008.
 *  Copyright (c) 2008 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <string>
#include <sstream>
#include "jstreamsrc.h"
#include "tiffin.h"
#include "pimage.h"
#include "exif.h"

// Get Exif thumbnail from IFD1
bool
GetExifThumbnail(ImgStruct *thm, TIFFin *ti)
{
	if (thm == NULL || ti == NULL)
		return false;
	if (!ti->SetIFD(1))			// get IFD1 tags
		return false;
	int		tag;
	unsigned short	comp = 0;
	unsigned long	JPEGOffset = 0;
	unsigned long	JPEGSize = 0;
	unsigned long	xres = 0, yres = 0;
	for (tag = ti->FirstTag(); tag >= 0; tag = ti->NextTag())
		switch (tag) {
		case 0x0100:			// Image Width
			xres = ti->GetInteger();
			break;
		case 0x0101:			// Image Length
			yres = ti->GetInteger();
			break;
		case 0x0103:			// Compression
			ti->GetData(&comp, TIFFunsigned_short);
			break;
		case 0x0201:			// JPEG IFO Offset
			ti->GetData(&JPEGOffset, TIFFunsigned_long);
			break;
		case 0x0202:			// JPEG Byte Count
			ti->GetData(&JPEGSize, TIFFunsigned_long);
			break;
		}
	if ((comp != 6) | (JPEGOffset == 0) | (JPEGSize == 0))
		return false;			// only deal with JPEG
	if (sizeof(JSAMPLE) != 1)
		return false;			// code error!
	istream *	istr = ti->GrabStream(JPEGOffset);
	if (istr == NULL)
		return false;			// couldn't position pointer
	char *	jpbuf = new char [JPEGSize];
	istr->read(jpbuf, JPEGSize);
	ti->ReleaseStream(istr);
	if ((off_t)istr->gcount() != JPEGSize) {
		delete [] jpbuf;
		return false;
	}
						// read JPEG thumbnail
	istringstream			isstr(string(jpbuf,JPEGSize));
	struct jpeg_decompress_struct	cinfo;
	struct jpeg_error_mgr		jerr;
	jumper_struct			myjumper;
	myjumper.w = JMSG_NOMESSAGE;
	cinfo.err = jpeg_std_error(&jerr);
	jerr.output_message = jpeg_error_output;
	jerr.emit_message = jpeg_emit_message;
	jerr.error_exit = jpeg_error_jump;
	if (setjmp(myjumper.b))			// JPEG error jumps to here
		goto fail;
	cinfo.client_data = (void *)&myjumper;
	jpeg_create_decompress(&cinfo);
	jpeg_stream_src(&cinfo, &isstr);
	if (jpeg_read_header(&cinfo, TRUE) != JPEG_HEADER_OK)
		goto fail;
	cinfo.out_color_space = JCS_RGB;
	if ((thm->xres > 0) & (thm->yres > 0)) { // constrain output size
		if (cinfo.image_width*thm->yres > cinfo.image_height*thm->xres)
			while (cinfo.image_width > cinfo.scale_denom*thm->xres)
				cinfo.scale_denom <<= 1;
		else
			while (cinfo.image_height > cinfo.scale_denom*thm->yres)
				cinfo.scale_denom <<= 1;
		if (cinfo.scale_denom > 8)
			cinfo.scale_denom = 8;
	}
	jpeg_start_decompress(&cinfo);
	if (cinfo.out_color_components != 3)
		goto fail;
	if (!PmatchColorSpace(thm->csp, &ICS_sRGB, PICMptype))
		goto fail;
	if (thm->img != NULL && 
			(int)cinfo.output_width*(int)cinfo.output_height >
					thm->xres*thm->yres)
		PfreeImage(thm);
	thm->xres = (int)cinfo.output_width;
	thm->yres = (int)cinfo.output_height;
	if (!PnewImage(thm, .0))
		goto fail;
	{					// decompress thumbnail
		JSAMPLE *	ptr = (JSAMPLE *)thm->img;
		while (cinfo.output_scanline < cinfo.output_height) {
			jpeg_read_scanlines(&cinfo, &ptr, 1);
			ptr += thm->rowsize;
		}
	}
	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);
	delete [] jpbuf;
	return true;				// return success
fail:
	jpeg_destroy_decompress(&cinfo);	// something went wrong...
	delete [] jpbuf;
	PfreeImage(thm);
	return false;
}
