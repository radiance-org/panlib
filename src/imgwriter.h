/*
 *  imgwriter.h
 *  panlib
 *
 *  Depends on "imgio.h"
 *
 *  Declarations for image writing routines
 *
 *  Created by gward on Fri Jun 08 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */
#ifndef _IMGWRITER_H_
#define _IMGWRITER_H_

#ifdef __cplusplus
extern "C" {
#endif

/* ImgWriteBuf:
 *  Information and buffer pointers for writing out an image to a file.
 *  Pixels are interleaved in scanline order (i.e., left-to-right then
 *  top-to-bottom).  The PixAspect member gives the ratio of the vertical
 *  to horizontal density, which should be 1.0 for most formats.
 *  The rowsize member indicates how many bytes in each scanline,
 *  which must be greater than or equal to xres*ImgPixelSize(cs).
 *  The color definition is specified by csp, and additional
 *  information may be given in the info structure by setting
 *  the values and corresponding flags defined in "imgio.h".
 *  In particular, the info.quality member will be used to control
 *  output quality for compressed formats if set (0-100 range).
 */
typedef struct {
	int		xres, yres;	/* image size */
	int		rowsize;	/* row size (bytes) */
	float		pixAspect;	/* ratio of vertical/horiz. density */
	const ImgColorSpace *
			csp;		/* image color space */
	ImgInfo		info;		/* additional information */
	const uby8 *	img;		/* interleaved scanline data */
} ImgWriteBuf;

/* Prepare write buffer based on image structure */
extern int		PsetWriteBuf(ImgWriteBuf *iw, const ImgStruct *isrc);

/* ImgWriteInterface:
 *  The image writer interface is quite simple -- just one call
 *  to verify support of a given color space, and one more
 *  to write the image and return total size in bytes, or zero if
 *  there was an error.  The SupportedCS call returns a short text
 *  description of the format that will be output if the color space
 *  is supported, or NULL if it is not.  The second parameter to this
 *  function is the desired compression quality, or -1 for the default
 *  value.  The suffix member of the interface struct contains a text
 *  string with the recommended suffix for this format, which should
 *  be 3 characters in length.  If the file name passed to WriteImage
 *  has no suffix, beware that one will NOT be added.
 */
typedef struct {
	char		suffix[8];	/* recommended suffix for images */
	const char *	(*SupportedCS)(const ImgColorSpace *csp, int qual);
	long		(*WriteImage)(const char *fname, const ImgWriteBuf *wb);
} ImgWriterInterface;

extern const ImgWriterInterface	IWInterfaceJPEG;
extern const ImgWriterInterface	IWInterfaceTIFF;
extern const ImgWriterInterface	IWInterfaceRad;
extern const ImgWriterInterface IWInterfaceEXR;
extern const ImgWriterInterface	IWInterfaceDPT;
extern const ImgWriterInterface IWInterfaceNRM;
extern const ImgWriterInterface IWInterfaceMTX;

#ifdef __cplusplus
}
#endif

#endif	/* ! _IMGWRITER_H_ */
