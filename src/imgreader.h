/*
 *  imgreader.h
 *  panlib
 *
 *  Standard interface for image file readers in Pancine.
 *
 *  Created by gward on Tue May 01 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#ifndef _IMGREADER_H_
#define _IMGREADER_H_

#include "imgio.h"

#ifdef __cplusplus
#define EXTRN(type)	extern "C" type
#else
#define EXTRN(type)	extern type
#endif

/* Error Codes:  none, unsupported, unexpected end of file, format error,
 *		read/seek error, out of memory, unknown error
 */
typedef enum { IREnone, IREunsupported, IREtruncated, IREformat,
			IREread, IREmemory, IREunknown }	ImgReadErr;

/* Frame Types:  lone image, layered image, animated sequence, circular loop,
 *		loop back-and-forth, 3-D slice, unknown type
 */
typedef enum { IRFnone, IRFlayer, IRFanim, IRFloop,
			IRFbackforth, IRFslice, IRFunknown }	ImgFrameType;

/* Seek Modes:	absolute, relative, advance
 */
typedef enum { IRSabs, IRSrel, IRSadv }	ImgSeekMode;

struct _ImgReaderInterface;		/* forward declaration */

/* ImgReader:
 *  Allocated structure for an open image file.
 *  The file member contains the file name or specification.
 *  The current frame number (starts at 0) and total number
 *  of frames is given in frame and nframes, respectively.
 *  If the number of frames is unknown, nframes may be zero.
 *  The sequence type is specified in frameType, and
 *  the recommended frame rate (frames per second)
 *  is given in frameRate.  This can be negative for
 *  reverse sequences or zero to disable animation.
 *  The horizontal image dimension is given in xres,
 *  and the vertical resolution is returned in yres.
 *  Image orientation is always such that the (0,0)
 *  coordinate is in the upper-left corner, with the
 *  x-coordinate proceeding to the right, and the
 *  y-coordinate proceeding down in the image.
 *  (A fallback is provided for readers that cannot reorient
 *  in the "orientation" member of the ImgInfo struct.)
 *  The pixel aspect ratio (vert_density/horiz_density)
 *  is given in the pixAspect member, usually 1.0.
 *  The most natural color space and data type for
 *  this image are stored in the cs structure.
 *  A more specific description of the file encoding
 *  may be indicated with the "encoding" string.
 *  The fr rectangle specifies the position and
 *  dimensions of the first, smallest read buffer
 *  in the natural ordering of the input image.
 *  All readers must support reading this first
 *  scanline after any sequence of operations.
 *  The nr rectangle specifies the position and
 *  dimensions of the smallest next read buffer
 *  based on the current read pointer position,
 *  and will be outside the image at end-of-file.
 *  The errCode member is assigned by each call
 *  to one of the codes given in the ImgReadErr enum above
 *  (IREnone if the call was successful).  Additional
 *  error information may be stored in the nul-terminated
 *  string, errMsg, with no final newline.  The base
 *  structure below is designed to be extended to
 *  include additional (private) members needed by
 *  the reader, such as the open file pointer.
 */
#define BASE_IMGREADER	const struct _ImgReaderInterface *ri; \
			char		file[1024]; \
			int		frame, nframes; \
			ImgFrameType	frameType; \
			float		frameRate; \
			int		xres, yres; \
			float		pixAspect; \
			ImgColorSpace	cs; \
			const char *	encoding; \
			ImgRect		fr; \
			ImgRect		nr; \
			ImgReadErr	errCode; \
			char		errMsg[256]

typedef struct {
	BASE_IMGREADER;			/* first in derived struct's */
} ImgReader;

/* ImgReadBuf:
 *  Passed structure for reading rectangle from an image layer.
 *  The r member specifies the rectangle to be read.
 *  Scanline data must be packed (e.g., RGBRGBRGB),
 *  with each row proceeding in left-to-right order.
 *  If subsample is greater than 1, then the caller is
 *  requesting that the reader subsample the image using
 *  the given denominator for both the horizontal and
 *  vertical dimensions.  The boundary rectangle coordinates
 *  MUST BE exact multiples of subsample or bad things happen!
 *  The calling application owns the structure and may set a
 *  buffer with the rectangle and format it most desires,
 *  which is usually taken from the nr and cs members of ImgReader.
 *  If the buffer is initially NULL, it will be allocated by the
 *  reader using a call to malloc(), and should be free'd by caller.
 *  The image reader is free to change any unexpected
 *  parameters to ones it can deal with, and the calling
 *  application must make do.  If the caller assigned a buffer and
 *  any of the original parameters are changed, a new buffer will
 *  be allocated with malloc.  (Reader does NOT deallocate the
 *  passed buffer!  If it allocates and assigns a new buffer,
 *  the caller should free both the new one and the old one.)
 *  The basic requirement is that the reader must return a buffer
 *  having a useful degree of overlap with the requested rectangle
 *  (i.e., all of it plus some), or set the ImgReader errCode
 *  and errMsg members indicating what went wrong.  Readers
 *  _must_ support reading the next valid rectangle in the nr
 *  member of the base structure, as well as reading the first
 *  rectangle in the fr member at any point.  They should also
 *  support reading the entire image in the native color space.
 *  A good citizen also supports conversion to sRGB and Y8 for
 *  native color and grayscale spaces, respectively.
 *  Use the IRmoreRec() macro to detect the end of file.
 */
typedef struct {
	ImgColorSpace	cs;		/* color encoding */
	ImgRect		r;		/* image rectangle */
	int		subsample;	/* subsampling rate */
	uby8 *		buf;		/* read buffer */
} ImgReadBuf;

#define ImgReadBufLen(rb)	( (size_t)ImgPixelSize(&(rb)->cs) * \
				((rb)->r.xright/(rb)->subsample - \
					(rb)->r.xleft/(rb)->subsample) * \
				((rb)->r.ybottom/(rb)->subsample - \
					(rb)->r.ytop/(rb)->subsample) )

/* Call to make sure sampling is divisible */
EXTRN(void)	ImgFixSampling(ImgReadBuf *rb);

/* ImgReaderInterface:
 *  The reader interface function table.  The only required calls are
 *  Open(), ReadRec(), and Close().  Optional calls provide additional
 *  information about the image and allow multiple images and thumbnails
 *  to be loaded from layers or linked files.
 */
typedef struct _ImgReaderInterface {
	const char *	suffixList;
	/*  The suffixList pointer is assigned to a statically allocated,
	 *  nul-terminated string with a list of suffixes commonly
	 *  associated with the file format(s) this reader can handle.
	 *  Words are separated by a period ('.'), and the first word
	 *  is used as the main identifier for this image type as well
	 *  as a suffix.  Each alternate suffix begins with a period.
	 *  Suffixes may be any length, and are case-insensitive.
	 *  For example, a TIFF reader might use "TIFF.tif" as its
	 *  suffix list.
	 */
	ImgReader *	(*Open)(const char *fname);
	/*  The Open call allocates an ImgReader structure and
	 *  attempts to open the named image file, assigning the color
	 *  space parameters to those most natural for the first image
	 *  layer.  If the file cannot be opened or is clearly the wrong type,
	 *  a NULL pointer may be returned.  If there is some more subtle
	 *  problem, an error code should be set in the returned structure with
	 *  an explanation stored in errMsg to aid the user.
	 *  If all is well, errCode should be set to IREnone.  See the
	 *  description of the ImgReader structure above for more details.
	 */
	ImgReadErr	(*SeekFrame)(ImgReader *ir, int offs, ImgSeekMode sm);
	/*  The optional SeekFrame call advances to the given offset
	 *  in a sequence of image frames or data layers.  This call
	 *  is needed to access alpha channels and multi-spectral
	 *  image channels as well as associated depth layers.
	 *  The ImgReader struct is passed and modified with
	 *  a possibly changed file specification, size, and color
	 *  space.  A zero (IREnone) value should be returned unless
	 *  there is some file error, when an ImgReadErr code should
	 *  be returned and set in the ImgReader struct.  If the seek
	 *  is outside the range of images or layers  available, the
	 *  call should return IREtruncated.  (If the seek mode is
	 *  IRSadv and the frameType is IRFloop or IRFbackforth, then
	 *  the reader should continue advancing in a loop.)
	 */
	ImgReadErr	(*GetInfo)(ImgReader *ir, ImgInfo *info);
	/*  The optional GetInfo call returns metadata about the current
	 *  image layer, such as recorded date and time, exposure mode,
	 *  f-stop, and shutter speed.  Information the reader can
	 *  extract is set in the ImgInfo structure along with flags
	 *  indicating which parameters were set.  (See the description
	 *  of the ImgInfo struct in imgio.h.)  The ImgInfo structure
	 *  is assumed to be empty before the call.  A value of IREnone will
	 *  be returned normally, even if no metadata was extracted.
	 *  If there is a serious error, an ImgReadErr code should
	 *  be returned and set in the ImgReader struct.
	 */
	ImgReadErr	(*GetThumbnail)(ImgReader *ir, ImgStruct *tn);
	/*  The optional GetThumbnail call provides a means to quickly
	 *  read an image's thumbnail into the passed sRGB image.
	 *  Zero (IREnone) should be returned on success, and one
	 *  of the other ImgReadErr codes should be returned on failure,
	 *  setting this value also in the passed ImgReader structure.
	 *  If there is no thumbnail, then IREunsupported should be
	 *  returned.
	 */
	void		(*ToneMapping)(ImgReader *ir, double LdDyn,
					double LdMax, int humanVis);
	/*  The optional ToneMapping call is used to advise about tone-
	 *  mapping requirements for a specific display application.
	 *  On readers that perform floating-point to sRGB conversion
	 *  for display, this call helps decide how to map pixel
	 *  values appropriately.  The LdDyn parameter gives the
	 *  ratio between the maximum and minimum displayable luminance.
	 *  The LdMax parameter gives the maximum displayable luminance
	 *  in candelas/sq.meter (nits).  The humanVis boolean is non-zero
	 *  if the caller wishes to match world to display visibility.
	 */
	ImgReadErr	(*ReadRec)(ImgReader *ir, ImgReadBuf *rb);
	/*  The ReadRec call is the basic mechanism for reading a
	 *  block of pixel data from the current image layer.
	 *  The call should return zero (IREnone) on success, an
	 *  ImgReadErr value on failure, and should always set errCode.
	 *  See the lengthy description of the ImgReadBuf structure
	 *  above for details.
	 */
	void		(*Close)(ImgReader *ir);
	/*  The Close call closes the image if open, cleans up and
	 *  frees associated data and the ImgReader structure itself.
	 *  No error return value is permitted, and this call must
	 *  clean up and free its resources regardless of any 
	 *  problems encountered in the process.
	 */
} ImgReaderInterface;

#define IRseekFrame(ir,offs,sm)	((ir)->ri->SeekFrame==NULL ? IREunsupported \
					: (*(ir)->ri->SeekFrame)(ir,offs,sm))
#define IRgetInfo(ir,info)	((ir)->ri->GetInfo==NULL ? IREunsupported \
					: (*(ir)->ri->GetInfo)(ir,info))
#define IRgetThumbnail(ir,tn)	((ir)->ri->GetThumbnail==NULL ? IREunsupported \
					: (*(ir)->ri->GetThumbnail)(ir,tn))
#define IRtoneMapping(ir,d,m,h)	if ((ir)->ri->ToneMapping!=NULL) \
					(ir)->ri->ToneMapping(ir,d,m,h); else
#define IRreadRec(ir,rb)	(*(ir)->ri->ReadRec)(ir,rb)
#define IRmoreRec(ir)		PlegalRect(&(ir)->nr,(ir)->xres,(ir)->yres)
#define IRclose(ir)		if ((ir)!=NULL) (*(ir)->ri->Close)(ir); else

#endif /* ! _IMGREADER_H_ */
