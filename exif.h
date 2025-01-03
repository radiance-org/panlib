/*
 *  exif.h
 *  panlib
 *
 *  Include after "imgio.h" and "tiffin.h"
 *
 *  C++ function declarations for reading Exif header information
 *
 *  Created by Greg Ward on 8/19/06.
 *  Copyright 2006 Anyhere Software. All rights reserved.
 *
 */

// Read Exif color space
extern bool	GetExifCS(ImgColorSpace *ics, TIFFin *ti);

// Get Exif header information
extern bool	GetExifInfo(ImgInfo *info, TIFFin *ti);

// Read Exif thumbnail image
extern bool	GetExifThumbnail(ImgStruct *thm, TIFFin *ti);
