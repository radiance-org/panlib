/*
 *  panwriter.cpp
 *  panlib
 *
 *  Load standard image writers for PanImage class.
 *
 *  Created by Greg Ward on 6/24/06.
 *  Copyright 2006 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include "panimage.h"

// Load standard image writers
bool PloadStandardWriters()
{
    bool ok = true;
    
    #ifdef HAVE_JPEG_SUPPORT
    ok &= PanImage::AddImageWriter(&IWInterfaceJPEG);
    #endif

    #ifdef HAVE_TIFF_SUPPORT
    ok &= PanImage::AddImageWriter(&IWInterfaceTIFF);
    #endif

    ok &= PanImage::AddImageWriter(&IWInterfaceRad);

    #ifdef HAVE_OPENEXR_SUPPORT
    ok &= PanImage::AddImageWriter(&IWInterfaceEXR);
    #endif
    
    return ok;
}
