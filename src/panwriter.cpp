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
bool
PloadStandardWriters()
{
	bool	ok = true;
	
	ok &= PanImage::AddImageWriter(&IWInterfaceJPEG);
	ok &= PanImage::AddImageWriter(&IWInterfaceTIFF);
	ok &= PanImage::AddImageWriter(&IWInterfaceRad);
	ok &= PanImage::AddImageWriter(&IWInterfaceEXR);
	// ok &= PanImage::AddImageWriter(&IWInterfaceDPT);
	// ok &= PanImage::AddImageWriter(&IWInterfaceNRM);
	// ok &= PanImage::AddImageWriter(&IWInterfaceMTX);
	
	return ok;
}
