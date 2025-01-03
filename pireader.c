/*
 *  pireader.c
 *  panlib
 *
 *  Separate call to load known image readers so pulling them in is optional.
 *
 *  Created by Gregory Ward on Feb 15, 2006.
 *  Copyright (c) 2006 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include "pimage.h"

/* Add canonical reader interfaces to global table */
int
PloadStandardReaders()
{
	int	ok = 1;

	ok &= (PaddIReaderI(&IRInterfaceTIFF) >= 0);
	ok &= (PaddIReaderI(&IRInterfaceJPEG) >= 0);
	ok &= (PaddIReaderI(&IRInterfaceRad) >= 0);
	ok &= (PaddIReaderI(&IRInterfaceEXR) >= 0);
	ok &= (PaddIReaderI(&IRInterfaceBMP) >= 0);
	/* ok &= (PaddIReaderI(&IRInterfaceDPT) >= 0); */
	/* ok &= (PaddIReaderI(&IRInterfaceNRM) >= 0); */
	/* ok &= (PaddIReaderI(&IRInterfaceMTX) >= 0); */

	return ok;
}
