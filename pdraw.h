/*
 *  pdraw.h
 *  panlib
 *
 *  Drawing and stamping function declarations.
 *
 *  Created by Greg Ward on 1/24/19.
 *  Copyright 2019 Anyhere Software. All rights reserved.
 *
 */

#ifndef _PDRAW_H_
#define _PDRAW_H_

#include "panimage.h"
#include "abitmap.h"

/// Circle sections (quadrants) to draw
enum QuadSect { QSupperLeft=1, QSupperRight=2, QSlowerLeft=4, QSlowerRight=8,
		QSupper=3, QSlower=12, QSleft=5, QSright=10, QSfull=15 };

/// Draw (partial) circle into bitmap
extern int	BdrawCircle(ABitMap2 *map, int xc, int yc, int rad,
					int sc=QSfull, bool val=true);

/// Draw a line into bitmap
extern int	BdrawLine(ABitMap2 *map, int x0, int y0, int x1, int y1,
							bool val=true);

/// Draw a fat line into bitmap
extern int	BfillLine(ABitMap2 *map, int x0, int y0, int x1, int y1, int rad,
							bool val=true);

/// Draw a (rounded) rectangle into bitmap
extern int	BdrawRectangle(ABitMap2 *map, int x0, int y0, int x1, int y1,
						int rad=0, bool val=true);

/// Flood-fill area starting from seed point which must be inside closed region
extern int	BfloodFill(ABitMap2 *map, int xs, int ys, int val=-1);

/// Draw a polygon as a (thick) outline
extern int	BdrawPolygon(ABitMap2 *map, const int vxy[][2], int nv,
						int rad=0, bool val=true);

/// Fill a polygon (may have holes connected by seams)
extern int	BfillPolygon(ABitMap2 *map, const int vxy[][2], int nv, bool val=true);

/// Stamp the given color according to a bitmap with optional offset and downsample
extern int	PstampInk(PanImage *rimp, const ABitMap2 &stencil, PixelVal pv,
					int xleft=0, int ytop=0, const int overs = 1);

/// Convert image to bitmap(s) based on threshold color
extern bool	PthreshMap(ABitMap2 map[], const PanImage &im, PixelVal pv,
						PHistoChan chan=PHCrgb);

/// Stencil operator to overlay selected image pixels, no offset
extern bool	PstencilImage(PanImage *rimp, const ABitMap2 &stencil,
						const PanImage &sim);

/// Perform operation on selected image pixels; compare with PblendImageOpCB()
extern bool	PstencilImageOpCB(PanImage *idstp, const ABitMap2 mask,
						PimageOp *op, void *udp=0);

/// Resize a bitmap -- use a lower threshold to be more conservative with 1's
extern bool	BresizeBitmap(ABitMap2 *map, const ABitMap2 &src, float thresh=.5f);

/// Factor same-sized bitmap into image (i.e., paint "true" positions black)
inline PanImage
operator/=(PanImage &imLeft, const ABitMap2 &bmRight)
{
	if (imLeft.Ready() != PICready ||
			(imLeft.Width() != bmRight.Width()) |
			(imLeft.Height() != bmRight.Height()) )
		return imLeft;
	PstampInk(&imLeft, bmRight, Pblack);
	return imLeft;
}

inline PanImage
operator/(PanImage imLeft, const ABitMap2 &bmRight)
{
	return imLeft /= bmRight;
}

/// Multiply same-sized bitmap into image (i.e., paint "false" positions black)
inline PanImage
operator*=(PanImage &imLeft, ABitMap2 bmRight)
{
	bmRight.Invert();
	return imLeft /= bmRight;
}

inline PanImage
operator*(PanImage imLeft, const ABitMap2 &bmRight)
{
	return imLeft *= bmRight;
}

inline PanImage
operator*(const ABitMap2 &bmLeft, PanImage imRight)
{
	return imRight *= bmLeft;
}

#endif	// _PDRAW_H_
