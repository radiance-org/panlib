/*
 *  imgio.h
 *  panlib
 *
 *  Common header for image i/o in Pancine.
 *
 *  Created by gward on Sat May 05 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#ifndef _IMGIO_H_
#define _IMGIO_H_

#ifndef _MEMOBJECT_H_
#include "memobject.h"
#endif

#ifndef DEF_IQUALITY
#define DEF_IQUALITY		90		/* default quality setting */
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 *  Pixel data must be transferred using a type and
 *  format chosen from the following options.
 */

#undef uby8
#define uby8	unsigned char			/* basic image data type */

/* Data Types:  unsigned 8-bit, unsigned 16-bit, 32-bit integer, IEEE float */
typedef enum { IDTubyte, IDTushort, IDTint, IDTfloat } ImgDataType;

/* Pixel Formats:  RGB, CMY, CMYK, CIE XYZ, grayscale, RGBalpha,
	alpha, distance, parametric coords, 3-D vector, encoded normal */
typedef enum { IPFrgb, IPFcmy, IPFcmyk, IPFxyz, IPFy, IPFrgba,
		IPFa, IPFd, IPFuv, IPFvec3, IPFn } ImgPixelFormat;

/* ImgColorSpace:
 *  This structure completely defines the
 *  color space being used for communication.
 *  The gamma value is the exponent to which
 *  each component is raised to get to a linear
 *  response (approximate).  If logorig is positive,
 *  gamma is interpreted as the dynamic range of a
 *  logarithmic encoding in orders of magnitude, where
 *  each sample p evaluates to logorig*10^(gamma*p).
 *  The chroma pairs are the x-y chromaticity coordinates
 *  in 1931 CIE xy color space for each of the (up to) four
 *  primaries.  If there are three primaries, the fourth
 *  specification is white.  The white chromaticity is the
 *  reference color corresponding to equal primary values
 *  (e.g., RGB=[255,255,255]).
 */
typedef struct {
	ImgDataType	dtype;		/* data type */
	ImgPixelFormat	format;		/* pixel format */
	short		cstatic;	/* is struct static? */
	float		gamma;		/* gamma value or dynamic range */
	float		logorig;	/* origin of log encoding if > 0 */
	float		chroma[4][2];	/* primary xy pairs */
} ImgColorSpace;

extern const int	ImgDataSize[];	/* data size (bytes) */
extern const int	ImgPixelLen[];	/* component count */

#define MaxDataSize	sizeof(float)	/* maximum data size */
#define MaxPixelLen	4		/* maximum # components */
#define MaxPixelSize	(MaxDataSize*MaxPixelLen)

#define ImgPixelSize(c)	(ImgDataSize[(c)->dtype]*ImgPixelLen[(c)->format])

		/* sRGB:  24-bit RGB, 709 primaries, D65 white, 2.2 gamma */
extern const ImgColorSpace	ICS_sRGB;
		/* AdoobeRGB:	24-bit RGB, Adobe primaries, D65 white, 2.2 gamma */
extern const ImgColorSpace	ICS_AdobeRGB;
		/* RGB98Adobe:  96-bit linear RGB, Adobe 1998, D65 white */
extern const ImgColorSpace	ICS_RGB98Adobe;
		/* ProPhotoRGB:	96-bit ProPhoto RGB (Reference Output Medium Metric) */
extern const ImgColorSpace	ICS_ProPhotoRGB;
		/* P3:	floating point with P3 primaries & D50 white point */
extern const ImgColorSpace	ICS_P3;
		/* RGB2020:	floating point with BT.2020 primaries */
extern const ImgColorSpace	ICS_RGB2020;
		/* RGB709:  96-bit linear RGB, 709 primaries, D65 white */
extern const ImgColorSpace	ICS_RGB709;
		/* RGB709dB:  logarithmic RGB floating point (dB) */
extern const ImgColorSpace	ICS_RGB709dB;
		/* RGBradiance:	96-bit linear RGB, Radiance, neutral white */
extern const ImgColorSpace	ICS_RGBradiance;
		/* SharpRGB:  Sharpened RGB color space for white balancing */
extern const ImgColorSpace	ICS_SharpRGB;
		/* XYZ:  96-bit XYZ (CIE 1931 2-degree) */
extern const ImgColorSpace	ICS_XYZ;
		/* XYZdB:  96-bit floating point CIE XYZ color space (dB) */
extern const ImgColorSpace	ICS_XYZdB;
		/* Y8:  8-bit luminance, 2.2 gamma */
extern const ImgColorSpace	ICS_Y8;
		/* Y16:  16-bit linear luminance */
extern const ImgColorSpace	ICS_Y16;
		/* LogL16:  16-bit logarithmic grayscale matching LogLuv */
extern const ImgColorSpace	ICS_LogL16;
		/* Y:  32-bit floating-point linear luminance */
extern const ImgColorSpace	ICS_Y;
		/* Alpha:  32-bit floating-point linear alpha */
extern const ImgColorSpace	ICS_Alpha;
		/* Log10Y:  32-bit floating-point log10(Y)/10 */
extern const ImgColorSpace	ICS_YdB;
		/* D:  32-bit floating-point depth map */
extern const ImgColorSpace	ICS_D;
		/* VEC3:  3-D floating-point vector map */
extern const ImgColorSpace	ICS_VEC3;
		/* ENORM:  32-bit encoded surface normal */
extern const ImgColorSpace	ICS_ENORM;

/* Copy color space, unsetting cstatic field */
#define PcopyCS(d,s)	(*(d) = *(s), (d)->cstatic=0)

/* Change color space to linear float */
#define PrealCS(d)	((d)->dtype=IDTfloat, (d)->gamma=1, (d)->logorig=0)

/* Flags to determine what to match between two color spaces
 */
typedef enum {
	PICMok=0x0,		/* just check pointers != NULL */
	PICMdtype=0x1,
	PICMformat=0x2,
	PICMgamma=0x4,
	PICMchroma=0x8,
	PICMwhite=0x10,
	PICMptype=0x3,		/* match pixel type and format */
	PICMprims=0x18,		/* match chroma and white point */
	PICMall=0x1f,		/* match everything */
	PICMequiv=0x20		/* color conversion would be no-op */
} PICMatch;

/* Return string with description of color space (str can be NULL) */
char *			PdescribeCS(const ImgColorSpace *csp, char *str);

/* Check two color spaces to see if they match */
int			PmatchColorSpace(const ImgColorSpace *cs1,
					const ImgColorSpace *cs2, int cm);

/* ImgRect:
 *  A subimage rectangle for i/o buffering operations.
 *  Subimage ordering is always left-to-right scanlines
 *  proceeding top to bottom.  X coordinate goes from xleft
 *  at left of subimage to xright-1 at right.  Y coordinate
 *  goes from ytop at top of image to ybottom-1 at the bottom.
 */
typedef struct {
	int		xleft, xright, ytop, ybottom;
} ImgRect;

/* Check that the given rectangle is non-empty and within (0,0) to (xm,ym) */
#define PlegalRect(r,xm,ym)	((r) != NULL && ((r)->xleft < (r)->xright) & \
				((r)->ytop < (r)->ybottom) && \
				((r)->xleft >= 0) & ((r)->xright <= (xm)) && \
				((r)->ytop >= 0) & ((r)->ybottom <= (ym)))

/* Check that a point is within the given rectangle */
#define PinsideRect(r,x,y)	(((r)->xleft <= (x)) & ((x) < (r)->xright) && \
				((r)->ytop <= (y)) & ((y) < (r)->ybottom))

/* Number of pixels in rectangle (assumes legal) */
#define PrectArea(r)		( (long)((r)->xright - (r)->xleft) * \
				  ((r)->ybottom - (r)->ytop) )

/* Check if two rectangles are the same */
#define PequalRect(r1,r2)	!memcmp(r1, r2, sizeof(ImgRect))

/* Assign rectangle based on inclusive corners */
ImgRect		*PsetRect(ImgRect *dr, int x1, int y1, int x2, int y2);

/* Clip one rectangle to boundaries of another, returning 0 if nothing left */
int		PclipRect(ImgRect *dr, const ImgRect *cr);

/* ImgStruct:
 *  Passed structure for returning an image or thumbnail.
 *  The desired image dimensions are passed in xres and yres,
 *  and these will be altered to the returned image dimensions.
 *  The img pointer points to the upper left pixel,
 *  with scanlines proceeding left to right, top to bottom.
 *  The pixels are assumed to be square (pixAspect==1).
 *  The mbase member points to the underlying memory object if one.
 *  The rowsize member indicates how many bytes to skip to
 *  get to the next row.  This may be greater than xres*PixelSize
 *  for subimage areas or negative for reverse-ordered images.
 */
typedef struct {
	int		xres, yres;	/* image resolution */
	const ImgColorSpace *
			csp;		/* color space reference */
	int		rowsize;	/* row size (bytes) */
	MemObject *	mbase;		/* base memory object */
	uby8 *		img;		/* image pixel data */
} ImgStruct;

/* PixelVal:
 *  A single pixel value, with a pointer to the defining color space.
 */
typedef struct {
	const ImgColorSpace *
			csp;		/* color space reference */
	union {
		uby8		b[MaxPixelLen];
		unsigned short	s[MaxPixelLen];
		float		f[MaxPixelLen];
	} v;				/* encoded pixel value */
} PixelVal;

/** Image pixel pointers */
#define ProwPtr(ib,y)	((ib)->img + (long)(y)*(ib)->rowsize)
#define PpixP(ib,x,y,s)	(ProwPtr(ib,y) + (x)*(s))
#define PpixPtr(ib,x,y)	PpixP(ib,x,y,ImgPixelSize((ib)->csp))

/** Get position from pointer */
extern int		PpixPos(int *xp, int *yp,
				const ImgStruct *ia, const uby8 *ptr);

/** Extract a single pixel from the given image */
extern PixelVal		PgetPixel(const ImgStruct *ia, int x, int y);

/** Assign a single pixel in the given image */
extern int		PsetPixel(ImgStruct *ia, int x, int y, PixelVal pv);

/** Check pixel values for equality */
extern int		PequalVal(PixelVal p1, PixelVal p2);

extern const PixelVal	Pblack;		/* zero in any color space */

/** Convenient routine to set 24-bit sRGB pixel */
extern PixelVal		PsRGB(int r, int g, int b);

/** Convenient routine to set floating-point XYZ pixel */
extern PixelVal		PXYZ(float X, float Y, float Z);

/** Convenient routine to set floating-point Y pixel */
extern PixelVal		PgraY(float Y);

/* Image Information Flags:  image source, calibration factor, capture date,
 *		horizontal density, horizontal view angle, focus distance,
 *		GMT date, latitude+longitude, altitude,
 *		flash mode, exposure time, aperture,
 *		ASA, white balance, view parameters,
 *		owner, parameters, comments,
 *		image quality, crop rectangle,
 *		image orientation, sight bearing, focal length
 */
typedef enum	{	IIFsource=0x1, IIFstonits=0x2, IIFcapdate=0x4,
			IIFhdensity=0x8, IIFhvangle=0x10, IIFfocus=0x20,
			IIFgmtdate=0x40, IIFlatlong=0x80, IIFaltitude=0x100,
			IIFflash=0x200, IIFexptime=0x400, IIFaperture=0x800,
			IIFasa=0x1000, IIFwhitebal=0x2000, IIFview=0x4000,
			IIFowner=0x8000, IIFparams=0x10000, IIFcomments=0x20000,
			IIFquality=0x40000, IIFcrop=0x80000, IIForientation=0x100000,
			IIFbearing=0x200000, IIFflength=0x400000
		} ImgInfoFlag;

/* Image Source Types:  digital camera, film scanner, flatbed scanner,
 *		laser rangefinder, image editor, rendering software
 */
typedef enum	{ ISTdigicam, ISTfilm, ISTflatbed,
		ISTrangefinder, ISTeditor, ISTrender }	ImgSourceType;

/* ImageSourceInfo:
 *  Subrecord to hold source type and details.
 */
typedef struct {
	int		stype;		/* source type (ImgSrcType) */
	char		make[64];	/* manufacturer */
	char		model[64];	/* model */
	char		vers[64];	/* version */
} ImgSourceInfo;

/* Flash modes:  off, on, strobe detected, not detected */
typedef enum	{ IFMoff=0, IFMon=1,
		IFMnostrobe=5, IFMstrobe=7, IFMmask=0x7 }	ImgFlashMode;

/* White balance modes:  auto, daylight, daylight1, 
 *		fluorescent, fluorescent1,
 *		incandescent, overcast, flash, shadow,
 *		Illuminant A, Illuminant B, Illuminant C,
 *		5500K, 5000K, incandescent1,
 *		6500K daylight, 7500K, other
 */
typedef enum	{ 	IWBauto=0, IWBdaylight=1, IWBdaylight1=9,
			IWBfluor=2, IWBfluor1=14,
			IWBincand=3, IWBovercast=10,
			IWBflash=4, IWBshadow=11,
			IWBstilla=17, IWBstillb=18, IWBstillc=19,
			IWBd55=20, IWBd50=23, IWBincand1=24,
			IWBd65=21, IWBd75=22, IWBother=255
		} ImgWhiteBal;

/* Image Orientations:  same as TIFF orientations
 * NOTE:  No assignment implies IOtopleft (i.e., normal scanline ordering).
 */
typedef enum	{	IOtopleft=1, IOtopright=2, IObotright=3,
			IObotleft=4, IOlefttop=5, IOrighttop=6,
			IOrightbot=7, IOleftbot=8
		} ImgOrientation;

/* ImgInfo:
 *  Passed structure for returning image metadata, such as
 *  recorded time and exposure.  The flags indicate which
 *  parameters are set (taken from ImgInfoFlag enum above).
 */
typedef struct {
	long		flags;		/* what is set below */
	ImgSourceInfo	source;		/* image source */
	float		hdensity;	/* horizontal pixels per inch */
	float		hvangle;	/* horizontal view angle (degrees) */
	float		stonits;	/* sample-to-nits (or meters) */
	char		capdate[20];	/* YYYY:MM:DD HH:MM:SS */
	char		gmtdate[20];	/* precise capture date/time in GMT */
	float		latlong[2];	/* degrees N latitude, E longitude */
	float		altitude;	/* meters above sea level */
	float		bearing;	/* sighting degrees East of North */
	ImgFlashMode	flash;		/* flash mode */
	float		exptime;	/* exposure time (seconds) */
	float		aperture;	/* aperture (f-stop) */
	float		asa;		/* sensitivity (ASA) */
	float		flength;	/* focal length (millimeters) */
	float		focus;		/* focus distance (meters) */
	ImgWhiteBal	whitebal;	/* white balance mode */
	int		quality;	/* recorded image quality (0-100) */
	ImgRect		crop;		/* crop rectangle for display */
	ImgOrientation	orientation;	/* image orientation in file */
	char		view[256];	/* view specification (Radiance) */
	char		owner[256];	/* name or user id of owner/originator */
	char		params[2048];	/* other parameters:  "key=value\n..." */
	char		comments[2048];	/* additional comments */
} ImgInfo;

extern const ImgInfo	defImgInfo;

/*  Convenience routines for setting and getting params key-value pairs.
 *  The key string may end in '=' but must not contain any white space.
 *  The value may contain white space but only one newline at the end.
 *  Routines return NULL or 0 on failure.  GetImgInfoParam() returns length.
 */
extern const char *	FindImgInfoParam(const ImgInfo *ip, const char *key);
extern int		SetImgInfoParam(ImgInfo *ip,
					const char *key, const char *val);
extern int		GetImgInfoParam(const ImgInfo *ip,
					const char *key, char val[256]);

#ifdef __cplusplus
}
#endif

#endif /* ! _IMGIO_H_ */
