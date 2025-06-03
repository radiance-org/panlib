/*
 *  pimage.c
 *  panlib
 *
 *  Pancine image rendering routines.
 *
 *  Created by gward on Wed Jun 13 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "system.h"
#include "dmessage.h"
#include "pimage.h"
#include "imgwriter.h"
#include "tiff.h"		/* needed for uint16, etc. */
#include "paths.h"

#ifndef RGB24COLORTRANS	
#define RGB24COLORTRANS		1	/* perform RGB24 color transforms? */
#endif

#ifndef ABS
#define ABS(x)		((x) >= 0 ? (x) : -(x))
#endif

#ifndef true
#define true	1
#define false	0
#endif

#define	FEQ(a,b)	((a) < (b)+1e-4 && (b) < (a)+1e-4)

/* Matching definitions for constant declarations in "imgio.h" */
const int		ImgDataSize[] = { sizeof(uby8), sizeof(uint16),
						sizeof(int32), sizeof(float) };
const int		ImgPixelLen[] = {3, 3, 4, 3, 1, 4, 1, 1, 2, 3, 1};

/* 24-bit sRGB color space */
const ImgColorSpace	ICS_sRGB = {
		IDTubyte, IPFrgb, true, 2.2, 0,
		{{.640, .330}, {.300, .600}, {.150, .060}, {.3127, .3290}}
	};

/* 24-bit Adobe 1998 RGB color space */
const ImgColorSpace	ICS_AdobeRGB = {
		IDTubyte, IPFrgb, true, 2.19921875, 0,
		{{.640, .330}, {.210, .710}, {.150, .060}, {.3127, .3290}}
	};

/* 96-bit linear RGB floating point with Adobe 1998 primaries */
const ImgColorSpace	ICS_RGB98Adobe = {
		IDTfloat, IPFrgb, true, 1, 0,
		{{.640, .330}, {.210, .710}, {.150, .060}, {.3127, .3290}}
	};

/* 96-bit ProPhoto RGB (Reference Output Medium Metric) wide-gamut */
const ImgColorSpace	ICS_ProPhotoRGB = {
		IDTfloat, IPFrgb, true, 1, 0,
		{{.734699, .265301}, {.159597, .840403}, {.036598, .000105},
		{.345704, .358540}}
	};

/* 96-bit linear RGB floating point with P3 primaries & D50 white point */
const ImgColorSpace	ICS_P3 = {
		IDTfloat, IPFrgb, true, 1, 0,
		{{.680, .320}, {.265, .690}, {.150, .060}, {.3140, .3510}}
	};

/* 96-bit linear RGB floating point with BT.2020 primaries */
const ImgColorSpace	ICS_RGB2020 = {
		IDTfloat, IPFrgb, true, 1, 0,
		{{.708, .292}, {.170, .797}, {.131, .046}, {.3127, .3290}}
	};

/* 96-bit linear RGB floating point with CCIR-709 primaries */
const ImgColorSpace	ICS_RGB709 = {
		IDTfloat, IPFrgb, true, 1, 0,
		{{.640, .330}, {.300, .600}, {.150, .060}, {.3127, .3290}}
	};

/* 96-bit logarithmic RGB floating point with CCIR-709 primaries (dB) */
const ImgColorSpace	ICS_RGB709dB = {
		IDTfloat, IPFrgb, true, .1, 1,
		{{.640, .330}, {.300, .600}, {.150, .060}, {.3127, .3290}}
	};
/* 96-bit linear RGB floating point with Radiance primaries and neutral white */
const ImgColorSpace	ICS_RGBradiance = {
		IDTfloat, IPFrgb, true, 1, 0,
		{{.640, .330}, {.290, .600}, {.150, .060}, {1./3., 1./3.}}
	};

/* Sharpened RGB color space for white balancing */
const ImgColorSpace	ICS_SharpRGB = {
		IDTfloat, IPFrgb, true, 1, 0,
		{{.6898,.3206}, {.0736,.9003}, {.1166,.0374}, {1./3.,1./3.}}
	};

/* 96-bit floating point CIE XYZ color space */
const ImgColorSpace	ICS_XYZ = {
		IDTfloat, IPFxyz, true, 1, 0,
		{{1, 0}, {0, 1}, {0, 0}, {1./3., 1./3.}}
	};

/* 96-bit floating point CIE XYZ color space (dB) */
const ImgColorSpace	ICS_XYZdB = {
		IDTfloat, IPFxyz, true, .1, 1,
		{{1, 0}, {0, 1}, {0, 0}, {1./3., 1./3.}}
	};

/* 8-bit grayscale */
const ImgColorSpace	ICS_Y8 = {
		IDTubyte, IPFy, true, 2.2, 0,
		{{1./3., 1./3.}}
	};

/* 16-bit linear grayscale */
const ImgColorSpace	ICS_Y16 = {
		IDTushort, IPFy, true, 1, 0,
		{{1./3., 1./3.}}
	};

/* 16-bit logarithmic grayscale matching LogLuv encoding */
const ImgColorSpace	ICS_LogL16 = {
		IDTushort, IPFy, true, 77.0636788, 5.42835481e-20,
		{{1./3., 1./3.}}
};

/* 32-bit floating point grayscale */
const ImgColorSpace	ICS_Y = {
		IDTfloat, IPFy, true, 1, 0,
		{{1./3., 1./3.}}
	};

/* 32-bit floating-point log10(Y)/10 */
const ImgColorSpace	ICS_YdB = {
		IDTfloat, IPFy, true, .1, 1,
		{{1./3., 1./3.}}
	};

/* 32-bit floating-point linear alpha */
const ImgColorSpace	ICS_Alpha = {
		IDTfloat, IPFa, true, 1, 0,
		{{1./3., 1./3.}}
	};

/* 32-bit floating-point depth map */
const ImgColorSpace	ICS_D = {
		IDTfloat, IPFd, true, 1, 0
	};
	
/* 3-D floating-point vector map */
const ImgColorSpace	ICS_VEC3 = {
		IDTfloat, IPFvec3, true, 1, 0
	};

/* 32-bit encoded surface normal */
const ImgColorSpace	ICS_ENORM = {
		IDTint, IPFn, true, 1, 0
	};

/* Zero in any (non-log) color space */
const PixelVal		Pblack;

/* Default image information data */
const ImgInfo		defImgInfo = {
	0L,
	{0, "", "", ""},
	72.f,
	45.f,
	1.f,
	"",
	"",
	{.0f, .0f},
	.0f,
	.0f,
	0,
	1.f,
	1.f,
	100.f,
	50.f,
	1.f,
	0,
	100,
	{0,0,0,0},
	IOtopleft,
	"",
	"",
	"",
	"",
};

/* ImgReader error messages and classes for DMESG(cls,msg) */
static const char *	IREMsg[] = {
	"No error",
	"Unsupported feature",
	"Truncated file",
	"Format error",
	"Read error",
	"Out of memory",
	"Unknown error"
};

static const DMsgClass	IREClass[] = {
	DMCnerr,
	DMCparameter,
	DMCdata,
	DMCdata,
	DMCsystem,
	DMCmemory,
	DMCparameter
};

/* The image reader interface table must be filled explicitly */
const ImgReaderInterface *	PIReaderI[P_MAXREADER+1];

/* Add reader interface to list */
int
PaddIReaderI(const ImgReaderInterface *iri)
{
	int	i;
	if (iri == NULL)
		return -1;
	for (i = 0; PIReaderI[i] != NULL; i++)
		if (PIReaderI[i] == iri)
			return i;
	if (i >= P_MAXREADER)
		return -1;		/* too many interfaces */
	PIReaderI[i] = iri;
	PIReaderI[i+1] = NULL;
	return i;
}

/* Target image readers based on file suffix, returning number of candidates */
int
PtargetReader(const ImgReaderInterface *irilist[P_MAXREADER], const char *sfx)
{
	int	i, n = 0;
	if (irilist == NULL)
		return 0;
					/* match suffix to reader(s) */
	for (i = 0; PIReaderI[i] != NULL; i++)
		if (PmatchSuffix(sfx, PIReaderI[i]->suffixList))
			irilist[n++] = PIReaderI[i];
	return n;
}

/* Open an image file and return reader struct */
ImgReader *
PopenImageF(const char *fname, int quiet, const ImgReaderInterface *iri)
{
	int				noError = true;
	const ImgReaderInterface *	irilist[P_MAXREADER];
	int				nri;
	ImgReader *			ir;

	if (fname == NULL)
		return NULL;
	if (access(fname, R_OK) < 0)	/* check for readable file, first */
		goto fail;
					/* get reader list */
	if (iri != NULL) {		/* use given reader or check suffix? */
		irilist[0] = iri;
		nri = 1;
	} else if (!(nri = PtargetReader(irilist, PgetSuffix(fname)))) {
		noError = false;
		if (!quiet)
			DMESGF(DMCparameter, "Unsupported image format '%s'", fname);
	}
	while (nri-- > 0) {		/* try each reader interface */
		iri = irilist[nri];
		ir = (*iri->Open)(fname);
		if (ir == NULL)
			continue;	/* wrong format? */
		if (ir->errCode == IREnone)
			return ir;	/* successful open! */
		noError = false;
		if (!quiet)
			PreportReaderErr(ir);
		(*iri->Close)(ir);	/* close on error */
	}
fail:
	if (noError & !quiet)
		DMESGF(DMCresource, "Cannot open '%s'", fname);
	return NULL;			/* could not open it */
}

/* Call to make sure sampling is divisible */
void
ImgFixSampling(ImgReadBuf *rb)
{
	if (rb->subsample <= 0)
		rb->subsample = 1;
	if (rb->subsample == 1)
		return;
	rb->r.xleft -= rb->r.xleft % rb->subsample;
	rb->r.xright -= rb->r.xright % rb->subsample;
	if (rb->r.xright <= rb->r.xleft)
		rb->r.xright = rb->r.xleft + rb->subsample;
	rb->r.ytop -= rb->r.ytop % rb->subsample;
	rb->r.ybottom -= rb->r.ybottom % rb->subsample;
	if (rb->r.ybottom <= rb->r.ytop)
		rb->r.ybottom = rb->r.ytop + rb->subsample;
}

/* Assign rectangle based on inclusive corners */
ImgRect *
PsetRect(ImgRect *dr, int x1, int y1, int x2, int y2)
{
	if (dr == NULL)
		return NULL;
	if (x1 < x2) {
		dr->xleft = x1; dr->xright = x2+1;
	} else {
		dr->xleft = x2; dr->xright = x1+1;
	}
	if (y1 < y2) {
		dr->ytop = y1; dr->ybottom = y2+1;
	} else {
		dr->ytop = y2; dr->ybottom = y1+1;
	}
	return dr;
}

/* Clip one rectangle to boundaries of another, returning false if nothing left */
int
PclipRect(ImgRect *dr, const ImgRect *cr)
{
	if ((dr == NULL) | (cr == NULL))
		return false;
	if (dr->xleft < cr->xleft)
		dr->xleft = cr->xleft;
	if (dr->ytop < cr->ytop)
		dr->ytop = cr->ytop;
	if (dr->xright > cr->xright)
		dr->xright = cr->xright;
	if (dr->ybottom > cr->ybottom)
		dr->ybottom = cr->ybottom;
	return (dr->xleft < dr->xright) & (dr->ytop < dr->ybottom);
}

/* Find key-value pair in info parameter string (key may end with '=') */
const char *
FindImgInfoParam(const ImgInfo *ip, const char key[])
{
	const char	*cp, *cp1;
	
	if ((ip == NULL) | (key == NULL))
		return NULL;
	if (!(ip->flags & IIFparams) | !key[0])
		return NULL;
	cp = ip->params;
	while (*cp) {				/* search for key */
		const char	*cp0 = cp;
		for (cp1 = key; *cp1; cp1++) {
			if (*cp1 != *cp)
				break;
			if (*cp1 == '=')	/* key includes '=' ? */
				return cp0;
			if (isspace(*cp1))	/* spaces not allowed */
				return NULL;
			if (!*++cp)		/* unexpected end? */
				return NULL;
		}
		if (!*cp1 & (*cp == '='))	/* matched key? */
			return cp0;
		while (*cp && *cp++ != '\n')	/* else skip to next */
			;
	}
	return NULL;
}

/* Convenience routine for setting params key-value pair */
int
SetImgInfoParam(ImgInfo *ip, const char *key, const char *val)
{
	char	*cp;
	int	vlen, klen;

	if ((ip == NULL) | (key == NULL))
		return false;
						/* verify key */
	for (klen = 0; (key[klen] != '\0') & (key[klen] != '='); klen++)
		if (isspace(key[klen]))
			return false;		/* white space not allowed */
	if (!klen)
		return false;
		
	if (!(ip->flags & IIFparams)) {
		ip->params[0] = '\0';		/* paranoia */
	} else {
		while ((cp = (char *)FindImgInfoParam(ip, key)) != NULL) {
			char	*cp1 = cp;	/* remove previous setting */
			while (*cp1 && *cp1++ != '\n')
				;
			memmove(cp, cp1, strlen(cp1)+1);
		}
		if (!*ip->params)
			ip->flags &= ~IIFparams;
	}
	if (val == NULL) val = "";
						/* get value length */
	for (vlen = 0; (val[vlen] != '\0') & (val[vlen] != '\n'); vlen++)
		if (vlen >= 256)
			return false;
	if (!vlen)				/* unset parameter? */
		return true;
	cp = ip->params;			/* append new value if room */
	while (*cp) cp++;
	if (klen+vlen >= ip->params+(sizeof(ip->params)-3) - cp)
		return false;			/* won't fit! */
	while (klen--)
		*cp++ = *key++;			/* copy key */
	*cp++ = '=';
	while (vlen--)
		*cp++ = *val++;			/* copy value */
	*cp++ = '\n';
	*cp = '\0';
	ip->flags |= IIFparams;			/* made new setting */
	return true;
}

/* Convenience routine for getting params key-value pair */
int
GetImgInfoParam(const ImgInfo *ip, const char *key, char val[256])
{
	const char	*cp;
	char		*cp1;

	if (val == NULL)
		return 0;
	val[0] = '\0';
	cp = FindImgInfoParam(ip, key);		/* index "key=value\n" */
	if (cp == NULL)
		return 0;
	while (*cp++ != '=')
		;
	for (cp1 = val; (cp1-val < 255) & (*cp != '\n'); cp1++)
		if (!(*cp1 = *cp++))
			return cp1 - val;
	*cp1 = '\0';
	return cp1 - val;
}

/* Compute absolute frame from current frame, # frames, animation, and seek value */
int
PabsFrame(int *fmp, int nf, ImgFrameType ft, int offs, ImgSeekMode sm)
{
	int	goal;
	if (sm == IRSabs)			/* absolute positioning? */
		goal = offs;
	else
		goal = *fmp + offs;
	if (goal == *fmp)
		return true;
	if (nf <= 0)
		return false;
						/* no animation? */
	if (sm != IRSadv || (ft != IRFloop) & (ft != IRFbackforth)) {
		if (goal < 0) {
			*fmp = 0;
			return false;
		}
		if (goal >= nf) {
			*fmp = nf - 1;
			return false;
		}
		*fmp = goal;
		return true;
	}
	if (ft == IRFloop) {			/* loop animation */
		if (goal < 0)
			goal = nf - (-goal % nf);
		else
			goal = goal % nf;
		*fmp = goal;
		return true;
	}
						/* back-and-forth */
	if (goal < 0)
		goal = -goal;
	goal = goal % (2*nf);
	if (goal >= nf)
		goal = nf - goal;
	*fmp = goal;
	return true;
}

/* Pancine reader error reporting function */
void
PreportReaderErr(const ImgReader *ir)
{
	if (ir == NULL || ir->errCode == IREnone)
		return;
	if (ir->errCode == IREunknown)
		sprintf(dmessage_buf, "%s: %s", ir->file, ir->errMsg);
	else
		sprintf(dmessage_buf, "%s: %s: %s", ir->file,
				IREMsg[ir->errCode], ir->errMsg);
	dmessage(IREClass[ir->errCode], dmessage_buf, NULL, 0);
}

/* Return pointer to final component "fname.ext" in "/directory/fname.ext" */
const char *
PgetFilename(const char *path)
{
	const char *	cp = path;
	while (*cp)
		cp++;
	while (cp > path && !IS_DIRSEP(cp[-1]))
		cp--;
	return cp;
}

/* Check if suffix matches the given suffix list (e.g., "TIFF.tif") */
int
PmatchSuffix(const char *sfx, const char *sfx_list)
{
	const char *	cp1;

	if ((sfx_list == NULL) | (sfx == NULL) || !*sfx)
		return false;		/* don't attempt w/o suffix */

	for (cp1 = sfx_list; *cp1; cp1++) {
		const char *	cp2;
		for (cp2 = sfx; *cp2; cp2++, cp1++)
			if (tolower(*cp1) != tolower(*cp2))
				break;
		if (!*cp2 && !*cp1 | (*cp1 == '.'))
			return true;
		while (*cp1 != '.')
			if (!*cp1++)
				return false;
	}
	return false;
}

/* Return pointer to suffix in fname, or NULL if none */
const char *
PgetSuffix(const char *fname)
{
	const char *	sfx = fname;
	if (sfx == NULL)
		return "NO FILE";
	while (*sfx)
		sfx++;
	while (sfx > fname && sfx[-1] != '.') {
		--sfx;
		if (IS_DIRSEP(*sfx))
			return NULL;
	}
	if (!*sfx)
		return NULL;
	if (sfx <= fname+1 || IS_DIRSEP(sfx[-2]))
		return NULL;
	return sfx;
}

/* Check pixel values for equality */
int
PequalVal(PixelVal p1, PixelVal p2)
{
	int	p1len;
	int	i;

	if (p1.csp == NULL)
		return (p2.csp == NULL || (p2.csp->logorig == 0 &&
			!memcmp(Pblack.v.b, p2.v.b, ImgPixelSize(p2.csp))));

	p1len = ImgPixelLen[p1.csp->format];

	if (p2.csp == NULL)
		return (p1.csp->logorig == 0 &&
			!memcmp(Pblack.v.b, p1.v.b, p1len*ImgDataSize[p1.csp->dtype]));

	if (p1.csp != p2.csp) {
		const int	p2len = ImgPixelLen[p2.csp->format];
		if (p1len > p2len || (p1len == p2len &&
					p1.csp->dtype > p2.csp->dtype)) {
			p2 = PconvertPixel(p2, p1.csp);
		} else {
			p1 = PconvertPixel(p1, p2.csp);
			p1len = p2len;
		}
	}
	if (p1.csp->dtype != IDTfloat)
		return !memcmp(p1.v.b, p2.v.b, p1len*ImgDataSize[p1.csp->dtype]);

	for (i = p1len; i--; )
		if (!FEQ(p1.v.f[i], p2.v.f[i]))
			return false;
	return true;
}

/* Convenient routine to assign 24-bit sRGB pixel */
PixelVal
PsRGB(int r, int g, int b)
{
	PixelVal	p;
	
	p.csp = &ICS_sRGB;
	p.v.b[0] = (r < 0) ? 0 : (r > 255) ? 255 : r;
	p.v.b[1] = (g < 0) ? 0 : (g > 255) ? 255 : g;
	p.v.b[2] = (b < 0) ? 0 : (b > 255) ? 255 : b;
	
	if (!p.v.b[1] && !p.v.b[0] & !p.v.b[2])
		return Pblack;

	return p;
}

/* Convenient routine to assign floating-point XYZ pixel */
PixelVal
PXYZ(float X, float Y, float Z)
{
	PixelVal	p;
	
	if (Y == .0f)
		return Pblack;

	p.csp = &ICS_XYZ;
	p.v.f[0] = X;
	p.v.f[1] = Y;
	p.v.f[2] = Z;

	return p;
}

/* Convenient routine to assign floating-point Y pixel */
PixelVal
PgraY(float Y)
{
	PixelVal	p;

	if (Y == .0f)
		return Pblack;

	p.csp = &ICS_Y;
	p.v.f[0] = Y;

	return p;
}

/* Extract a single pixel from the given image */
PixelVal
PgetPixel(const ImgStruct *ia, int x, int y)
{
	PixelVal	p;
	int		n;
	const uby8 *	bp;
	
	if (ia == NULL || (ia->img == NULL) | (ia->csp == NULL) ||
			(x < 0) | (x >= ia->xres) || (y < 0) | (y >= ia->yres))
		return Pblack;
	p.csp = ia->csp;
	n = ImgPixelSize(p.csp);
	bp = PpixP(ia, x, y, n) + n;
	while (n--)
		p.v.b[n] = *--bp;
	return p;
}

/* Assign a single pixel in the given image */
int
PsetPixel(ImgStruct *ib, int x, int y, PixelVal p)
{
	int		n;
	uby8 *		bp;

	if (ib == NULL || (ib->img == NULL) | (ib->csp == NULL))
		return false;
	if ((x < 0) | (x >= ib->xres) || (y < 0) | (y >= ib->yres))
		return false;
	if (p.csp != ib->csp)
		p = PconvertPixel(p, ib->csp);
	n = ImgPixelSize(ib->csp);
	bp = PpixP(ib, x, y, n) + n;
	while (n--)
		*--bp = p.v.b[n];
	return true;
}

/* Struct and callback for image averaging */
typedef struct {
	double			sum[MaxPixelLen];
	double			wsum;
	int			ncomp;
	PweightingMethod	*wcb;
	void			*udp;
} ImgAvg;

static int
avgSL(const uby8 *scan, int len, int y, void *dp)
{
	ImgAvg *	iap = (ImgAvg *)dp;
	const float *	fsp = (const float *)scan;
	int		i, x;

	for (x = 0; x < len; x++) {
		double	wt = 1;
		if (iap->wcb != NULL)
			wt = (*iap->wcb)(x, y, iap->udp);
		for (i = 0; i < iap->ncomp; i++)
			iap->sum[i] += *fsp++ * wt;
		iap->wsum += wt;
	}
	return 0;
}

/* Compute weighted average pixel value (in specified CS) from image ia */
PixelVal
PweightedAverageCB(const ImgStruct *ia, const ImgColorSpace *csp,
					PweightingMethod *wcb, void *udp)
{
	ImgColorSpace	fltCS;
	ImgAvg		avg;
	PixelVal	res;

	if (ia == NULL || (ia->img == NULL) | (ia->csp == NULL))
		return Pblack;
	if (csp == NULL)
		csp = ia->csp;		/* use image color space as default */

	PcopyCS(&fltCS, csp);
	fltCS.dtype = IDTfloat;		/* always average in floating point */
	if (fltCS.logorig <= 0)
		fltCS.gamma = 1;	/* force linear averaging if not log */
	memset(&avg, 0, sizeof(avg));
	avg.ncomp = ImgPixelLen[fltCS.format];
	avg.wcb = wcb;
	avg.udp = udp;

	PconvertImageCB(ia, &fltCS, 1.f, avgSL, &avg);

	if (avg.wsum == 0)
		return Pblack;

	res.csp = &fltCS;		/* compute result */
	while (avg.ncomp--)
		res.v.f[avg.ncomp] = (float)(avg.sum[avg.ncomp]/avg.wsum);

	return PconvertPixel(res, csp);	/* convert to destination color space */
}

/* Struct used for image weighting callback */
typedef struct {
	const ImgStruct	*ib;
	void		*bconv;
	float		*bscan;
	int		cury;
} ImgWeightScan;

/* Callback for image weighting (assumes working a scanline at a time) */
static double
grayPixV(int x, int y, void *udp)
{
	ImgWeightScan	*iwp = (ImgWeightScan *)udp;

	if (iwp->cury != y) {
		if (iwp->bconv != NULL)
			PmapPixels((uby8 *)iwp->bscan, ProwPtr(iwp->ib,y),
					iwp->ib->xres, iwp->bconv);
		else
			iwp->bscan = (float *)ProwPtr(iwp->ib,y);
		iwp->cury = y;
	}
	return (double)iwp->bscan[x];
}

/* Compute average of image ia weighted by gray values from ib */
PixelVal
PweightedAverage(const ImgStruct *ia, const ImgStruct *ib,
					const ImgColorSpace *csp)
{
	ImgWeightScan	iws;
	PixelVal	pres;

	if ((ia == NULL) | (ib == NULL) ||
			(ia->csp == NULL) | (ib->csp == NULL))
		return Pblack;
	if ((ia->xres != ib->xres) | (ia->yres != ib->yres)) {
		DMESG(DMCparameter, "Weight image dimensions do not match");
		return Pblack;
	}
	iws.ib = ib;			/* set up call */
	if (PmatchColorSpace(ib->csp, &ICS_Alpha, PICMequiv)) {
		iws.bconv = NULL;	/* no conversion necessary */
		iws.bscan = NULL;
	} else {
		iws.bconv = PcreateColorConv(&ICS_Y, ib->csp, 1.f, 1);
		if (iws.bconv == NULL) {
			DMESG(DMCparameter, "Cannot convert weight image to gray");
			return Pblack;
		}
		iws.bscan = (float *)malloc(sizeof(float)*ib->xres);
		if (iws.bscan == NULL) {
			DMESG(DMCmemory, "malloc() failed in PweightedAverage");
			PfreeColorConv(iws.bconv);
			return Pblack;
		}
	}
	iws.cury = -1;			/* make the call */
	pres = PweightedAverageCB(ia, csp, grayPixV, &iws);
	if (iws.bconv != NULL) {
		free(iws.bscan);
		PfreeColorConv(iws.bconv);
	}
	return pres;			/* return after cleanup */
}

/* Return string with description of color space */
char *
PdescribeCS(const ImgColorSpace *csp, char *str)
{
	static const char	typNam[][16] = {"8-bit/sample", "16-bit/sample",
					"32-bit/sample", "32-bit/float"};
	static const char	fmtNam[][8] = {"RGB", "CMY", "CMYK", "XYZ", "gray",
					"RGBA", "alpha", "depth", "UV",
					"vector", "normal"};
	static char		myBuf[64];
	char			*s;

	if (str == NULL) str = myBuf;
	if (csp == NULL)
		return strcpy(str, "NULL colorspace");
	if ((csp->dtype < 0) | (csp->dtype >= 4) |
			(csp->format < 0) | (csp->format >= 11))
		return strcpy(str, "BAD colorspace");
	s = str;
	if (csp->logorig > 0) {
		strcpy(s, "log ");
		s += 4;
	} else if (FEQ(csp->gamma, 1.0)) {
		strcpy(s, "linear ");
		s += 7;
	} else {
		sprintf(s, "gamma %.2f ", csp->gamma);
		s += 6;
		while (*s) s++;
	}
	strcpy(s, typNam[csp->dtype]);
	while (*s) s++;
	*s++ = ' ';
	strcpy(s, fmtNam[csp->format]);

	return str;
}

/* Check two color spaces to see if they match */
int
PmatchColorSpace(const ImgColorSpace *cs1, const ImgColorSpace *cs2, int cm)
{
	if ((cs1 == NULL) | (cs2 == NULL))
		return false;
	if (cs1 == cs2)
		return true;
	if (cm == PICMequiv) {			/* special case for equivalent CS */
		cm = PICMdtype|PICMgamma;
		cm |= ((ImgPixelLen[cs1->format] > 1) | (ImgPixelLen[cs2->format] > 1)) *
				(PICMformat|PICMprims);
	}
	if (cm&PICMdtype && cs1->dtype != cs2->dtype)
		return false;
	if (cm&(PICMchroma|PICMwhite) &&
			!((cs1->format == IPFrgb) & (cs2->format == IPFrgba)) &&
			!((cs1->format == IPFrgba) & (cs2->format == IPFrgb)))
		cm |= PICMformat;
	if (cm&PICMformat && cs1->format != cs2->format)
		return false;
	if (cm&PICMgamma && !FEQ(cs1->gamma, cs2->gamma) |
				(cs1->logorig != cs2->logorig))
		return false;
	if (cm&PICMchroma) {
		switch (cs1->format) {
		case IPFcmyk:
			if (!FEQ(cs1->chroma[3][0], cs2->chroma[3][0]) ||
				!FEQ(cs1->chroma[3][1], cs2->chroma[3][1]))
				return false;
		/* FALL THROUGH */
		case IPFrgb:
		case IPFrgba:
		case IPFcmy:
			if (!FEQ(cs1->chroma[2][0], cs2->chroma[2][0]) ||
				!FEQ(cs1->chroma[2][1], cs2->chroma[2][1]))
				return false;
			if (!FEQ(cs1->chroma[1][0], cs2->chroma[1][0]) ||
				!FEQ(cs1->chroma[1][1], cs2->chroma[1][1]))
				return false;
			if (!FEQ(cs1->chroma[0][0], cs2->chroma[0][0]) ||
				!FEQ(cs1->chroma[0][1], cs2->chroma[0][1]))
				return false;
			break;
		default:
			break;
		}
	}
	if (cm&PICMwhite)
		switch (cs1->format) {
		case IPFrgb:
		case IPFrgba:
		case IPFcmy:
			if (!FEQ(cs1->chroma[3][0], cs2->chroma[3][0]) ||
				!FEQ(cs1->chroma[3][1], cs2->chroma[3][1]))
				return false;
			break;
		/* white point is meaningless for gray images...
		case IPFy:
			if (!FEQ(cs1->chroma[0][0], cs2->chroma[0][0]) ||
				!FEQ(cs1->chroma[0][1], cs2->chroma[0][1]))
				return false;
			break;
		*/
		default:		/* we consider XYZ to be absolute */
			break;
		}
	return true;
}

/* Convert from 8-bit gray to 24-bit RGB (basically a copy function) */
void
PgetRGB24fromY8(uby8 *dptr, const uby8 *sptr, int n)
{
	sptr += n;			/* go from end in case dptr == sptr */
	dptr += n*3;
	while (n-- > 0) {
		*--dptr = *--sptr;
		*--dptr = *sptr;
		*--dptr = *sptr;
	}
}

/* Convert from 32-bit RGBA to 24-bit RGB (copy function) */
void
PgetRGB24fromRGBA32(uby8 *dptr, const uby8 *sptr, int n)
{
	while (n-- > 0) {
		*dptr++ = *sptr++;
		*dptr++ = *sptr++;
		*dptr++ = *sptr++;
		sptr++;
	}
}

/* Add alpha channel to 24-bit RGB (copy function) */
void
PgetRGBA32fromRGB24(uby8 *dptr, const uby8 *sptr, int n)
{
	sptr += n*3;			/* go from end in case dptr == sptr */
	dptr += n*4;
	while (n-- > 0) {
		*--dptr = 255;
		*--dptr = *--sptr;
		*--dptr = *--sptr;
		*--dptr = *--sptr;
	}
}

/* Convert from 8-bit gray to 32-bit RGBA (copy function) */
void
PgetRGBA32fromY8(uby8 *dptr, const uby8 *sptr, int n)
{
	sptr += n;			/* go from end in case dptr == sptr */
	dptr += n*4;
	while (n-- > 0) {
		*--dptr = 255;
		*--dptr = *--sptr;
		*--dptr = *sptr;
		*--dptr = *sptr;
	}
}

/* Link subimage -- assumes ib is empty or equals ia */
int
PlinkSubimage(ImgStruct *ib, const ImgStruct *ia, const ImgRect *r)
{
						/* check arguments */
	if ((ib == NULL) | (ia == NULL) || (ia->img == NULL) | (ia->csp == NULL))
		return false;
	if (r == NULL) {			/* link entire image */
		if (ib == ia)
			return true;
		*ib = *ia;
		if (ib->mbase != NULL)
			MobjRetain(ib->mbase);
		return true;
	}
	if (!PlegalRect(r, ia->xres, ia->yres)) {
		ib->img = NULL;
		return false;
	}
						/* else link to region */
	if (ib != ia) {
		ib->csp = ia->csp;
		ib->rowsize = ia->rowsize;
		if ((ib->mbase = ia->mbase) != NULL)
			MobjRetain(ib->mbase);
	}
	ib->img = PpixPtr(ia, r->xleft, r->ytop);
	ib->xres = r->xright - r->xleft;
	ib->yres = r->ybottom - r->ytop;
	return true;
}

/* Link corresponding overlap, placing upper left corner of img2 in img1 */
int
PlinkCoverage(ImgStruct *cover1, ImgStruct *cover2,
			const ImgStruct *img1, const ImgStruct *img2,
			int xleft, int ytop)
{
	ImgRect		rect1, rect2;
						/* check arguments */
	if ((img1 == NULL) | (img2 == NULL))
		return false;
	if ((img1->img == NULL) | (img2->img == NULL))
		return false;
						/* assign image regions */
	rect2.xleft = 0; rect2.xright = img2->xres;
	rect2.ytop = 0; rect2.ybottom = img2->yres;
	rect1.xleft = xleft;
	rect1.xright = xleft + img2->xres;
	rect1.ytop = ytop;
	rect1.ybottom = ytop + img2->yres;
	if (rect1.xleft < 0) {
		rect2.xleft -= rect1.xleft;
		rect1.xleft = 0;
	}
	if (rect1.xright > img1->xres) {
		rect2.xright += img1->xres - rect1.xright;
		rect1.xright = img1->xres;
	}
	if (rect1.ytop < 0) {
		rect2.ytop -= rect1.ytop;
		rect1.ytop = 0;
	}
	if (rect1.ybottom > img1->yres) {
		rect2.ybottom += img1->yres - rect1.ybottom;
		rect1.ybottom = img1->yres;
	}
	if (!PlinkSubimage(cover2, img2, &rect2))
		return false;
	if (!PlinkSubimage(cover1, img1, &rect1)) {
		if (cover2 != img2) PfreeImage(cover2);
		return false;
	}
	return true;
}

/* Get position from pointer */
int
PpixPos(int *xp, int *yp, const ImgStruct *ia, const uby8 *ptr)
{
	ssize_t	offset;

	if ((xp == NULL) | (yp == NULL) | (ia == NULL) | (ptr == NULL))
		return false;
	if ((ia->csp == NULL) | (ia->img == NULL))
		return false;
	offset = ptr - ia->img;
	if (ia->rowsize < 0) {
		if ((offset += ia->rowsize + 1) > 0)
			return false;
	} else if (offset < 0)
		return false;
	*yp = offset / ia->rowsize;
	if (*yp >= ia->yres)
		return false;
	*xp = (ptr - ProwPtr(ia,*yp)) / ImgPixelSize(ia->csp);
	return (*xp < ia->xres);
}

/* Check if image can possibly fit memory */
int
PsizeOK(int xr, int yr, int psiz)
{
	return (log((double)xr*yr*psiz) < 5.54517743*sizeof(size_t));
}

/* Allocate an image buffer for the given color space, adjusting dimensions */
int
PnewImage(ImgStruct *ib, double imgAspect)
{
	int	cspace;
	int	psize;
						/* argument check */
	if (ib == NULL || (ib->xres <= 0) | (ib->yres <= 0) | (ib->csp == NULL))
		return false;
						/* get pixel size */
	psize = ImgPixelSize(ib->csp);
						/* adjust dimensions */
	if (imgAspect > .01) {
		int	nres;
		if (ib->yres >= ib->xres*imgAspect) {
			nres = (int)(ib->xres*imgAspect + .5);
			if (ib->yres > nres + 1)
				ib->yres = nres;
		} else {
			nres = (int)(ib->yres/imgAspect + .5);
			if (ib->xres > nres + 1)
				ib->xres = nres;
		}
	}
	if (ib->img != NULL)			/* assume it's big enough */
		return (ABS(ib->rowsize) >= ib->xres*psize);

	if (!PsizeOK(ib->xres, ib->yres, psize)) {
		sprintf(dmessage_buf, "%dx%d image too big for address space",
				ib->xres, ib->yres);
		DMESG(DMCparameter, dmessage_buf);
		return false;
	}
						/* space for color if needed */
	cspace = !ib->csp->cstatic * sizeof(ImgColorSpace);
						/* get new image buffer */
	ib->rowsize = ib->xres*psize;
	CacheMakeRoom(cspace + (size_t)ib->yres*ib->rowsize);
	ib->mbase = MOalloc(cspace + (size_t)ib->yres*ib->rowsize);
	if (ib->mbase == NULL)
		return false;			/* MOalloc reported error */
	if (cspace) {
		ImgColorSpace *	newcs = (ImgColorSpace *)MobjMem(ib->mbase);
		PcopyCS(newcs, ib->csp);	/* keep our own copy of CS */
		ib->csp = newcs;
		ib->img = (uby8 *)(newcs + 1);
	} else
		ib->img = (uby8 *)MobjMem(ib->mbase);
	return true;				/* ready for image data */
}

/* Make unique copy of image and eliminate dead space */
int
PdelinkImage(ImgStruct *ib)
{
	int		cspace;
	ImgStruct	orig_img;
						/* sanity check */
	if (ib == NULL || ib->img == NULL)
		return false;
						/* temporary holder */
	orig_img = *ib;
						/* will remove dead space */
	ib->rowsize = ib->xres*ImgPixelSize(ib->csp);
						/* special cases */
	if ((orig_img.rowsize < 0) | (ib->mbase == NULL) ||
			(ib->mbase->nref > 1) | (ib->mbase->freo != SysFree)) {
		ib->img = NULL;			/* need fresh allocation */
		if (!PmapImage(ib, &orig_img, 1.f)) {
			*ib = orig_img;		/* insufficient memory? */
			return false;
		}
		PfreeImage(&orig_img);		/* unlink original */
		return true;
	}
						/* private color space? */
	cspace = (ib->csp == (ImgColorSpace *)MobjMem(ib->mbase))*sizeof(ImgColorSpace);
						/* reset image pointer */
	ib->img = (uby8 *)MobjMem(ib->mbase) + cspace;

	if (!PmapImage(ib, &orig_img, 1.f))	/* overlapping copy */
		return false;
						/* reallocate to save space */
	ib->mbase = (MemObject *)realloc(orig_img.mbase, sizeof(MemObject) +
					cspace + (size_t)ib->yres*ib->rowsize);
	if (ib->mbase == NULL) {
		ib->mbase = orig_img.mbase;	/* realloc() screw-up */
		return true;			/* irrelevant to result */
	}
						/* reassign pointers */
	if (cspace)
		ib->csp = (ImgColorSpace *)MobjMem(ib->mbase);
	ib->img = (uby8 *)MobjMem(ib->mbase) + cspace;
	return true;
}

/* Unlink an image buffer */
void
PfreeImage(ImgStruct *ib)
{
	if (ib == NULL || ib->img == NULL)
		return;
	ib->img = NULL;
	if (ib->mbase == NULL)
		return;
	if (ib->csp == (ImgColorSpace *)MobjMem(ib->mbase))
		ib->csp = NULL;			/* color space w/ image */
	MobjRelease(ib->mbase);
	ib->mbase = NULL;
}

/* Determine if two images have overlapping pixel regions */
int
PimagesOverlap(const ImgStruct *ia, const ImgStruct *ib)
{
	ImgStruct	ia_rev, ib_rev, ibase;
	ImgRect		recta, rectb;
	int		psize;

	if ((ia == NULL) | (ib == NULL) ||
			(ia->img == NULL) | (ib->img == NULL) |
			(ia->csp == NULL) | (ib->csp == NULL))
		return false;
	if (ia->img == ib->img)
		return true;
	if ((ia->mbase == NULL) ^ (ib->mbase == NULL) || ia->mbase != ib->mbase)
		return false;			/* reasonable assumption #1 */
	psize = ImgPixelSize(ia->csp);
	if (psize != ImgPixelSize(ib->csp))
		return false;			/* reasonable assumption #2 */
	if (ia->rowsize < 0) {
		ia_rev = *ia;			/* need to flip image 'a' */
		ia_rev.img = ia->img + (ssize_t)(ia->yres - 1)*ia->rowsize;
		ia_rev.rowsize = -ia->rowsize;
		ia = &ia_rev;
	}
	if (ib->rowsize < 0) {
		ib_rev = *ib;			/* need to flip image 'b' */
		ib_rev.img = ib->img + (ssize_t)(ib->yres - 1)*ib->rowsize;
		ib_rev.rowsize = -ib->rowsize;
		ib = &ib_rev;
	}
	if (ia->img < ib->img + (ssize_t)ib->yres*ib->rowsize ^
			ib->img < ia->img + (ssize_t)ia->yres*ia->rowsize)
		return false;			/* scanlines do not overlap */
	if (ia->rowsize != ib->rowsize)
		return (ia->mbase != NULL);	/* reasonable assumption #3 */
	if ((ia->xres + ib->xres)*psize > ia->rowsize)
		return true;			/* columns must overlap */
	if (ia->mbase == NULL || ia->mbase->freo != SysFree)
		return false;			/* XXX not ours -- unknown */
	/*
	 * The following is a bit brute-force, but I gave up on other
	 * approaches because they were too clever for me.  We derive
	 * a facsimile of the original base image, then determine the
	 * subareas covered by our two subimages.  All that's left is
	 * a quick check to see if the subrectangles overlap.
	 */
	ibase = *ia;				/* recreate base image */
	ibase.img = (uby8 *)MobjMem(ibase.mbase);
	if (((const uby8 *)ia->csp == ibase.img) |
			((const uby8 *)ib->csp == ibase.img))
		ibase.img += sizeof(ImgColorSpace);
	ibase.xres = ibase.rowsize / psize;
	if (ProwPtr(ia,ia->yres) > ProwPtr(ib,ib->yres))
		ibase.yres = (ProwPtr(ia,ia->yres) - ibase.img)/ibase.rowsize;
	else
		ibase.yres = (ProwPtr(ib,ib->yres) - ibase.img)/ibase.rowsize;
						/* determine subareas */
	if (!PpixPos(&recta.xleft, &recta.ytop, &ibase, ia->img))
		return false;
	if (!PpixPos(&recta.xright, &recta.ybottom, &ibase,
				PpixP(ia,ia->xres-1,ia->yres-1,psize)))
		return false;
	recta.xright++; recta.ybottom++;
	if (!PpixPos(&rectb.xleft, &rectb.ytop, &ibase, ib->img))
		return false;
	if (!PpixPos(&rectb.xright, &rectb.ybottom, &ibase,
				PpixP(ib,ib->xres-1,ib->yres-1,psize)))
		return false;
	rectb.xright++; rectb.ybottom++;
						/* do rectangles overlap? */
	return PclipRect(&recta, &rectb);
}

/* Assign entire image to the given pixel value */
int
PsetImage(ImgStruct *ib, PixelVal pfv)
{
	if (!PnewImage(ib, .0))			/* make sure we're allocated */
		return false;
	if (pfv.csp == NULL) {			/* recognize black */
		if (ib->csp->logorig > 0)
			DMESG(DMCwarning, "Clearing log image to zeroes");
		PclearImage(ib, NULL);
		return true;
	}
	pfv = PconvertPixel(pfv, ib->csp);
	if (pfv.csp == NULL) {
		DMESG(DMCparameter, "Cannot convert pixel to destination space");
		return false;
	}
	PclearImage(ib, pfv.v.b);
	return true;
}

/* Check if run of bytes are all the same */
static int
isByteRun(const uby8 *bp, int n)
{
	while (n > 1)
		if (bp[--n] != bp[0])
			return 0;
	return 1;
}

/* Clear an image to the given pixel (black if pf==NULL) */
void
PclearImage(ImgStruct *ib, const void *pf)
{
	int		psize;
	int		nrows;
	uby8 *		bp;

	if (ib == NULL || (ib->img == NULL) | (ib->csp == NULL))
		return;
	psize = ImgPixelSize(ib->csp);
	if (ib->rowsize == ib->xres*psize) {
		if (pf == NULL) {
			memset(ib->img, 0, (size_t)ib->yres*ib->rowsize);
			return;
		}
		if (isByteRun((const uby8 *)pf, psize)) {
			memset(ib->img, *(const uby8 *)pf, (size_t)ib->yres*ib->rowsize);
			return;
		}
	}
	for (bp = ib->img, nrows = ib->yres; nrows-- > 0; bp += ib->rowsize) {
		int	ncols;
		if (pf == NULL) {
			memset(bp, 0, ib->xres*psize);
			continue;
		}
		if (bp != ib->img) {
			memcpy(bp, ib->img, ib->xres*psize);
			continue;
		}
		bp += psize * (ncols = ib->xres);
		while (ncols-- > 0)
			memcpy((bp -= psize), pf, psize);
	}
}

/* Clear the indicated rectangle */
void
PclearRect(ImgStruct *ib, const ImgRect *r, void *pf)
{
	ImgStruct	subimg;

	if (PlinkSubimage(&subimg, ib, r)) {
		PclearImage(&subimg, pf);
		PfreeImage(&subimg);
	}
}

/* Shift image contents, filling with the given callback function */
int
PshiftImageCB(ImgStruct *ib, int dx, int dy, PfillRectMethod *fcb, void *udp)
{
	ImgRect		drct;
	ImgStruct	simg, dimg;
	int		ok;

	if (ib == NULL || (ib->img == NULL) | (ib->csp == NULL))
		return false;
	if (!dx & !dy)
		return true;
						/* link src & dest regions */
	if (!PlinkCoverage(&dimg, &simg, ib, ib, dx, dy)) {
		if (fcb == NULL)		/* nothing left */
			return false;
		drct.xleft = 0; drct.xright = ib->xres;
		drct.ytop = 0; drct.ybottom = ib->yres;
		(*fcb)(ib, &drct, udp);
		return true;
	}
						/* perform move */
	ok = PmapImage(&dimg, &simg, 1.f);
	PfreeImage(&simg);
	PfreeImage(&dimg);
	if (!ok | (fcb == NULL))		/* no fill requested? */
		return ok;
						/* fill callback(s) */
	if (dy) {
		drct.xleft = 0;
		drct.xright = ib->xres;
		if (dy > 0) {
			drct.ytop = 0; drct.ybottom = dy;
		} else {
			drct.ytop = ib->yres + dy; drct.ybottom = ib->yres;
		}
		(*fcb)(ib, &drct, udp);
	}
	if (dx) {
		if (dx > 0) {
			drct.xleft = 0; drct.xright = dx;
		} else {
			drct.xleft = ib->xres + dx; drct.xright = ib->xres;
		}
		if (dy > 0) {
			drct.ytop = dy; drct.ybottom = ib->yres;
		} else {
			drct.ytop = 0; drct.ybottom = ib->yres + dy;
		}
		(*fcb)(ib, &drct, udp);
	}
	return true;				/* all done */
}

/* Render (resample) an image from reader */
int
PrenderImageR(ImgStruct *ib, ImgReader *ir, int quiet)
{
	uby8 *		mybuf = NULL;
	ImgStruct	isrc;
	ImgReadBuf	rdb;
	int		xread, yread;

	if ((ib == NULL) | (ir == NULL) || ib->csp == NULL)
		return false;
						/* for error reporting */
	ir->errCode = IREnone;
						/* set output resolution */
	if ((ib->xres <= 0) | (ib->yres <= 0)) {
		DMESG(DMCparameter, "Illegal output resolution");
		return false;
	}
						/* allocate destination */
	if (!PnewImage(ib, ir->yres/(ir->xres*(double)ir->pixAspect)))
		return false;
						/* not RGB24? */
	if (!PmatchColorSpace(ib->csp, &ICS_sRGB, PICMptype)) {
		int	csxlate = !PmatchColorSpace(&ir->cs,ib->csp,PICMall);
		PcopyCS(&rdb.cs, ib->csp);
	reread1:
		rdb.r.xleft = 0; rdb.r.xright = ir->xres;
		rdb.r.ytop = 0; rdb.r.ybottom = ir->yres;
		rdb.subsample = 1;
		if ( (ib->xres == ir->xres) & (ib->yres == ir->yres) &&
				ib->rowsize == ib->xres*ImgPixelSize(&rdb.cs) ) {
			rdb.buf = ib->img;
		} else if ((rdb.buf = mybuf =
				(uby8 *)malloc(ImgReadBufLen(&rdb))) == NULL) {
			DMESG(DMCmemory, "No buffer memory in PrenderImageR()");
			PfreeImage(ib);
			return false;
		}
		if (IRreadRec(ir, &rdb) != IREnone)
			goto fail;
		if ((rdb.buf != ib->img) & (rdb.buf != mybuf)) {
			free(rdb.buf);
			if (mybuf != NULL)
				free(mybuf);
			mybuf = NULL;
			if (csxlate > 0) {
				csxlate = -1;
				goto reread1;
			}
			DMESGF(DMCparameter, "IRreadRec() misbehaving for %s",
						ir->encoding);
			PfreeImage(ib);
			return false;
		}
		if (rdb.buf != ib->img) {	/* copy/resample image */
			isrc.csp = &rdb.cs;
			isrc.xres = ir->xres; isrc.yres = ir->yres;
			isrc.rowsize = ir->xres*ImgPixelSize(&rdb.cs);
			isrc.mbase = NULL; isrc.img = rdb.buf;
			if (!PsizeImage(ib, &isrc, PSbest))
				isrc.img = NULL;
			free(rdb.buf);
			return (isrc.img != NULL);
		}
		if (csxlate < 0) {		/* in situ color space conv. */
			const ImgColorSpace *	cstarget = ib->csp;
			ib->csp = &rdb.cs;	/* as a matter of fact */
			return PconvertColorSpace(ib, cstarget, 1.f);
		}
		return true;			/* direct read ends here */
	}
	/* Tone-map & resample in reader when destination is RGB24... */
						/* set up read buffer */
	if (ir->cs.format == IPFy)		/* if it's grayscale... */
		PcopyCS(&rdb.cs, &ICS_Y8);	/* request 8-bit */
	else if (PmatchColorSpace(&ir->cs, ib->csp, PICMptype))
		PcopyCS(&rdb.cs, &ir->cs);	/* close enough? */
	else
		PcopyCS(&rdb.cs, ib->csp);	/* request full conversion */
	/*
	 *  Here we resort to reading the entire image at a resolution
	 *  of as much as 4 times the final output, then downsample
	 *  (or upsample) using a Gaussian (or bicubic) filter.  It would be
	 *  better to read strips, but I nearly went insane trying to write
	 *  that little piece of code....
	 */
	if (PmatchSuffix("jpg", ir->ri->suffixList)) {
		rdb.subsample = ir->yres/ib->yres;
		if (rdb.subsample >= 8)		/* JPEGlib downsamples */
			rdb.subsample = 8;
		else if (rdb.subsample >= 4)
			rdb.subsample = 4;
		else if (rdb.subsample >= 2)
			rdb.subsample = 2;
		else
			rdb.subsample = 1;
	} else {
		rdb.subsample = ir->pixAspect >= 1.f ? ir->xres/(2*ib->xres)
						: ir->yres/(2*ib->yres) ;
		rdb.subsample += !rdb.subsample;
	}
	rdb.r.xleft = 0;
	rdb.r.xright = ir->xres - (ir->xres % rdb.subsample);
	rdb.r.ytop = 0;
	rdb.r.ybottom = ir->yres - (ir->yres % rdb.subsample);
	xread = rdb.r.xright / rdb.subsample;
	yread = rdb.r.ybottom / rdb.subsample;
	if (ib->xres == xread && ib->yres == yread && ib->rowsize == 3*xread) {
		rdb.buf = ib->img;		/* read directly if we can */
	} else {				/* else buffer the read */
		mybuf = (uby8 *)malloc(yread*xread*3);
		if (mybuf == NULL) {
			DMESG(DMCmemory, "No buffer memory in PrenderImageR()");
			goto fail;
		}
		rdb.buf = mybuf;
	}
	if (IRreadRec(ir, &rdb) != IREnone)	/* read (and subsample) */
		goto fail;
	if ((rdb.buf != ib->img) & (rdb.buf != mybuf)) {
		DTESTF(!quiet, DMCparameter, "Unsupported image conversion for %s",
					ir->encoding);
		free(rdb.buf);
		goto fail;
	}
	if (rdb.cs.format == IPFy) {		/* convert from grayscale */
		PgetRGB24fromY8(rdb.buf, rdb.buf, xread*yread);
		PcopyCS(&rdb.cs, &ICS_sRGB);
	}
	if (mybuf != NULL) {			/* resample result? */
		isrc.csp = &rdb.cs;
		isrc.xres = xread; isrc.yres = yread;
		isrc.rowsize = xread*3;
		isrc.mbase = NULL; isrc.img = mybuf;
		if (!PsizeImage(ib, &isrc, PSbest))
			goto fail;
		free(mybuf);
	} else if (!PmatchColorSpace(&rdb.cs, ib->csp, PICMall)) {
#if RGB24COLORTRANS				/* final color conversion? */
		const ImgColorSpace *	cstarget = ib->csp;
		ib->csp = &rdb.cs;		/* as a matter of fact */
		if (!PconvertColorSpace(ib, cstarget, 1.f))
			goto fail;
#else
		DTESTF(!quiet, DMCwarning, "Unconverted color space for %s",
					ir->encoding);
#endif
	}
	return true;				/* all done! */
fail:
	if (!quiet)
		PreportReaderErr(ir);
	if (mybuf != NULL)
		free(mybuf);
	PfreeImage(ib);
	return false;
}

/* Render (resample) an image from another image */
int
PrenderImageI(ImgStruct *ib, const ImgStruct *ia)
{
	if ((ib == NULL) | (ia == NULL))	/* check & allocate dest. */
		return false;
	if ((ia->xres <= 0) | (ia->yres <= 0))
		return false;
	if ((ib->xres <= 0) | (ib->yres <= 0)) {
		DMESG(DMCparameter, "Illegal output resolution");
		return false;
	}
	if (!PnewImage(ib, (double)ia->yres/ia->xres))
		return false;
						/* filter/copy image */
	return PsizeImage(ib, ia, PSbest);
}

/* Copy one image into another at indicated position */
int
PcopyImage(ImgStruct *ib, const ImgStruct *ia, int xleft, int ytop)
{
	ImgStruct	srcImg, dstImg;
	int		ok;
						/* link our overlap */
	if (!PlinkCoverage(&dstImg, &srcImg, ib, ia, xleft, ytop))
		return false;
						/* copy/convert pixels */
	ok = PmapImage(&dstImg, &srcImg, 1.f);
						/* clean up and return */
	PfreeImage(&srcImg);
	PfreeImage(&dstImg);
	return ok;
}

/* Extract alpha channel from RGBA image */
int
PseparateAlpha(ImgStruct *ib, ImgStruct *ialpha, const ImgStruct *ia)
{
	ImgStruct	tempImg, *dimp;
	ImgColorSpace	mySpace;
	int		ok;

	if (ia == NULL || ia->img == NULL || ia->csp->format != IPFrgba)
		return false;
	tempImg.img = NULL;
	if (ib != NULL) {			/* copy to RGB destination image */
		const int	dsiz = ImgDataSize[ia->csp->dtype];
		int		y = ia->yres;
		dimp = ib;			/* check destination space */
		if (ib->csp->format != IPFrgb) {
			DMESG(DMCparameter, "Bad RGB destination in PseparateAlpha()");
			return false;
		}
		PcopyCS(&mySpace, ia->csp);
		mySpace.format = IPFrgb;
		if (!PmatchColorSpace(ib->csp, &mySpace, PICMall)) {
			tempImg.csp = &mySpace;
			dimp = &tempImg;
		}
		if (dimp->img == NULL) {	/* fresh destination? */
			dimp->xres = ia->xres;
			dimp->yres = ia->yres;
			if (!PnewImage(dimp, .0))
				return false;
		} else if ((dimp->xres != ia->xres) | (dimp->yres != ia->yres)) {
			DMESG(DMCparameter, "Image size mismatch in PseparateAlpha()");
			return false;
		}
		if (dsiz == 1)
			while (y-- > 0)
				PgetRGB24fromRGBA32(ProwPtr(dimp,y), ProwPtr(ia,y), ia->xres);
		else
			while (y-- > 0) {
				int		n = ia->xres;
				uby8		*dp = ProwPtr(dimp,y);
				const uby8	*sp = ProwPtr(ia,y);
				while (n-- > 0) {
					memcpy(dp, sp, 3*dsiz);
					dp += 3*dsiz;
					sp += 4*dsiz;
				}
			}
		if (dimp == &tempImg) {		/* convert to destination space? */
			ok = PmapImage(ib, &tempImg, 1.f);
			PfreeImage(&tempImg);
			if (!ok) return false;
		}
	}
	if (ialpha == NULL)			/* discarding alpha? */
		return true;
	if (ialpha->csp == NULL || ImgPixelLen[ialpha->csp->format] != 1) {
		DMESG(DMCparameter, "Need single-channel alpha destination in PseparateAlpha()");
		PfreeImage(ib);
		return false;
	}
	dimp = ialpha;				/* check color space */
	PcopyCS(&mySpace, ia->csp);
	mySpace.format = ialpha->csp->format;
	if (!PmatchColorSpace(ialpha->csp, &mySpace, PICMall)) {
		tempImg.csp = &mySpace;
		dimp = &tempImg;
	}
	ok = PcopyComponent(dimp, 0, ia, 3);	/* copy alpha channel */
	if (ok && dimp == &tempImg)		/* convert to destination space? */
		ok = PmapImage(ialpha, &tempImg, 1.f);
	PfreeImage(&tempImg);
	if (!ok) PfreeImage(ib);
	return ok;
}

/* Integrate alpha channel into RGBA destination */
int
PmarryAlpha(ImgStruct *ib, const ImgStruct *ia, const ImgStruct *ialpha)
{
	ImgStruct	tempImg;
	ImgColorSpace	mySpace;
	int		y, dsiz;
	
	if (ia == NULL || ia->img == NULL || ia->csp->format != IPFrgb)
		return false;
	if (ib == NULL || ib->csp == NULL || ib->csp->format != IPFrgba)
		return false;
	tempImg.img = NULL;			/* need to convert input color? */
	PcopyCS(&mySpace, ib->csp);
	mySpace.format = IPFrgb;
	if (!PmatchColorSpace(&mySpace, ia->csp, PICMall)) {
		tempImg.csp = &mySpace;
		if (!PmapImage(&tempImg, ia, 1.f))
			return false;
		ia = &tempImg;
	}
	if (ib->img == NULL) {
		ib->xres = ia->xres;
		ib->yres = ia->yres;
	} else if ((ib->xres != ia->xres) | (ib->yres != ia->yres)) {
		DMESG(DMCparameter, "Image resolution mismatch in PmarryAlpha()");
		return false;
	}
	if (ialpha != NULL && ialpha->img != NULL) {
		ImgStruct	myAlpha;	/* use supplied alpha channel */
		myAlpha.img = NULL;
		if (ialpha->csp == NULL || ImgPixelLen[ialpha->csp->format] != 1) {
			DMESG(DMCparameter,
				"Source alpha must be single channel in PmarryAlpha()");
			return false;
		}
		PcopyCS(&mySpace, ib->csp);
		mySpace.format = ialpha->csp->format;
		if (!PmatchColorSpace(&mySpace, ialpha->csp, PICMall)) {
			myAlpha.csp = &mySpace;
			if (!PmapImage(&myAlpha, ialpha, 1.f))
				return false;
			ialpha = &myAlpha;
		}
		y = PcopyComponent(ib, 3, ialpha, 0);
		PfreeImage(&myAlpha);
		if (!y) return false;
	} else {
		PixelVal	alphaPix;	/* else assign default alpha */
		alphaPix.csp = ib->csp;
		switch (alphaPix.csp->dtype) {
		case IDTubyte:
			if (!PnewImage(ib, .0))
				return false;
			y = ia->yres;
			while (y-- > 0)
				PgetRGBA32fromRGB24(ProwPtr(ib,y), ProwPtr(ia,y), ia->xres);
			PfreeImage(&tempImg);
			return true;
		case IDTushort:
			alphaPix.v.s[0] = alphaPix.v.s[1] = alphaPix.v.s[2] = 0;
			alphaPix.v.s[3] = 0xffff;
			break;
		case IDTfloat:
			alphaPix.v.f[0] = alphaPix.v.f[1] = alphaPix.v.f[2] = .0f;
			alphaPix.v.f[3] = 1.f;
			break;
		default:
			DMESG(DMCparameter, "Unhandled data type in PmarryAlpha()");
			return false;
		}
		if (!PsetImage(ib, alphaPix))
			return false;
	}
	dsiz = ImgDataSize[ia->csp->dtype];
	y = ia->yres;
	while (y--) {				/* copy RGB channels */
		uby8		*dp = ProwPtr(ib,y);
		const uby8	*sp = ProwPtr(ia,y);
		int		n = ia->xres;
		while (n--) {
			memcpy(dp, sp, 3*dsiz);
			dp += 3*dsiz;
			sp += 4*dsiz;
		}
	}
	PfreeImage(&tempImg);
	return true;
}

/* Copy the indicated source component to destination image */
int
PcopyComponent(ImgStruct *ib, const int cb, const ImgStruct *ia, const int ca)
{
	int	nca, ncb;
	int	y;

	if ((ib == NULL) | (ia == NULL) || (ib->csp == NULL) | (ia->csp == NULL))
		return false;
	if (ib->img == NULL) {
		ib->xres = ia->xres;
		ib->yres = ia->yres;
	} else if ((ib->xres != ia->xres) | (ib->yres != ia->yres)) {
		DMESG(DMCparameter, "Mismatched dimensions in PcopyComponent()");
		return false;
	}
	nca = ImgPixelLen[ia->csp->format];
	ncb = ImgPixelLen[ib->csp->format];
	if ((ca < 0) | (ca >= nca) | (cb < 0) | (cb >= ncb)) {
		DMESG(DMCparameter, "Illegal channel in PcopyComponent()");
		return false;
	}
	if (ca == cb) {
		if (ia->img == ib->img)
			return true;
		if (PimagesOverlap(ia, ib)) {
			DMESG(DMCparameter, "Copy on same channel to self");
			return false;
		}
	}
	if (!PmatchColorSpace(ia->csp, ib->csp, PICMdtype|PICMgamma)) {
		DMESG(DMCparameter, "Mismatched data types in PcopyComponent()");
		return false;
	}
	if (!PnewImage(ib, .0))			/* make sure dest. ready */
		return false;
	for (y = 0; y < ia->yres; y++) {
		const int	dsz = ImgDataSize[ia->csp->dtype];
		const int	ssrc = (nca-1)*dsz;
		const int	sdst = (ncb-1)*dsz;
		const uby8	*psrc = ProwPtr(ia,y) + ca*dsz;
		uby8		*pdst = ProwPtr(ib,y) + cb*dsz;
		int		nx = ia->xres;
		while (nx--) {
			int	nb = dsz;
			while (nb--)
				*pdst++ = *psrc++;
			psrc += ssrc;
			pdst += sdst;
		}
	}
	return true;
}

/* Prepare write buffer based on image structure */
int
PsetWriteBuf(ImgWriteBuf *iw, const ImgStruct *isrc)
{
	if (isrc == NULL || (isrc->csp == NULL) | (isrc->img == NULL))
		return false;
	if (iw == NULL)
		return false;
	iw->xres = isrc->xres;
	iw->yres = isrc->yres;
	iw->rowsize = isrc->rowsize;
	iw->csp = isrc->csp;
	iw->pixAspect = 1.f;
	iw->info = defImgInfo;
	iw->img = isrc->img;
	return true;
}
