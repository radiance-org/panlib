/*
 *  readrad.c
 *  panlib
 *
 *  Radiance picture file reader for Pancine.
 *
 *  Created by gward on Tue May 22 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#include <stdlib.h>
#include <math.h>
#include "system.h"
#include "imgreader.h"
#include "radheader.h"
#include "resolu.h"
#include "color.h"
#include "tonemap.h"

#ifndef RR_MAXHISTO
#define RR_MAXHISTO	(512L*512L)	/* maximum pixels to histogram */
#endif

typedef struct {
	BASE_IMGREADER;			/* base struct members (first!) */
	FILE		*fp;		/* picture file pointer */
	long		*frpos;		/* frame positions */
	long		frsiz;		/* size of frpos array */
	RESOLU		rs;		/* Radiance resolution spec. */
	RGBPRIMP	pr;		/* RGB primaries (NULL if XYZE) */
	int		ncs;		/* number of color samples (usu. 3) */
	float		wpt[4];		/* wavelength partitions */
	COLR		*sl;		/* scanline input buffer */
	long		*sp;		/* scanline seek positions */
	int		sn;		/* next scanline number */
	double		stonits;	/* conversion to cd/m^2 */
	double		LdDyn;		/* dynamic range of display */
	double		LdMax;		/* maximum display luminance */
	int		humanVis;	/* match human visibility? */
	TMstruct	*ts;		/* tone-mapping structure */
	TMbright	*lbuf;		/* luminance buffer for TM */
	RGBPRIMS	myprims;	/* custom RGB input primaries */
	RGBPRIMS	monprims;	/* display monitor primaries */
} RadianceReader;

#define RRcopyPrims(p1,p2)	memcpy(p1, p2, sizeof(RGBPRIMS))

/* Free allocated read buffers */
static void
RRfreeBuffers(RadianceReader *rr)
{
	if (rr->sl != NULL) {
		free(rr->sl);
		rr->sl = NULL;
	}
	if (rr->sp != NULL) {
		free(rr->sp);
		rr->sp = NULL;
	}
	if (rr->ts != NULL) {
		tmDone(rr->ts);
		rr->ts = NULL;
	}
	if (rr->lbuf != NULL) {
		free(rr->lbuf);
		rr->lbuf = NULL;
	}
}

/* Allocate new read buffers, returning 0 on success or -1 on failure */
static int
RRallocBuffers(RadianceReader *rr)
{
	RRfreeBuffers(rr);
	rr->sp = (long *)calloc(numscans(&rr->rs)+1, sizeof(long));
	rr->sl = (COLR *)malloc(scanlen(&rr->rs)*sizeof(COLR));
	if ((rr->sp == NULL) | (rr->sl == NULL)) {
		RRfreeBuffers(rr);
		return -1;
	}
	return 0;
}

/* Close open file and free allocated reader struct members */
static void
RRcloseFile(RadianceReader *rr)
{
	if (rr->fp != NULL) {
		fclose(rr->fp);
		rr->fp = NULL;
	}
	if (rr->frpos != &rr->frsiz) {
		free(rr->frpos);
		rr->frpos = &rr->frsiz;
	}
	RRfreeBuffers(rr);
}

/* Compare two sets of primaries and return 0 if they are close enough */
static int
RRcomparePrims(RGBPRIMP rp1, RGBPRIMP rp2)
{
	int     n = 8;
	float	*p1 = (float *)rp1, *p2 = (float *)rp2;
	if (p1 == p2)
		return 0;
	if ((p1 == NULL) | (p2 == NULL))
		return 1;
	while (n--) {
		if (*p1 - *p2 + 1e-4 < 0)
			return 1;
		if (*p2 - *p1 + 1e-4 < 0)
			return 1;
		p1++; p2++;
	}
	return 0;
}

/* Check header line, modifying RadianceReader struct appropriately */
static int
RRcheckHeadline(char *hl, void *p)
{
	RadianceReader *	rr = (RadianceReader *)p;
	char			fmt[MAXFMTLEN];
	if (isheadid(hl))
		return 0;		/* accept any header ID */
	if (formatval(fmt, hl)) {	/* check pixel format */
		if (!strcmp(fmt, COLRFMT)) {
			rr->encoding = "RLE RGBE";
			return 0;
		}
		if (!strcmp(fmt, CIEFMT)) {
			rr->pr = TM_XYZPRIM;
			rr->stonits /= WHTEFFICACY;
			rr->encoding = "RLE XYZE";
			return 0;
		}
		if (!strcmp(fmt, SPECFMT)) {
			rr->encoding = "Hyperspectral Radiance";
			return 0;
		}
		rr->errCode = IREformat;
		sprintf(rr->errMsg, "unrecognized data format: %s", fmt);
		return -1;
	}
	if (isaspect(hl)) {		/* Radiance pixel aspect (invert) */
		double	av = aspectval(hl);
		if ((av < .1) | (av > 10.))
			return 0;	/* should report error? */
		rr->pixAspect /= av;
		return 0;
	}
	if (isprims(hl)) {		/* explicit RGB primaries */
		primsval(rr->myprims, hl);
		if (RRcomparePrims(rr->myprims, stdprims))
			rr->pr = rr->myprims;
		else
			rr->pr = stdprims;
		return 0;
	}
	if (isncomp(hl)) {		/* # components */
		rr->ncs = ncompval(hl);
		if (rr->ncs < 3) {
			sprintf(rr->errMsg, "illegal number of components: %d", rr->ncs);
			rr->errCode = IREformat;
			return -1;
		}
		return 0;
	}
	if (iswlsplit(hl)) {		/* wavelength partition */
		wlsplitval(rr->wpt, hl);
		return 0;
	}
	if (isexpos(hl)) {		/* modify stonits */
		double	ev = exposval(hl);
		if ((ev <= 1e-10) | (ev >= 1e10))
			return 0;	/* should report error? */
		rr->stonits /= ev;
		return 0;
	}
	if (!strncmp(hl,"FRAME=",6)) {	/* may be animated sequence */
		if (rr->nframes == 1)
			rr->nframes = 0;
		rr->frameType = IRFanim;
		rr->frameRate = 0.1;
		return 0;
	}
	return 0;			/* ignore the rest */
}

static int	RRnextScanline(RadianceReader *rr);

/* Set indicated scanline, returning 0 on success and -1 on failure */
static int
RRsetScanline(RadianceReader *rr, int n)
{
	if ((n < 0) | (n > numscans(&rr->rs)))
		return -1;
					/* record current scanline position */
	if ((rr->sn >= 0) & (rr->sn <= numscans(&rr->rs)) &&
			rr->sp[rr->sn] <= 0L)
		rr->sp[rr->sn] = ftell(rr->fp);
					/* is file pointer positioned? */
	if (n != rr->sn) {
		if (rr->sp[n] > 0L) {	/* seek to known position */
			if (fseek(rr->fp, rr->sp[n], 0) < 0) {
				strcpy(rr->errMsg, "seek error");
				rr->errCode = IREread;
				return rr->sn = -1;
			}
			rr->sn = n;
		} else {		/* scan to unknown position */
			if (n > rr->sn)
				if (rr->sp[rr->sn = 0] <= 0L)
					return rr->sn = -1;
			while (rr->sp[rr->sn+1] > 0L)
				rr->sn++;
			if (fseek(rr->fp, rr->sp[rr->sn], 0) < 0) {
				strcpy(rr->errMsg, "seek error");
				rr->errCode = IREread;
				return rr->sn = -1;
			}
			while (rr->sn < n)
				if (RRnextScanline(rr) < 0)
					return rr->sn = -1;
		}
	}
					/* set next scanline rectangle */
	switch (rr->rs.rt) {
	case YMAJOR|YDECR:
	case YMAJOR|YDECR|XDECR:
		rr->nr.xleft = 0; rr->nr.xright = rr->xres;
		rr->nr.ytop = rr->sn;
		rr->nr.ybottom = rr->sn + 1;
		break;
	case YMAJOR:
	case YMAJOR|XDECR:
		rr->nr.xleft = 0; rr->nr.xright = rr->xres;
		rr->nr.ytop = rr->yres - 1 - rr->sn;
		rr->nr.ybottom = rr->yres - rr->sn;
		break;
	case YDECR:
	case 0:
		rr->nr.xleft = rr->sn;
		rr->nr.xright = rr->sn + 1;
		rr->nr.ytop = 0; rr->nr.ybottom = rr->yres;
		break;
	case XDECR:
	case XDECR|YDECR:
		rr->nr.xleft = rr->xres - 1 - rr->sn;
		rr->nr.xright = rr->xres - rr->sn;
		rr->nr.ytop = 0; rr->nr.ybottom = rr->yres;
		break;
	default:
		return -1;
	}
	return 0;
}

/* Read scanline from current position and advance to next */
static int
RRnextScanline(RadianceReader *rr)
{
	if ((rr->sn < 0 ) | (rr->sn >= numscans(&rr->rs))) {
		strcpy(rr->errMsg, "attempt to read past end of file");
		rr->errCode = IREread;
		return -1;
	}
	if (fread2colrs(rr->sl, scanlen(&rr->rs), rr->fp, rr->ncs, rr->wpt) < 0) {
		sprintf(rr->errMsg, "unexpected EOF (scanline %d)", rr->sn);
		rr->errCode = IREtruncated;
		return -1;
	}
	RRsetScanline(rr, ++rr->sn);
					/* check for ending frame */
	if (!rr->nframes && rr->sn == numscans(&rr->rs)) {
		int	c = getc(rr->fp);
		if (c == EOF)
			rr->nframes = rr->frame + 1;
		else
			ungetc(c, rr->fp);
	}
	return 0;
}

/* Load header and begin frame, returning 0 on success or -1 on failure */
static int
RRinitFrame(RadianceReader *rr, int fn)
{
	int	i;

	if ((rr->fp == NULL) | (fn < 0)) {
		rr->errCode = IREunknown;
		return -1;
	}
	if (fn == rr->frame) {		/* check for same frame */
		rr->nr = rr->fr;
		rr->errCode = IREnone;
		return 0;
	}
					/* check for standard image */
	if ((rr->frame == 0) & (rr->nframes == 1)) {
		rr->errCode = IREtruncated;
		strcpy(rr->errMsg, "single-frame picture");
		return -1;
	}
					/* make frpos array big enough */
	if (fn > 0 && (rr->frpos == &rr->frsiz || fn >= rr->frsiz)) {
		if (rr->frpos == &rr->frsiz) {
			rr->frpos = (long *)calloc(60, sizeof(long));
			if (rr->frpos == NULL)
				goto memerr;
			rr->frpos[0] = rr->frsiz;
			rr->frsiz = 60;
		} else {
			rr->frpos = (long *)realloc(rr->frpos,
					sizeof(long)*2*rr->frsiz);
			if (rr->frpos == NULL)
				goto memerr;
			i = rr->frsiz;
			rr->frsiz *= 2;
			while (i < rr->frsiz)
				rr->frpos[i++] = 0L;
		}
	}
	if (fn == 0) {			/* set first frame initial values */
		if (fseek(rr->fp, rr->frpos[0], 0) < 0) {
			rr->errCode = IREread;
			strcpy(rr->errMsg, "seek error");
			return -1;
		}
		rr->encoding = NULL;
		rr->pr = stdprims;
		rr->pixAspect = 1;
		rr->stonits = WHTEFFICACY;
		rr->ncs = 3;
		memcpy(rr->wpt, WLPART, sizeof(rr->wpt));
	} else if (fn == rr->frame+1) {	/* else seek to desired frame */
		if (RRsetScanline(rr, numscans(&rr->rs)) < 0)
			return -1;
		rr->frpos[fn] = rr->sp[numscans(&rr->rs)];
	} else if (rr->frpos[fn] <= 0L) {
		for (i = fn-1; i > 0 && rr->frpos[i] <= 0L; i--)
			;
		if (RRinitFrame(rr, i) < 0)
			return -1;
		while (rr->frame < fn)
			if (RRinitFrame(rr, rr->frame+1) < 0)
				return -1;
		return 0;
	} else if (fseek(rr->fp, rr->frpos[fn], 0) < 0) {
		rr->errCode = IREread;
		strcpy(rr->errMsg, "seek error");
		return -1;
	}
					/* attempt to load header */
	if (getheader(rr->fp, &RRcheckHeadline, rr) < 0) {
		if (rr->errCode == IREnone) {
			rr->nframes = fn;
			rr->errCode = IREtruncated;
			if (fn)
				strcpy(rr->errMsg, "no more frames");
			else
				strcpy(rr->errMsg, "end of file in header");
		}
		return -1;
	}
					/* header was OK */
	rr->frame = fn;
					/* set preferred color space */
	if (rr->pr == rr->myprims) {
		i = colorprimsOK(rr->myprims);
		if (!i) {
			rr->errCode = IREunknown;
			strcpy(rr->errMsg, "illegal color primitives");
			return -1;
		}
		if (i < 0)		/* flag for XYZ color space */
			rr->pr = TM_XYZPRIM;
	}
	if (rr->pr != TM_XYZPRIM) {
		PcopyCS(&rr->cs, &ICS_RGB709);
		RRcopyPrims((RGBPRIMP)rr->cs.chroma, rr->pr);
		if (!rr->encoding)
			rr->encoding = "RLE RGBE";
	} else {
		PcopyCS(&rr->cs, &ICS_XYZ);
		if (!rr->encoding)
			rr->encoding = "RLE XYZE";
	}
					/* set up buffers */
	if (!fgetsresolu(&rr->rs, rr->fp)) {
		rr->errCode = IREformat;
		strcpy(rr->errMsg, "missing image resolution");
		return -1;
	}
	if (RRallocBuffers(rr) < 0)
		goto memerr;
	rr->sp[0] = ftell(rr->fp);
					/* initialize first scanline */
	rr->xres = rr->rs.xr;
	rr->yres = rr->rs.yr;
	if (RRsetScanline(rr, rr->sn = 0) < 0) {
		rr->errCode = IREformat;
		strcpy(rr->errMsg, "cannot initialize first scanline");
		return -1;
	}
	rr->fr = rr->nr;
					/* all is well */
	rr->errCode = IREnone;
	return 0;
memerr:
	rr->errCode = IREmemory;
	strcpy(rr->errMsg, "cannot allocate input buffers");
	return -1;
}

/* Open a Radiance picture file for reading */
static ImgReader *
RRopen(const char *fname)
{
	extern const ImgReaderInterface	IRInterfaceRad;
	RadianceReader			*rr;
	if (!fname || !*fname)		/* see if we can open it, first */
		return NULL;
	rr = (RadianceReader *)calloc(1, sizeof(RadianceReader));
	if (rr == NULL)
		return NULL;
	rr->fp = fopen(fname, "rb");
	if (rr->fp == NULL) {
		free(rr);
		return NULL;
	}
	rr->ri = &IRInterfaceRad;
	strlcpy(rr->file, fname, sizeof(rr->file));
					/* default tone-mapping parameters */
	rr->LdDyn = 100.;
	rr->LdMax = 100.;
	rr->humanVis = 0;
					/* buffers are unallocated */
	rr->frpos = &rr->frsiz;
					/* start with still image assumption */
	rr->frpos[0] = 0L;
	rr->frame = -1; rr->nframes = 1;
	rr->frameType = IRFnone; rr->frameRate = 0;
					/* initialize frame */
	if (RRinitFrame(rr, 0) < 0) {
		RRcloseFile(rr);
		if ((rr->errCode == IREtruncated) | (rr->errCode == IREformat)) {
			free(rr);
			return NULL;
		}
	}
					/* return reader */
	return (ImgReader *)rr;
}

/* Advance the specified number of frames in a sequence */
static ImgReadErr
RRseekFrame(ImgReader *ir, int offs, ImgSeekMode sm)
{
	RadianceReader	*rr = (RadianceReader *)ir;
	int		fn = 0;
	if (rr == NULL)
		return IREunknown;
	switch (sm) {
	case IRSabs:
		fn = offs; break;
	case IRSrel:
	case IRSadv:
		fn = rr->frame+offs; break;
	}
	RRinitFrame(rr, fn);
	return rr->errCode;
}

/* Get header information and bundle for calling application */
static ImgReadErr
RRgetInfo(ImgReader *ir, ImgInfo *info)
{
	RadianceReader	*rr = (RadianceReader *)ir;
	if (rr == NULL || rr->fp == NULL || info == NULL)
		return IREunknown;
	*info = defImgInfo;
	rr->sn = -1;			/* rewind to header */
	if (fseek(rr->fp, rr->frpos[0], 0) < 0) {
		strcpy(rr->errMsg, "cannot rewind to picture header");
		return rr->errCode = IREread;
	}
	info->flags |= IIFstonits;	/* a priori information */
	info->stonits = rr->stonits;
					/* scan main info. header */
	if (RHgetInfo(rr->fp, info) < 0) {
		strcpy(rr->errMsg, "error encountered reading info header");
		return rr->errCode = IREunknown;
	}
					/* scan header for this frame */
	if (rr->frame > 0 && rr->frpos[rr->frame] > 0L &&
			(fseek(rr->fp, rr->frpos[rr->frame], 0) < 0 ||
			RHgetInfo(rr->fp, info) < 0)) {
		strcpy(rr->errMsg, "error encountered reading frame header");
		return rr->errCode = IREunknown;
	}
	return rr->errCode = IREnone;	/* all done */
}

/* Tone-mapping advisory call */
static void
RRtoneMapping(ImgReader *ir, double LdDyn, double LdMax, int humanVis)
{
	RadianceReader	*rr = (RadianceReader *)ir;
	if (rr == NULL)
		return;
	rr->LdDyn = LdDyn;		/* just store parameters */
	rr->LdMax = LdMax;
	rr->humanVis = humanVis;
	if (rr->ts != NULL) {		/* reset tone-mapping */
		tmDone(rr->ts);
		rr->ts = NULL;
	}
}

/* Try to set up tone-mapping for this picture (if not set already) */
static int
RRtoneMapOK(RadianceReader *rr, ImgColorSpace *cs)
{
	int	ns = numscans(&rr->rs);
	int	slen = scanlen(&rr->rs);
	int	i;
					/* check if we tried already */
	if (rr->ts != NULL && !RRcomparePrims(rr->monprims, (RGBPRIMP)cs->chroma))
		return (rr->lbuf != NULL);
	if (cs->logorig > 0)		/* we don't handle log output */
		return 0;
	RRcopyPrims(rr->monprims, (RGBPRIMP)cs->chroma);
	if (rr->ts != NULL)		/* start over */
		tmDone(rr->ts);
	rr->ts = tmInit((rr->humanVis?TM_F_HUMAN:TM_F_CAMERA)|TM_F_NOSTDERR,
			rr->monprims, cs->gamma);
	if (rr->ts == NULL)
		return 0;
	if (!RRcomparePrims(rr->monprims, rr->pr))
		rr->pr = rr->monprims;	/* input primaries same as output */
	if (tmSetSpace(rr->ts, rr->pr, rr->stonits) != TM_E_OK)
		goto tmerror;
	if (rr->lbuf == NULL) {		/* allocate luminance buffer */
		rr->lbuf = (TMbright *)malloc(slen*sizeof(TMbright));
		if (rr->lbuf == NULL)
			goto tmerror;
	}
					/* compute histogram */
	if ((long)rr->xres*rr->yres > RR_MAXHISTO) {
		int		j, step;
		step = (int)(sqrt((double)rr->xres*rr->yres/RR_MAXHISTO) + 1.);
		for (i = 0; i < ns; i += step) {
			if (RRsetScanline(rr, i) < 0)
				goto tmerror;
			if (RRnextScanline(rr) < 0)
				goto tmerror;
			for (j = slen/step; j--; )
				if (tmCvColrs(rr->ts, rr->lbuf+j, TM_NOCHROM,
						rr->sl+j*step, 1) != TM_E_OK)
					goto tmerror;
			if (tmAddHisto(rr->ts, rr->lbuf, slen/step, 1) != TM_E_OK)
				goto tmerror;
		}
	} else {
		if (RRsetScanline(rr, 0) < 0)
			goto tmerror;
		for (i = 0; i < ns; i++) {
			if (RRnextScanline(rr) < 0)
				goto tmerror;
			if (tmCvColrs(rr->ts, rr->lbuf, TM_NOCHROM, rr->sl, slen)
					!= TM_E_OK)
				goto tmerror;
			if (tmAddHisto(rr->ts, rr->lbuf, slen, 1) != TM_E_OK)
				goto tmerror;
		}
	}
					/* compute tone-mapping */
	if (tmComputeMapping(rr->ts, cs->gamma, rr->LdDyn, rr->LdMax) != TM_E_OK)
		goto tmerror;
	return 1;			/* all is well */
tmerror:
	if (rr->lbuf != NULL) {		/* deallocate conversion buffers */
		free(rr->lbuf);
		rr->lbuf = NULL;
	}
	return 0; 
}

/* Set up scanline limits for the given image rectangle */
static int
RRscanLim(RadianceReader *rr, int sl[2], ImgRect *rl)
{
	if (!PlegalRect(rl, rr->xres, rr->yres))
		return -1;
	switch (rr->rs.rt) {
	case YMAJOR|YDECR:
	case YMAJOR|YDECR|XDECR:
		sl[0] = rl->ytop;
		sl[1] = rl->ybottom;
		break;
	case YMAJOR:
	case YMAJOR|XDECR:
		sl[0] = rr->yres - rl->ybottom;
		sl[1] = rr->yres - rl->ytop;
		break;
	case 0:
	case YDECR:
		sl[0] = rl->xleft;
		sl[1] = rl->xright;
		break;
	case XDECR:
	case XDECR|YDECR:
		sl[0] = rr->xres - rl->xright;
		sl[1] = rr->xres - rl->xleft;
		break;
	default:
		return -1;
	}
					/* position at first scanline */
	if (rr->sn == sl[0])
		return 0;
	return RRsetScanline(rr, sl[0]);
}

/* Read scanline to display buffer, applying conversions and subsampling */
static int
RRscanConvert(RadianceReader *rr, ImgReadBuf *rb, int si)
{
	int	psiz = ImgPixelSize(&rb->cs);
	int	i;
	int	ec;
	int	sbeg;
	uby8	*bpos;
	int	bstep;
	ssize_t	npix;
					/* position input */
	if (rr->sn != si && RRsetScanline(rr, si) < 0)
		return -1;
	switch (rr->rs.rt) {		/* correct for orientation */
	case YMAJOR|YDECR:
	case YMAJOR:
		sbeg = rb->r.xleft;
		npix = (rb->r.xright - rb->r.xleft)/rb->subsample;
		bpos = rb->buf +
			(rr->nr.ytop - rb->r.ytop)/rb->subsample*npix*psiz;
		bstep = psiz;
		break;
	case YMAJOR|YDECR|XDECR:
	case YMAJOR|XDECR:
		sbeg = rr->xres - rb->r.xright;
		npix = (rb->r.xright - rb->r.xleft)/rb->subsample;
		bpos = rb->buf +
			(rr->nr.ytop - rb->r.ytop)/rb->subsample*npix*psiz
				+ (npix - 1)*psiz;
		bstep = -psiz;
		break;
	case XDECR:
	case 0:
		sbeg = rr->yres - rb->r.ybottom;
		npix = (rb->r.ybottom - rb->r.ytop)/rb->subsample;
		bpos = rb->buf + ImgReadBufLen(rb)
			- (rb->r.xright - rr->nr.xleft)/rb->subsample*psiz;
		bstep = -(rb->r.xright - rb->r.xleft)/rb->subsample*psiz;
		break;
	case YDECR:
	case XDECR|YDECR:
		sbeg = rb->r.ytop;
		npix = (rb->r.ybottom - rb->r.ytop)/rb->subsample;
		bpos = rb->buf +
			(rr->nr.xleft - rb->r.xleft)/rb->subsample*psiz;
		bstep = (rb->r.xright - rb->r.xleft)/rb->subsample*psiz;
		break;
	default:
		return -1;
	}
	if (RRnextScanline(rr) < 0)	/* read scanline */
		return -1;
					/* covert to output color space */
	if (rb->cs.dtype == IDTubyte && rb->cs.format == IPFrgb) {
		if (rr->lbuf == NULL)
			return -1;	/* should never happen! */
		if ((bstep == 3) & (rb->subsample == 1)) {
			ec = tmCvColrs(rr->ts, rr->lbuf, bpos, rr->sl+sbeg, npix);
			if (ec != TM_E_OK)
				goto tmerr;
			ec = tmMapPixels(rr->ts, bpos, rr->lbuf, bpos, npix);
			if (ec != TM_E_OK)
				goto tmerr;
			return 0;	/* transferred directly */
		}
					/* else map pixel by pixel */
		for (i = sbeg; npix-- > 0; bpos += bstep, i += rb->subsample) {
			ec = tmCvColrs(rr->ts, rr->lbuf, bpos, rr->sl+i, 1);
			if (ec != TM_E_OK)
				goto tmerr;
			ec = tmMapPixels(rr->ts, bpos, rr->lbuf, bpos, 1);
			if (ec != TM_E_OK)
				goto tmerr;
		}
		return 0;		/* all done */
	}
	if ((rb->cs.dtype != IDTfloat) | (rb->cs.format != rr->cs.format)) {
		strcpy(rr->errMsg, "pixel format mismatch in RRscanConvert");
		rr->errCode = IREunknown;
		return -1;
	}
					/* convert COLR to float[3] */
	for (i = sbeg; npix-- > 0; bpos += bstep, i += rb->subsample)
		colr_color((float *)bpos, rr->sl[i]);
	return 0;			/* finito */
tmerr:					/* tone-mapping error jumps here */
	strcpy(rr->errMsg, tmErrorMessage[ec]);
	rr->errCode = IREunknown;
	return -1;
}

/* Read a rectangle, in the specified color space if possible */
static ImgReadErr
RRreadRec(ImgReader *ir, ImgReadBuf *rb)
{
	RadianceReader	*rr = (RadianceReader *)ir;
	int		slim[2];
	int		i;
	if (rr == NULL || rr->fp == NULL || rb == NULL)
		return IREunknown;
					/* can we convert to requested CS? */
	if (rb->cs.dtype == IDTubyte && rb->cs.format == IPFrgb) {
		if (!RRtoneMapOK(rr, &rb->cs)) {
			rb->cs = rr->cs;
			rb->buf = NULL;
		}
	} else if (!PmatchColorSpace(&rb->cs, &rr->cs, PICMall)) {
		rb->cs = rr->cs;	/* otherwise, make caller work */
		rb->buf = NULL;
	}
	ImgFixSampling(rb);		/* align rectangle request */
	rr->errCode = IREnone;		/* set scanline limits for rect */
	if (RRscanLim(rr, slim, &rb->r) < 0) {
		if (rr->errCode == IREnone) {
			strcpy(rr->errMsg, "illegal read rectangle");
			rr->errCode = IREunknown;
		}
		return rr->errCode;
	}
	if (!rb->buf) {			/* need new buffer */
		int	buflen;
		rb->r = rr->nr;		/* use scanline-sized buffer */
		ImgFixSampling(rb);
		slim[1] = slim[0] + rb->subsample;
		buflen = ImgReadBufLen(rb);
		if (buflen <= 0) {
			strcpy(rr->errMsg, "internal buffer size error");
			return rr->errCode = IREunknown;
		}
		rb->buf = (uby8 *)malloc(buflen);
		if (rb->buf == NULL) {
			strcpy(rr->errMsg, "cannot allocate new read buffer");
			return rr->errCode = IREmemory;
		}
	}
					/* transfer each scanline */
	for (i = slim[0]; i < slim[1]; i += rb->subsample)
		if (RRscanConvert(rr, rb, i) < 0)
			break;
	return rr->errCode;
}

/* Close a Radiance picture and free reader struct */
static void
RRclose(ImgReader *ir)
{
	RadianceReader	*rr = (RadianceReader *)ir;
	if (rr == NULL)
		return;
	RRcloseFile(rr);
	free(rr);
}

/* Interface for Radiance picture reader */
const ImgReaderInterface IRInterfaceRad = {
	"Radiance.hdr.pic.rgbe.xyze.hsr",
	&RRopen, &RRseekFrame, &RRgetInfo, NULL,
	&RRtoneMapping, &RRreadRec, &RRclose
};
