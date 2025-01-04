/*
 *  readdpt.c
 *  panlib
 *
 *  Encoded 16-bit depth map reader for Pancine.
 *
 *  Created by gward on July 26 2019.
 *  Copyright (c) 2019 Anyhere Software. All rights reserved.
 *
 */

#include <stdlib.h>
#include <math.h>
#include "imgreader.h"
#include "radheader.h"
#include "fvect.h"
#include "depthcodec.h"

typedef struct {
	BASE_IMGREADER;			/* base struct members (first!) */
	DEPTHCODEC	dc;		/* depth codec */
} DepthReader;

/* Open a depth map and read header */
static ImgReader *
DMopen(const char *fname)
{
	extern const ImgReaderInterface	IRInterfaceDPT;
	DepthReader			*dr;
	if (!fname || !*fname)
		return NULL;		/* see if we can open it, first */
	dr = (DepthReader *)calloc(1, sizeof(DepthReader));
	if (!dr)
		return NULL;
	set_dc_defaults(&dr->dc);
	dr->dc.finp = fopen(fname, "rb");
	if (!dr->dc.finp) {
		free(dr);
		return NULL;
	}				/* read header & get resolution */
	strlcpy(dr->file, fname, sizeof(dr->file));
	dr->dc.inpname = dr->file;
	dr->dc.hdrflags = HF_HEADIN|HF_RESIN;
	if (!process_dc_header(&dr->dc, 0, NULL) ||
			!check_decode_depths(&dr->dc)) {
		fclose(dr->dc.finp);
		free(dr);
		return NULL;
	}
	dr->ri = &IRInterfaceDPT;
	PcopyCS(&dr->cs, &ICS_D);
	dr->encoding = "16-bit encoded depth";
	dr->xres = scanlen(&dr->dc.res);
	dr->yres = numscans(&dr->dc.res);
	if (dr->dc.gotview)
		dr->dc.gotview = (setview(&dr->dc.vw) == NULL);
	if (dr->dc.gotview)
		dr->pixAspect = viewaspect(&dr->dc.vw) * dr->xres / dr->yres;
	else
		dr->pixAspect = 1;
	dr->frame = 0; dr->nframes = 1;
	dr->frameType = IRFnone; dr->frameRate = 0;
	dr->fr.xleft = 0; dr->fr.xright = dr->xres;
	dr->fr.ytop = 0; dr->fr.ybottom = 1;
	dr->nr = dr->fr;		/* we're ready to roll */
	return (ImgReader *)dr;
}

/* Close depth map */
static void
DMcloseFile(DepthReader *dr)
{
	if (!dr || !dr->dc.finp)
		return;
	fclose(dr->dc.finp);
	dr->dc.finp = NULL;
}

/* Get header information and bundle for calling application */
static ImgReadErr
DMgetInfo(ImgReader *ir, ImgInfo *info)
{
	DepthReader	*dr = (DepthReader *)ir;
	if (!dr || !dr->dc.finp || !info)
		return IREunknown;
	*info = defImgInfo;		/* rewind to beginning */
	if (fseek(dr->dc.finp, 0, SEEK_SET) < 0) {
		strcpy(dr->errMsg, "cannot rewind to header");
		DMcloseFile(dr);
		return dr->errCode = IREread;
	}
					/* scan info. header */
	if (RHgetInfo(dr->dc.finp, info) < 0) {
		strcpy(dr->errMsg, "error encountered reading info header");
		DMcloseFile(dr);
		return dr->errCode = IREread;
	}
	if (dr->dc.gotview) {		/* copy view parameters? */
		strcpy(info->view, viewopt(&dr->dc.vw));
		info->flags |= IIFview;
		if (dr->dc.vw.type == VT_PER) {
			info->hvangle = dr->dc.vw.horiz;
			info->flags |= IIFhvangle;
		} else if (dr->dc.vw.type == VT_PAR) {
			info->hvangle = 0;
			info->flags |= IIFhvangle;
		}
	}
	if (dr->dc.res.rt != PIXSTANDARD) {
		int	i = 8;		/* relay map orientation */
		while (i--)
			if (dr->dc.res.rt == RHortab[i]) {
				info->orientation = i + 1;
				info->flags |= IIForientation;
				break;
			}
	}
					/* reset file position */
	if (fseek(dr->dc.finp, dr->dc.curpos, SEEK_SET) < 0) {
		strcpy(dr->errMsg, "seek error");
		DMcloseFile(dr);
		return dr->errCode = IREread;
	}
	return dr->errCode = IREnone;	/* all done */
}

/* Read and convert a scanline */
static int
DMscanConvert(DepthReader *dr, ImgReadBuf *rb, const int y)
{
	int	x;
	uby8	*bufp;
					/* align to start */
	x = seek_dc_pix(&dr->dc, rb->r.xleft, y);
	if (x < 0) {
		strcpy(dr->errMsg, "seek error for requested scanline");
		DMcloseFile(dr);
		dr->errCode = IREread;
		return -1;
	}
	if (x == 0) {
		strcpy(dr->errMsg, "bad read rectangle");
		dr->errCode = IREunknown;
		return -1;
	}
	bufp = rb->buf + (size_t)ImgPixelSize(&rb->cs) *
			((y - rb->r.ytop)/rb->subsample) *
			((rb->r.xright - rb->r.xleft)/rb->subsample);
					/* read each depth pixel */
	for (x = rb->r.xleft; x < rb->r.xright; x += rb->subsample) {
		long	pnext;
		int	c = getint(2, dr->dc.finp);
		if (c == EOF && feof(dr->dc.finp))
			goto hitEOF;
		dr->dc.curpos += 2;
		switch (rb->cs.dtype) {
		case IDTfloat:
			*(float *)bufp = code2depth(c, dr->dc.refdepth);
			bufp += sizeof(float);
			break;
		case IDTushort:
			*(unsigned short *)bufp = (1<<15) + c;
			bufp += sizeof(unsigned short);
			break;
		case IDTubyte:
			*bufp++ = ((1<<15) + c) >> 8;
			break;
		default:
			strcpy(dr->errMsg, "unsupported color space");
			dr->errCode = IREunsupported;
			return -1;
		}
		pnext = rb->subsample - 1;	/* skip to next sample/row */
		if (pnext > dr->nr.xright - x - 1)
			pnext = dr->nr.xright - x - 1;
		pnext = dr->dc.curpos + 2*pnext;
		while (dr->dc.curpos < pnext) {
			if (getc(dr->dc.finp) == EOF)
				goto hitEOF;
			dr->dc.curpos++;
		}
	}
	dr->nr.ytop = y+1;
	dr->nr.ybottom = y+2;
	return 0;
hitEOF:
	strcpy(dr->errMsg, "unexpected EOF");
	DMcloseFile(dr);
	dr->errCode = IREtruncated;
	return -1;
}

/* Read a rectangle, in the specified color space if possible */
static ImgReadErr
DMreadRec(ImgReader *ir, ImgReadBuf *rb)
{
	DepthReader	*dr = (DepthReader *)ir;
	int		i;
	if (!dr || !dr->dc.finp || !rb)
		return IREunknown;
					/* can we convert to requested CS? */
	if ((rb->cs.format != IPFd) & (rb->cs.format != IPFy)) {
		rb->cs = dr->cs;	/* make caller work */
		rb->buf = NULL;
	}
	dr->errCode = IREnone;
	if (!rb->buf) {			/* need new buffer */
		int	buflen;
		rb->r.xleft = 0;	/* use scanline-sized buffer */
		rb->r.xright = dr->xres;
		rb->r.ybottom = rb->r.ytop + rb->subsample;
		ImgFixSampling(rb);
		buflen = ImgReadBufLen(rb);
		if (buflen <= 0) {
			strcpy(dr->errMsg, "internal buffer size error");
			return dr->errCode = IREunknown;
		}
		rb->buf = (uby8 *)malloc(buflen);
		if (rb->buf == NULL) {
			strcpy(dr->errMsg, "cannot allocate new read buffer");
			return dr->errCode = IREmemory;
		}
	} else				/* else align rectangle request */
		ImgFixSampling(rb);
					/* transfer each scanline */
	for (i = rb->r.ytop; i < rb->r.ybottom; i += rb->subsample)
		if (DMscanConvert(dr, rb, i) < 0)
			break;

	return dr->errCode;
}

/* Close a depth map and free reader struct */
static void
DMclose(ImgReader *ir)
{
	DepthReader	*dr = (DepthReader *)ir;
	if (!dr)
		return;
	DMcloseFile(dr);
	free(dr);
}

/* Interface for encoded depth map reader */
const ImgReaderInterface IRInterfaceDPT = {
	"Depth.dpt",
	&DMopen, NULL, &DMgetInfo, NULL,
	NULL, &DMreadRec, &DMclose
};
