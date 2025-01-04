/*
 *  readnrm.c
 *  panlib
 *
 *  Encoded 32-bit normal map reader for Pancine.
 *
 *  Created by gward on July 26 2019.
 *  Copyright (c) 2019 Anyhere Software. All rights reserved.
 *
 */

#include <stdlib.h>
#include "rtmath.h"
#include "imgreader.h"
#include "radheader.h"
#include "fvect.h"
#include "normcodec.h"

typedef struct {
	BASE_IMGREADER;			/* base struct members (first!) */
	NORMCODEC	nc;		/* normal codec */
} NormReader;

/* Open a normal map and read header */
static ImgReader *
NMopen(const char *fname)
{
	extern const ImgReaderInterface	IRInterfaceNRM;
	NormReader			*nr;
	if (!fname || !*fname)
		return NULL;		/* see if we can open it, first */
	nr = (NormReader *)calloc(1, sizeof(NormReader));
	if (!nr)
		return NULL;
	set_nc_defaults(&nr->nc);
	nr->nc.finp = fopen(fname, "rb");
	if (!nr->nc.finp) {
		free(nr);
		return NULL;
	}				/* read header & get resolution */
	strlcpy(nr->file, fname, sizeof(nr->file));
	nr->nc.inpname = nr->file;
	nr->nc.hdrflags = HF_HEADIN|HF_RESIN;
	if (!process_nc_header(&nr->nc, 0, NULL) ||
			!check_decode_normals(&nr->nc)) {
		fclose(nr->nc.finp);
		free(nr);
		return NULL;
	}
	nr->ri = &IRInterfaceNRM;
	PcopyCS(&nr->cs, &ICS_ENORM);
	nr->encoding = "32-bit encoded direction";
	nr->xres = scanlen(&nr->nc.res);
	nr->yres = numscans(&nr->nc.res);
	nr->pixAspect = 1;		/* XXX assumption */
	nr->frame = 0; nr->nframes = 1;
	nr->frameType = IRFnone; nr->frameRate = 0;
	nr->fr.xleft = 0; nr->fr.xright = nr->xres;
	nr->fr.ytop = 0; nr->fr.ybottom = 1;
	nr->nr = nr->fr;		/* we're ready to roll */
	return (ImgReader *)nr;
}

/* Close normal map */
static void
NMcloseFile(NormReader *nr)
{
	if (!nr || !nr->nc.finp)
		return;
	fclose(nr->nc.finp);
	nr->nc.finp = NULL;
}

/* Get header information and bundle for calling application */
static ImgReadErr
NMgetInfo(ImgReader *ir, ImgInfo *info)
{
	NormReader	*nr = (NormReader *)ir;
	long		curpos;
	if (!nr || !nr->nc.finp || !info)
		return IREunknown;
	*info = defImgInfo;
	curpos = ftell(nr->nc.finp);	/* rewind to beginning */
	if (fseek(nr->nc.finp, 0, SEEK_SET) < 0) {
		strcpy(nr->errMsg, "cannot rewind to header");
		NMcloseFile(nr);
		return nr->errCode = IREread;
	}
					/* scan info. header */
	if (RHgetInfo(nr->nc.finp, info) < 0) {
		strcpy(nr->errMsg, "error encountered reading info header");
		NMcloseFile(nr);
		return nr->errCode = IREread;
	}
	if (nr->nc.res.rt != PIXSTANDARD) {
		int	i = 8;		/* relay map orientation */
		while (i--)
			if (nr->nc.res.rt == RHortab[i]) {
				info->orientation = i + 1;
				info->flags |= IIForientation;
				break;
			}
	}
					/* reset file position */
	if (fseek(nr->nc.finp, curpos, SEEK_SET) < 0) {
		strcpy(nr->errMsg, "seek error");
		NMcloseFile(nr);
		return nr->errCode = IREread;
	}
	return nr->errCode = IREnone;	/* all done */
}

/* Read and convert a scanline */
static int
NMscanConvert(NormReader *nr, ImgReadBuf *rb, const int y)
{
	int	x;
	uby8	*bufp;
					/* align to start */
	x = seek_nc_pix(&nr->nc, rb->r.xleft, y);
	if (x < 0) {
		strcpy(nr->errMsg, "seek error for requested scanline");
		NMcloseFile(nr);
		nr->errCode = IREread;
		return -1;
	}
	if (x == 0) {
		strcpy(nr->errMsg, "bad read rectangle");
		nr->errCode = IREunknown;
		return -1;
	}
	bufp = rb->buf + (size_t)ImgPixelSize(&rb->cs) *
			((y - rb->r.ytop)/rb->subsample) *
			((rb->r.xright - rb->r.xleft)/rb->subsample);
					/* read each normal pixel */
	for (x = rb->r.xleft; x < rb->r.xright; x += rb->subsample) {
		long	toskip;
		int	c;
		FVECT	nrm;
		float	frgb[3];
		int32	nc = getint(4, nr->nc.finp);
		if (nc == EOF && feof(nr->nc.finp))
			goto hitEOF;
		switch (rb->cs.format) {
		case IPFn:
			*(int32 *)bufp = nc;
			bufp += sizeof(int32);
			break;
		case IPFvec3:
			decodedir(nrm, nc);
			VCOPY((float *)bufp, nrm);
			bufp += 3*sizeof(float);
			break;
		case IPFrgb:
			decodedir(nrm, nc);
			for (c = 3; c--; )
				frgb[c] = .5 + .5*nrm[c];
			switch (rb->cs.dtype) {
			case IDTfloat:
				memcpy(bufp, frgb, sizeof(frgb));
				bufp += sizeof(frgb);
				break;
			case IDTubyte:
				for (c = 0; c < 3; c++)
					*bufp++ = (int)(frgb[c]*255.9);
				break;
			case IDTushort:
				for (c = 0; c < 3; c++) {
					*(unsigned short *)bufp = (int)(frgb[c]*65535.9);
					bufp += sizeof(unsigned short);
				}
				break;
			default:
				goto badColorSpace;
			}
			break;
		default:;
badColorSpace:
			strcpy(nr->errMsg, "unsupported color space");
			nr->errCode = IREunsupported;
			return -1;
		}
		toskip = rb->subsample - 1;	/* skip to next sample/row */
		if (toskip > nr->nr.xright - x - 1)
			toskip = nr->nr.xright - x - 1;
		toskip *= 4;
		while (toskip-- > 0)
			if (getc(nr->nc.finp) == EOF)
				goto hitEOF;
	}
	nr->nr.ytop = y+1;
	nr->nr.ybottom = y+2;
	return 0;
hitEOF:
	strcpy(nr->errMsg, "unexpected EOF");
	NMcloseFile(nr);
	nr->errCode = IREtruncated;
	return -1;
}

/* Read a rectangle, in the specified color space if possible */
static ImgReadErr
NMreadRec(ImgReader *ir, ImgReadBuf *rb)
{
	NormReader	*nr = (NormReader *)ir;
	int		i;
	if (!nr || !nr->nc.finp || !rb)
		return IREunknown;

	switch (rb->cs.format) {	/* can we convert to requested CS? */
	case IPFn:
		if (rb->cs.dtype != IDTint)
			goto makeCallerConvert;
		break;
	case IPFvec3:
		if (rb->cs.dtype != IDTfloat)
			goto makeCallerConvert;
		break;
	case IPFrgb:
		if (rb->cs.dtype == IDTint)
			goto makeCallerConvert;
		break;
	default:;
makeCallerConvert:
		rb->cs = nr->cs;
		rb->buf = NULL;
		break;
	}
	nr->errCode = IREnone;
	if (!rb->buf) {			/* need new buffer */
		int	buflen;
		rb->r.xleft = 0;	/* use scanline-sized buffer */
		rb->r.xright = nr->xres;
		rb->r.ybottom = rb->r.ytop + rb->subsample;
		ImgFixSampling(rb);
		buflen = ImgReadBufLen(rb);
		if (buflen <= 0) {
			strcpy(nr->errMsg, "internal buffer size error");
			return nr->errCode = IREunknown;
		}
		rb->buf = (uby8 *)malloc(buflen);
		if (rb->buf == NULL) {
			strcpy(nr->errMsg, "cannot allocate new read buffer");
			return nr->errCode = IREmemory;
		}
	} else				/* else align rectangle request */
		ImgFixSampling(rb);
					/* transfer each scanline */
	for (i = rb->r.ytop; i < rb->r.ybottom; i += rb->subsample)
		if (NMscanConvert(nr, rb, i) < 0)
			break;	
	return nr->errCode;
}

/* Close a normal map and free reader struct */
static void
NMclose(ImgReader *ir)
{
	NormReader	*nr = (NormReader *)ir;
	if (!nr)
		return;
	NMcloseFile(nr);
	free(nr);
}

/* Interface for encoded normal map reader */
const ImgReaderInterface IRInterfaceNRM = {
	"Normal.nrm",
	&NMopen, NULL, &NMgetInfo, NULL,
	NULL, &NMreadRec, &NMclose
};
