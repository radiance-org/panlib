/*
 *  readmtx.c
 *  panlib
 *
 *  Read a Radiance matrix as image, usually float but can be ASCII or double.
 *
 *  Created by Greg Ward on 3/30/20.
 *  Copyright 2020 Anyhere Software. All rights reserved.
 *
 */

#include <math.h>
#include "imgreader.h"
#include "radheader.h"
#include "resolu.h"
#include "color.h"

/* Element component types (first is default) */
typedef enum { MCTascii, MCTfloat, MCTdouble } MtxCompType;

static const size_t	MCsize[] = {0, sizeof(float), sizeof(double)};

static const char	*luminance_format[] = {
				"ASCII luminance",
				"32-bit-float luminance",
				"64-bit-double luminance"
			};
static const char	*xyz_format[] = {
				"ASCII CIE-XYZ",
				"32-bit-float CIE-XYZ",
				"64-bit-double CIE-XYZ"
			};
static const char	*rgb_format[] = {
				"ASCII RGB",
				"32-bit-float RGB",
				"64-bit-double RGB"
			};

typedef struct {
	BASE_IMGREADER;			/* base struct members (first!) */
	FILE		*fin;		/* input stream */
	MtxCompType	dt;		/* component data type */
	short		ncomp;		/* number of components */
	short		swapped;	/* swapped byte order? */
	float		stonits;	/* sample-to-nits conversion factor */
	int		start;		/* start of actual data */
} MatrixReader;

/* Check header line, modifying MatrixReader struct appropriately */
static int
MRcheckHeadline(char *hl, void *p)
{
	MatrixReader *	mr = (MatrixReader *)p;
	char		fmt[MAXFMTLEN];
	int		i;

	if (isheadid(hl))
		return 0;		/* accept any header ID */
					/* matrix dimensions */
	if (!strncmp(hl, "NCOMP=", 6)) {
		mr->ncomp = atoi(hl+6);
		return 0;
	}
	if (!strncmp(hl, "NCOLS=", 6)) {
		mr->xres = atoi(hl+6);
		if (mr->xres <= 0) {
			mr->errCode = IREformat;
			strcpy(mr->errMsg, "bad NCOLS in header");
			return -1;
		}
	}
	if (!strncmp(hl, "NROWS=", 6)) {
		mr->yres = atoi(hl+6);
		if (mr->yres <= 0) {
			mr->errCode = IREformat;
			strcpy(mr->errMsg, "bad NROWS in header");
			return -1;
		}
	}
	if (formatval(fmt, hl)) {	/* check component type */
		if (!strcmp(fmt, "ascii")) {
			mr->dt = MCTascii;
			return 0;
		}
		if (!strcmp(fmt, "float")) {
			mr->dt = MCTfloat;
			return 0;
		}
		if (!strcmp(fmt, "double")) {
			mr->dt = MCTdouble;
			return 0;
		}
		mr->errCode = IREformat;
		sprintf(mr->errMsg, "unsupported data format: %s", fmt);
		return -1;
	}
	if ((i = isbigendian(hl)) >= 0) {
		mr->swapped = (i != nativebigendian());
		return 0;
	}
	if (isaspect(hl)) {		/* Radiance pixel aspect */
		float	av = aspectval(hl);
		if ((av < .1f) | (av > 10.f))
			return 0;	/* should report error? */
		mr->pixAspect /= av;
		return 0;
	}
	if (isprims(hl)) {		/* explicit color primaries */
		int	ok;
		primsval(mr->cs.chroma, hl);
		ok = colorprimsOK(mr->cs.chroma);
		if (!ok) {
			mr->errCode = IREformat;
			strcpy(mr->errMsg, "illegal color primitives");
			return -1;
		}
		if (ok < 0)		/* flag as XYZ */
			mr->cs.format = IPFxyz;
		return 0;
	}
	if (isexpos(hl)) {		/* modify stonits */
		float	ev = exposval(hl);
		if ((ev <= 1e-10f) | (ev >= 1e10f))
			return 0;	/* should report error? */
		mr->stonits /= ev;
		return 0;
	}
	return 0;			/* ignore the rest */
}

/* Close matrix file */
static void
MRcloseFile(MatrixReader *mr)
{
	if (!mr || !mr->fin)
		return;
	fclose(mr->fin);
	mr->fin = NULL;
}

/* Open a matrix image and read header */
static ImgReader *
MRopen(const char *fname)
{
	extern const ImgReaderInterface	IRInterfaceMTX;
	MatrixReader			*mr;

	if (!fname || !*fname)
		return NULL;		/* see if we can open it, first */
	mr = (MatrixReader *)calloc(1, sizeof(MatrixReader));
	if (!mr)
		return NULL;
	mr->fin = fopen(fname, "rb");
	if (!mr->fin) {
		free(mr);
		return NULL;
	}
	mr->ri = &IRInterfaceMTX;
	strlcpy(mr->file, fname, sizeof(mr->file));
	mr->pixAspect = 1.f;		/* set up defaults and scan header */
	mr->stonits = 1.f;
	mr->cs.dtype = IDTfloat;
	mr->cs.gamma = 1.f;
	memcpy(mr->cs.chroma, stdprims, sizeof(stdprims));
	mr->nframes = 1;		/* no animation support */
	if (getheader(mr->fin, MRcheckHeadline, mr) < 0) {
		if (mr->errCode == IREnone) {
			mr->errCode = IREtruncated;
			strcpy(mr->errMsg, "end of file in header");
		}
		MRcloseFile(mr);
		return (ImgReader *)mr;
	}
					/* file has resolution string? */
	if ((mr->xres <= 0) | (mr->yres <= 0) &&
			!fscnresolu(&mr->xres, &mr->yres, mr->fin)) {
		mr->errCode = IREformat;
		strcpy(mr->errMsg, "missing/bad matrix resolution");
		MRcloseFile(mr);
		return (ImgReader *)mr;
	}
	mr->start = ftell(mr->fin);	/* record data origin */
retry_cspace:
	switch (mr->ncomp) {		/* sort out interpretation */
	case 1:
		mr->cs.format = IPFy;
		mr->cs.chroma[0][1] = mr->cs.chroma[0][1] = 1.f/3.f;
		mr->encoding = luminance_format[mr->dt];
		break;
	case 3:
		if (PmatchColorSpace(&mr->cs, &ICS_XYZ, PICMchroma)) {
			mr->cs.format = IPFxyz;
			mr->encoding = xyz_format[mr->dt];
		} else {
			mr->cs.format = IPFrgb;
			mr->encoding = rgb_format[mr->dt];
			mr->stonits *= WHTEFFICACY;
		}
		break;
	case 0:				/* see if we can guess #components */
		if (MCsize[mr->dt]) {
			const long	perplane = MCsize[mr->dt]*mr->xres*mr->yres;
			long		flen;
			if (fseek(mr->fin, 0, SEEK_END) < 0)
				goto seek_err;
			flen = ftell(mr->fin);
			mr->ncomp = (flen - mr->start) / perplane;
			if ((mr->ncomp > 0) & (flen - mr->start == perplane*mr->ncomp)) {
				if (fseek(mr->fin, mr->start, SEEK_SET) < 0)
					goto seek_err;
				goto retry_cspace;	/* try again */
			}
		}
		/* fall through */
	default:
		mr->errCode = IREformat;
		strcpy(mr->errMsg, "missing/unsupported number of components");
		MRcloseFile(mr);
		return (ImgReader *)mr;
	}
	mr->fr.xleft = 0; mr->fr.xright = mr->xres;
	mr->fr.ytop = 0; mr->fr.ybottom = 1;
	mr->nr = mr->fr;		/* we're ready to roll */
	return (ImgReader *)mr;
seek_err:
	mr->errCode = IREread;
	strcpy(mr->errMsg, "seek error");
	MRcloseFile(mr);
	return (ImgReader *)mr;
}

/* Get header information and bundle for calling application */
static ImgReadErr
MRgetInfo(ImgReader *ir, ImgInfo *info)
{
	MatrixReader	*mr = (MatrixReader *)ir;
	long		fpos;

	if (!mr || !mr->fin || !info)
		return IREunknown;
	fpos = ftell(mr->fin);
	*info = defImgInfo;		/* rewind to beginning */
	if (fseek(mr->fin, 0, SEEK_SET) < 0) {
		strcpy(mr->errMsg, "cannot rewind to header");
		MRcloseFile(mr);
		return mr->errCode = IREread;
	}
					/* scan info. header */
	if (RHgetInfo(mr->fin, info) < 0) {
		strcpy(mr->errMsg, "error encountered reading info header");
		MRcloseFile(mr);
		return mr->errCode = IREread;
	}
	info->stonits = mr->stonits;
	info->flags |= IIFstonits;
					/* restore file position */
	if (fseek(mr->fin, fpos, SEEK_SET) < 0) {
		strcpy(mr->errMsg, "seek error");
		MRcloseFile(mr);
		return mr->errCode = IREread;
	}
	return mr->errCode = IREnone;	/* all done */
}

/* Read and convert a scanline */
static int
MRscanConvert(MatrixReader *mr, ImgReadBuf *rb, const int y)
{
	char	vbuf[32];
	float	*bufp;
	int	x, i;
					/* our assumptions */
	if ((rb->cs.dtype != IDTfloat) |
			(ImgPixelLen[rb->cs.format] != mr->ncomp)) {
		strcpy(mr->errMsg, "bad parameters to MscanConvert");
		mr->errCode = IREunknown;
		return -1;
	}
					/* seek if indicated */
	if ((y != mr->nr.ytop) | (rb->r.xleft != mr->nr.xleft)) {
		long	seek2;
		if (mr->dt == MCTascii && (y != mr->fr.ytop) |
					(rb->r.xleft != mr->fr.xleft)) {
			seek2 = mr->ncomp*((y - mr->nr.ytop)*mr->xres +
						(rb->r.xleft - mr->nr.xleft));
			if (seek2 < 0) {	/* too expensive */
				strcpy(mr->errMsg, "cannot seek on ASCII input");
				mr->errCode = IREread;
				return -1;
			}
			while (seek2--)
				if (!fgetword(vbuf, sizeof(vbuf), mr->fin))
					goto hitEOF;
		} else {
			seek2 = mr->start + MCsize[mr->dt]*mr->ncomp*
						(y*mr->xres + rb->r.xleft);
			if (fseek(mr->fin, seek2, SEEK_SET) < 0) {
				MRcloseFile(mr);
				strcpy(mr->errMsg, "seek error for requested scanline");
				mr->errCode = IREread;
				return -1;
			}
		}
	}
	bufp = (float *)rb->buf + mr->ncomp *
			((y - rb->r.ytop)/rb->subsample) *
			((rb->r.xright - rb->r.xleft)/rb->subsample);
					/* handle standard case efficiently */
	if ((rb->subsample == 1) & (mr->dt == MCTfloat)) {
	    if (getbinary(bufp, sizeof(float)*mr->ncomp,
			rb->r.xright - rb->r.xleft, mr->fin) < rb->r.xright - rb->r.xleft)
		goto hitEOF;
	    if (mr->swapped)
		swap32((char *)bufp, rb->r.xright - rb->r.xleft);
	    i = (mr->nr.xright - rb->r.xright)*sizeof(float)*mr->ncomp;
	    while (i-- > 0)		/* skip rest of this row */
		if (getc(mr->fin) == EOF)
		    goto hitEOF;
	    goto done;
	}
					/* read each matrix element */
	for (x = rb->r.xleft; x < rb->r.xright; x += rb->subsample) {
	    if (mr->dt == MCTfloat) {
		if (getbinary(bufp, sizeof(float), mr->ncomp, mr->fin) < mr->ncomp)
			goto hitEOF;
		if (mr->swapped)
			swap32((char *)bufp, mr->ncomp);
		bufp += mr->ncomp;
	    } else if (mr->dt == MCTdouble) {
		double	dbuf[MaxPixelLen];
		if (getbinary(dbuf, sizeof(double), mr->ncomp, mr->fin) < mr->ncomp)
			goto hitEOF;
		if (mr->swapped)
			swap64((char *)dbuf, mr->ncomp);
		for (i = 0; i < mr->ncomp; i++)
			*bufp++ = dbuf[i];
	    } else /* mr->dt == MCTascii */ {
		for (i = 0; i < mr->ncomp; i++) {
		    if (!fgetword(vbuf, sizeof(vbuf), mr->fin))
			goto hitEOF;
		    if (!isflt(vbuf)) {
			MRcloseFile(mr);
			sprintf(mr->errMsg, "Mangled real value: '%s'", vbuf);
			mr->errCode = IREformat;
			return -1;
		    }
		    *bufp++ = atof(vbuf);
		}
	    }
					/* skip to next element (or row) */
	    if (x + rb->subsample < rb->r.xright)
		i = rb->subsample;
	    else
		i = mr->nr.xright - x;
	    if (--i <= 0)
		continue;
	    i *= mr->ncomp;
	    if (mr->dt == MCTascii) {
		while (i--)
		    if (!fgetword(vbuf, sizeof(vbuf), mr->fin))
			goto hitEOF;
	    } else {
		i *= MCsize[mr->dt];
		while (i--)
		    if (getc(mr->fin) == EOF)
			goto hitEOF;
	    }
	}
done:					/* aligned to next row at this point */
	mr->nr.ytop = y+1;
	mr->nr.ybottom = y+2;
	return 0;
hitEOF:
	MRcloseFile(mr);
	strcpy(mr->errMsg, "unexpected EOF");
	mr->errCode = IREtruncated;
	return -1;
}

/* Read a rectangle in our color space (not flexible at this point) */
static ImgReadErr
MRreadRec(ImgReader *ir, ImgReadBuf *rb)
{
	MatrixReader	*mr = (MatrixReader *)ir;
	int		i;

	if (!mr || !mr->fin || !rb)
		return IREunknown;
	mr->errCode = IREnone;
	if (!rb->buf || !PmatchColorSpace(&rb->cs, &mr->cs, PICMall)) {
		int	buflen;		/* need new buffer (don't free old) */
		rb->cs = mr->cs;	/* force our color space on reader */
		rb->r.xleft = 0;	/* use scanline-sized buffer */
		rb->r.xright = mr->xres;
		rb->r.ybottom = rb->r.ytop + rb->subsample;
		ImgFixSampling(rb);
		buflen = ImgReadBufLen(rb);
		if (buflen <= 0) {
			strcpy(mr->errMsg, "internal buffer size error");
			return mr->errCode = IREunknown;
		}
		rb->buf = (uby8 *)malloc(buflen);
		if (rb->buf == NULL) {
			strcpy(mr->errMsg, "cannot allocate new read buffer");
			return mr->errCode = IREmemory;
		}
	} else {			/* else align rectangle request */
		ImgFixSampling(rb);
		if (!PlegalRect(&rb->r, mr->xres, mr->yres)) {
			strcpy(mr->errMsg, "bad read rectangle");
			return mr->errCode = IREread;
		}
	}
					/* transfer each scanline */
	for (i = rb->r.ytop; i < rb->r.ybottom; i += rb->subsample)
		if (MRscanConvert(mr, rb, i) < 0)
			break;	
	return mr->errCode;
}

/* Close a normal map and free reader struct */
static void
MRclose(ImgReader *ir)
{
	MatrixReader	*mr = (MatrixReader *)ir;
	if (!mr)
		return;
	MRcloseFile(mr);
	free(mr);
}

/* Interface for matrix image reader */
const ImgReaderInterface IRInterfaceMTX = {
	"Matrix.mtx.mfx.mdx",
	&MRopen, NULL, &MRgetInfo, NULL,
	NULL, &MRreadRec, &MRclose
};
