/*
 *  pimage.h
 *  panlib
 *
 *  C header for Pancine image rendering.
 *
 *  Created by gward on Fri Sept. 7 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#ifndef _PIMAGE_H_
#define _PIMAGE_H_

#include "imgreader.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Image reader interfaces we know about */
extern const ImgReaderInterface	IRInterfaceTIFF;
extern const ImgReaderInterface	IRInterfaceJPEG;
extern const ImgReaderInterface	IRInterfaceRad;
extern const ImgReaderInterface	IRInterfaceEXR;
extern const ImgReaderInterface	IRInterfaceBMP;
extern const ImgReaderInterface IRInterfaceDPT;
extern const ImgReaderInterface IRInterfaceNRM;
extern const ImgReaderInterface IRInterfaceMTX;

/** NULL-terminated list of Pancine image reader interfaces */
#define P_MAXREADER			15
extern const ImgReaderInterface *	PIReaderI[P_MAXREADER+1];

/** Add reader interface to our list, return entry or -1 on error */
extern int		PaddIReaderI(const ImgReaderInterface *iri);

/** Adds all known readers to PIReaderI; returns 1 on success */
extern int		PloadStandardReaders(void);

/** Pancine image reader error reporting function */
extern void		PreportReaderErr(const ImgReader *ir);

/** Move absolute frame, # frames, animation and seek value */
extern int		PabsFrame(int *fmp, int nf, ImgFrameType ft,
					int offs, ImgSeekMode sm);

/** Open an image from the named file, returning an image reader structure
 *  If quiet is true, PopenImageF won't report any errors
 *  If iri is not NULL, PopenImageF will forego the usual search for readers
 */
extern ImgReader *	PopenImageF(const char *fname, int quiet,
					const ImgReaderInterface *iri);

/** Return pointer to final component "fname.ext" in "/directory/fname.ext" */
extern const char *	PgetFilename(const char *path);

/** Return pointer to suffix "ext" in "/directory/fname.ext", NULL if none */
extern const char *	PgetSuffix(const char *fname);

/** Check if suffix matches the given suffix list (e.g., "TIFF.tif") */
extern int		PmatchSuffix(const char *sfx, const char *sfx_list);

/** Find image readers based on file suffix, returning number of candidates
 *  Returns 0 if sfx is NULL
 */
extern int		PtargetReader(const ImgReaderInterface
						*irilist[P_MAXREADER],
					const char *sfx);

/** Allocate with CacheMalloc, which makes room in cache then calls malloc()
 */
extern void		CacheMakeRoom(size_t nbytes);
extern void *		CacheMalloc(size_t nbytes, const char *file, int line);

#define Pmalloc(n)	CacheMalloc(n, __FILE__, __LINE__)
#define Pfree		free

/** Check if image can possibly fit memory */
extern int		PsizeOK(int xr, int yr, int psiz);

/** Allocate an image buffer, adjusting to fit image height/width
 *  If (ib->img != NULL), just adjust image dimensions to match imgAspect
 */
extern int		PnewImage(ImgStruct *ib, double imgAspect);

/** Deallocate an image buffer using member MOMrelease() call */
extern void		PfreeImage(ImgStruct *ib);

/** Link subimage -- assumes ib is uninitialized (or equals ia) */
extern int		PlinkSubimage(ImgStruct *ib, const ImgStruct *ia,
					const ImgRect *r);

#define PlinkImage(ib,ia)	PlinkSubimage(ib,ia,NULL)

/** Link corresponding overlap, placing upper left corner of img2 in img1 */
extern int		PlinkCoverage(ImgStruct *cover1, ImgStruct *cover2,
					const ImgStruct *img1,
					const ImgStruct *img2,
					int xleft, int ytop);

/** Make unique copy of image and eliminate dead space */
extern int		PdelinkImage(ImgStruct *ib);

/** Determine if two images have overlapping pixel memory */
extern int		PimagesOverlap(const ImgStruct *ia,
					const ImgStruct *ib);

/** Assign to entire image the given pixel value */
extern int		PsetImage(ImgStruct *ib, PixelVal pfv);

/** Clear an image with the given raw pixel (black if pf==NULL) */
extern void		PclearImage(ImgStruct *ib, const void *pf);

/** Type definition for rectangle fill callback */
typedef void	PfillRectMethod(ImgStruct *ib, const ImgRect *r, void *udp);

/** Clear the indicated rectangle (black if pf==NULL) */
extern PfillRectMethod	PclearRect;

/** Shift image contents, filling with the given callback function */
extern int		PshiftImageCB(ImgStruct *ib, int dx, int dy,
				PfillRectMethod *fcb, void *udp);

/** Shift image contents, filling with given pixel (no fill if pf==NULL) */
#define PshiftImage(ib,dx,dy,pf)	( (pf) ? \
				PshiftImageCB(ib,dx,dy,PclearRect,pf) : \
				PshiftImageCB(ib,dx,dy,NULL,NULL) )

/** Shift image contents, filling with black */
#define PshiftImageB(ib,dx,dy)	PshiftImageCB(ib,dx,dy,PclearRect,NULL)

/** Copy one image into another at indicated position */
extern int		PcopyImage(ImgStruct *ib, const ImgStruct *ia,
					int xleft, int ytop);

/** Copy the indicated source component to destination image */
extern int		PcopyComponent(ImgStruct *ib, int cb,
					const ImgStruct *ia, int ca);

/** Histogram channels */
typedef enum { PHCred=0, PHCgreen, PHCblue,
		PHCluminance, PHCrgb, PHCrgb3 } PHistoChan;

/** Compute histogram for image:
 *  If minmax[0] < minmax[1] on call, then these are taken as the
 *  historgram range and the histogram is only appended (not cleared).
 *  Otherwise, the minmax values are set based on the image and the
 *  histogram is cleared to all zeroes before the tally begins.
 *  If hlen <= 0, then the minmax calculation is all that happens, and
 *  and the extrema will be updated even if starting minmax[0] < minmax[1].
 *  The actual type of the minmax array is determined by ia->csp->dtype.
 *  If chan == PHCrgb3, the extrema array is dimensioned minmax[3][2],
 *  the histogram is hist[3][histlen].  If chan == PHCluminance and
 *  ia->csp->dtype == IDTfloat, then the histogram is partitioned using
 *  logarithmic steps, i.e.:
 *	hist_bin_low[i] = min*pow(max/min, (double)i/hlen)
 *  and the minimum is set to the smallest positive Y value.
 *  The total count of pixels added into the histogram during the
 *  call is returned, or zero if there was an error.  This total will be
 *  less than ia->xres*ia->yres if pixels are outside the minmax range.
 */
extern unsigned long	PcomputeHisto(void *minmax, unsigned long *hist,
					int hlen, const ImgStruct *ia,
					PHistoChan chan);

/** Compute pixel percentiles:
 *  The pctl[] array specifies the desired percentiles between 0 and 100.
 *  Results are returned in the res[] array, which is dimensioned res[3][n]
 *  if chan == PHCrgb3.
 */
extern int		PcomputePercentiles(void *res, const float *pctl, int n,
					const ImgStruct *ia, PHistoChan chan);

/** Modify image ib to match the histogram of image ia */
extern int		PmatchHisto(ImgStruct *ib, const ImgStruct *ia);

/** Type definition for weighted average callback */
typedef double	PweightingMethod(int x, int y, void *udp);

/** Compute weighted average pixel value (in specified CS) from image ia */
extern PixelVal		PweightedAverageCB(const ImgStruct *ia,
					const ImgColorSpace *csp,
					PweightingMethod *wcb, void *udp);

/** Compute mean pixel value (in the given CS) */
#define PcomputeAverage(ia,csp)		PweightedAverageCB(ia,csp,NULL,NULL)

/** Compute average of image ia weighted by gray values from ib */
extern PixelVal		PweightedAverage(const ImgStruct *ia,
					const ImgStruct *ib,
					const ImgColorSpace *csp);

/*
 *  The rendering routines below only allocate a destination buffer
 *  if ib->img is initialized to NULL before the call.  The dimensions of
 *  the image will be less than or equal to the xres and yres members.
 *  In PsizeImage(), PblurImage(), PwarpImage(), ProtateImage(), and
 *  PsampleImage(), the output image will always match the requested size.
 *  In PconvolveImage(), PdequantizeImage(), PsumImage(), and PmapImage(),
 *  the output image dimensions must match the input dimensions exactly,
 *  and will be set to match if ib->img is NULL.
 */

/** Render (resample) an image from an active reader */
extern int		PrenderImageR(ImgStruct *ib, ImgReader *ir, int quiet);

/** Render (resample) an image from another image */
extern int		PrenderImageI(ImgStruct *ib, const ImgStruct *ia);

/** Resampling basis functions */
enum { PSfast, PSbest, PSnearest, PSbox, PSgaussian, PSlinear, PScubic };

/** Similar to PrenderImageI(), but may alter pixel aspect ratio */
extern int		PsizeImage(ImgStruct *ib, const ImgStruct *ia,
					int basis);

/** Blur image using Gaussian kernel of given radius */
extern int		PblurImage(ImgStruct *ib, const ImgStruct *ia,
					float rad);

/** Structure to hold weights & measures for PsampleImage() */
typedef struct {
	float		wt;		/* input pixel weight (0==end) */
	short		mx, my;		/* x & y measure (input offsets) */
} PweightMeasure;

/** Apply sampling weights & measures from given upper-left start point.
 *  The wma[] array of lists gets tiled over destination image.  The tile's
 *  array dimensions are whres by wvres, and the uppper-left of the first
 *  tile is placed at the axleft and aytop coordinates in input ia.  Steps in
 *  the input for increments in output pixel are given by ahstep and avstep.
 *  Use axleft=aytop=0 & ahstep=avstep=NULL for "onto" image mapping.
 *  Can be used to sample one channel at a time using PHCred, etc.
 *  Fill border as needed with given pixel pfv.
 */
extern int		PsampleImage(ImgStruct *ib, const ImgStruct *ia,
					PweightMeasure *wma[],
					const int whres, const int wvres,
					const double axleft, const double aytop,
					const double ahstep[2],
					const double avstep[2],
					PHistoChan chan, PixelVal pfv);
			
/** Convolve source image with filter kernel */
extern int		PconvolveImage(ImgStruct *ib, const ImgStruct *ia,
					const ImgStruct *ikern);

/** Dilate (or erode) image with a disk of the given radius */
extern int		PdilateImage(ImgStruct *ib, const ImgStruct *ia,
					float rad);

#define	PerodeImage(ib,ia,rad)	PdilateImage(ib,ia,-(rad))

/** Warp an image based on the given grid of source positions in scanline order.
 *  Fill border as needed with given pixel pfv.
 */
extern int		PwarpImage(ImgStruct *ib, const ImgStruct *ia,
					float warpGrid[][2],
					const int whres, const int wvres,
					int basis, PixelVal pfv);

/** Rotate image clockwise about its center, filling border with pfv */
extern int		ProtateImage(ImgStruct *ib, const ImgStruct *ia,
					double degCW,
					int basis, PixelVal pfv);

/** Sum one image into another of the same size (but possibly different type) */
extern int		PsumImage(ImgStruct *ib, const ImgStruct *ia, float sf);

/** Blend between two images using alpha image */
extern int		PblendImages(ImgStruct *ib, const ImgStruct *ia0,
					const ImgStruct *ialpha,
					const ImgStruct *ia1);

/** Type definition for image operation callback (sets ib->csp & allocates) */
typedef int	PimageOp(ImgStruct *ib, const ImgStruct *ia, void *udp);

/** Blend image operation given by callback */
extern int		PblendImageOpCB(ImgStruct *ib, const ImgStruct *ialpha,
						PimageOp *op, void *udp);

/** Byte to float image conversion with debanding filter */
extern int		PdequantizeImage(ImgStruct *ib,
					const ImgStruct *ia, float sf);

/** Byte to float image conversion with variable-width debanding filter */
extern int		PdequantizeImage2(ImgStruct *ib, const ImgStruct *ia,
					float sf, int quant, int minrad);

/** Extract alpha channel from RGBA image */
extern int		PseparateAlpha(ImgStruct *ib, ImgStruct *ialpha,
					const ImgStruct *ia);

/** Integrate alpha channel into RGBA destination */
extern int		PmarryAlpha(ImgStruct *ib, const ImgStruct *ia,
					const ImgStruct *ialpha);

extern int		PrandCQuant;	/* random color quantization? */

/** In situ color space conversion (may delink image) */
extern int		PconvertColorSpace(ImgStruct *ib,
					const ImgColorSpace *dcs, float sf);
					
/** Color space conversion/copy routine (called by routines above) */
extern int		PmapImage(ImgStruct *ib, const ImgStruct *ia, float sf);

/** Type definition for scanline conversion callback (return < 0 to abort) */
typedef int	PscanlineMethod(const uby8 *scan, int len, int y, void *udp);

/** Convert image color space, sending scanlines to given callback function */
extern int		PconvertImageCB(const ImgStruct *ia,
					const ImgColorSpace *dcs, float sf,
					PscanlineMethod *scb, void *udp);

/** Allocate an opaque color conversion structure.
 *  Set cv_wp=0 to skip white point conversion (i.e., no VonKries transform).
 */
extern void *		PcreateColorConv(const ImgColorSpace *cdp,
					const ImgColorSpace *csp,
					float sf, int cv_wp);
#define PfreeColorConv	free

/** Get/add cached simple color conversion with WB (no free routine) */
extern const void *	PgetColorConv(const ImgColorSpace *cdp,
					const ImgColorSpace *csp);

/** Apply color conversion to an array of pixel values */
extern int		PmapPixels(uby8 *pout, const uby8 *pinp, int n,
					const void *cv);

/** Convert color for a single pixel */
extern PixelVal		PconvertPixel(PixelVal pval, const ImgColorSpace *cdp);

/** Get pixel Y value */
#define PgetY(p)	(PconvertPixel(p,&ICS_Y).v.f[0])

/** In situ (modified) gray-world white balance.
  * Set frac to the fraction of full compensation desired.
  * Set orig non-zero if final color space should use original white point.
  */
extern int		PautoWhiteBal(ImgStruct *ib, float frac, int orig);

/** Compute histogram adjustment tone-mapping operator (camera) */
extern int		PhistAdjTMO(ImgStruct *ib, const ImgStruct *ia, double LdDyn);

/** Convert from 8-bit gray to 24-bit RGB (basically a copy function) */
extern void		PgetRGB24fromY8(uby8 *dptr, const uby8 *sptr, int n);

/** Convert from 32-bit RGBA to 24-bit RGB (copy function) */
extern void		PgetRGB24fromRGBA32(uby8 *dptr, const uby8 *sptr, int n);

/** Add alpha channel to 24-bit RGB (copy function) */
extern void		PgetRGBA32fromRGB24(uby8 *dptr, const uby8 *sptr, int n);

/** Convert from 8-bit gray to 32-bit RGBA (copy function) */
extern void		PgetRGBA32fromY8(uby8 *dptr, const uby8 *sptr, int n);

/** Create a MIP map for fast and general resampling */
extern ImgStruct *	PcreateMIPmap(const ImgStruct *im, int linkOK);

/** Create a floating point MIP map from any source */
extern ImgStruct *	PcreateFMIPmap(const ImgStruct *im);

/** Free allocated MIP map */
extern void		PfreeMIPmap(ImgStruct *mip);

/** Render image from MIP map (region) */
extern int		PrenderMIPmap(ImgStruct *ib, const ImgStruct *mip,
					const ImgRect *r);

/** Sample floating point MIP map at a pixel using the given filter radius */
extern int		PsampleFMIPmap(float *res, const ImgStruct *mip,
					float x, float y, float rad);

/** Compute bilateral filter:
 *  Source image ia is modulated by control image ic, which
 *  may be the same as ia, and is assumed the same if set to NULL.
 *  The dimensions of ic must match the larger of ia or ib.
 *  Sigma_s controls the sampling radius measured in input pixels;
 *  sigma_r controls the sampling radius in the signal domain
 *  using normalized units (e.g., an 8-bit quantum is 1./256.).
 */
extern int		PbilatFilter(ImgStruct *ib, const ImgStruct *ia,
					float sigma_s, float sigma_r,
					const ImgStruct *ic);

/** Remove noise from image:
 *  Apply bilateral filter to remove noise from the given image (ia)
 *  which may be the same as destination image (ib).
 *  Input and output image dimensions must match in any case.
 *  The filter sampling radius in pixels is specified as sigma_s.
 *  A linear noise model is specified in noise_c, which is dimensioned
 *  noise_c[nc*2], where nc is the number of components in ia.
 *  The variance due to noise for component c at
 *  level v equals (noise_c[2*c] + v*noise_c[2*c+1]).
 *  A non-float image is considered to have a 0-1 range.
 */
extern int		PdenoiseImage(ImgStruct *ib, const ImgStruct *ia,
					float sigma_s, const float *noise_c);

/** Blend two sections of a panorama:
 *	Input is two overlapping images in matching color spaces.
 *	The rectangle "anch" has its left-top coordinate
 *	set to a feature in "inpa" that matches "inpb" when the
 *	right-bottom corner of "inpb" is at the right-bottom corner of
 *	"anch".  (This is not necessarily a legal rectangle.)
 *	On completion, PblendPano() creates the image "blnd" to cover
 *	the overlapping region between "inpa" and "inpb".
 */
extern int		PblendPano(ImgStruct *blnd, const ImgRect *anch,
				const ImgStruct *inpa, const ImgStruct *inpb);

#ifdef __cplusplus
}			/* extern "C" */
#endif

#endif /* ! _PIMAGE_H_ */
