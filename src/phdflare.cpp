/*
 *  phdflare.cpp
 *  panlib
 *
 *  Remove lens flare from high dynamic-range image.
 *
 *  Created by Gregory Ward on Thu Feb 06 2003.
 *  Copyright (c) 20017 Anyhere Software. All rights reserved.
 *
 */

#include <math.h>
#include "pancine.h"
#include "phdrimg.h"
#include "gaussjord.h"

#ifdef phdf_debug
#include "imgwriter.h"
#endif

#ifndef PHDFRES
#define PHDFRES		200			// working resolution
#endif

#ifndef PHDFTHRESH
#define PHDFTHRESH	500.			// hot:minimum threshold ratio
#endif

#ifndef PHDFR2MIN
#define PHDFR2MIN	8			// radius^2 for PSF plateau
#endif

#define PHDFR2MAX	(PHDFRES*PHDFRES/4)	// maximum radius^2

#ifndef PHDFORDER
#define PHDFORDER	3			// order of PSF polynomial fit
#endif

#ifndef PHDFRMV
#define PHDFRMV		0.98f			// fraction of flare to remove
#endif

// Our flare removal class
class FlareFilter {
	ImgStruct	*im;			// image being filtered
	float		(*yval)(const uby8 *);	// luminance conversion
	ImgStruct	csi;			// reduced color image
	ImgStruct	smi;			// reduced gray image
	ABitMap2	hot;			// hot pixel bitmap
	float		hotavg[8];		// average hot color
	PPolynomial	psfp;			// point spread function fit
	void		Free() {
				hot.NewBitMap(0, 0);
				PfreeImage(&csi);
				PfreeImage(&smi);
				im = NULL;
			}
	double		PSF(int r2) const {
				if (r2 < PHDFR2MIN) r2 = PHDFR2MIN;
				else if (r2 > PHDFR2MAX) r2 = PHDFR2MAX;
				double	pv = psfp.Eval(sqrt(1./(double)r2));
				if (pv < 0) pv = 0;
				return pv;
			}
	float		LocalMin(int x, int y, int *cp = NULL) const;
	bool		Resample();
	bool		FindSources();
	bool		ComputePSF();
	bool		ApplyPSF();
public:
			FlareFilter() {
				im = NULL;
				yval = NULL;
				csi.img = NULL;
				smi.img = NULL;
				psfp.coef[psfp.degree=0] = 0;
			}
			~FlareFilter() {
				Free();
			}
	bool		RemoveFlare(ImgStruct *ims,
				int (*pbar)(const char *, int) = NULL);
};

static float
rgb2y(const uby8 *p)
{
	const float *	rgb = (const float *)p;
	
	return 0.21264f*rgb[0] + 0.71517f*rgb[1] + 0.07219f*rgb[2];
}

static float
xyz2y(const uby8 *p)
{
	return ((const float *)p)[1];
}

static float
y2y(const uby8 *p)
{
	return *(const float *)p;
}

// Get minimum value at this reduced pixel location
float
FlareFilter::LocalMin(int x, int y, int *cp) const
{
	const int	nc = ImgPixelLen[im->csp->format];
	const int	x0 = x*im->xres/smi.xres;
	const int	y0 = y*im->yres/smi.yres;
	const int	x1 = (x+1)*im->xres/smi.xres - 1;
	const int	y1 = (y+1)*im->yres/smi.yres - 1;
	float		vmin = 1e10, vrmin = 1e10;
	int		xl, yl, i;
	const float *	ip;
					// get minimum component ratio
	for (yl = y0 ; yl < y1; yl++) {
		ip = (const float *)PpixPtr(im,x0,yl);
		for (xl = x0; xl < x1; xl++)
			for (i = 0; i < nc; i++, ip++)
				if (*ip < vrmin*hotavg[i]) {
					vmin = *ip;
					vrmin = vmin/hotavg[i];
					if (cp != NULL)
						*cp = i;
				}
	}
	return vmin;
}

// Reduce original to grayscale sample image
bool
FlareFilter::Resample()
{
	if ((im->xres < PHDFRES) | (im->yres < PHDFRES))
		return false;
					// compute reduced image
	PfreeImage(&csi);
	csi.csp = im->csp;
	csi.xres = csi.yres = PHDFRES;
	if (!PrenderImageI(&csi, im))
		return false;
					// convert to grayscale
	switch (im->csp->format) {
	case IPFrgb:
		yval = &rgb2y;
		break;
	case IPFxyz:
		yval = &xyz2y;
		break;
	case IPFy:
		yval = &y2y;
		break;
	default:
		DMESG(DMCparameter, "Unsupported image format in Resample");
		return false;
	}
	PfreeImage(&smi);
	if (yval == &y2y) {
		PlinkImage(&smi, &csi);
	} else {
		const int	psiz = ImgPixelSize(csi.csp);
		int		x, y;
		smi.csp = &ICS_Y;
		smi.xres = csi.xres; smi.yres = csi.yres;
		if (!PnewImage(&smi, .0))
			return false;
		for (y = smi.yres; y--; ) {
			float *		dp = (float *)ProwPtr(&smi,y);
			const uby8 *	sp = ProwPtr(&csi,y);
			for (x = smi.xres; x--; sp += psiz)
				*dp++ = (*yval)(sp);
		}
	}
#ifdef phdf_debug
{
ImgWriteBuf iwb;
PsetWriteBuf(&iwb, &smi);
IWInterfaceTIFF.WriteImage("/tmp/small.tif", &iwb);
}
#endif
	return true;
}

// Find flare sources
bool
FlareFilter::FindSources()
{
	float		thresh;
	int		x, y, i;
	const float	*lp, *ip;
					// find image minimum
	thresh = 1e20f;
	for (y = smi.yres; y--; ) {
		lp = (float *)ProwPtr(&smi,y);
		for (x = smi.xres; x--; lp++)
			if (*lp < thresh)
				thresh = *lp;
	}
					// use multiplier for threshold
	thresh *= (float)PHDFTHRESH;
					// sources are samples over threshold
	const int	nc = ImgPixelLen[csi.csp->format];
	for (i = nc; i--; )
		hotavg[i] = 0;
	hot.NewBitMap(smi.xres, smi.yres);
	for (y = smi.yres; y--; ) {
		lp = (const float *)ProwPtr(&smi,y);
		for (x = 0; x < smi.xres; x++, lp++)
			if (*lp >= thresh) {
				hot.Set(x, y);
				ip = (const float *)PpixPtr(&csi,x,y);
				for (i = nc; i--; )
					hotavg[i] += ip[i];
			}
	}
					// compute average hot color
	int	nhot = (int)hot.SumTotal();
	if (!nhot)
		return false;
	for (i = nc; i--; )
		hotavg[i] /= (float)nhot;
					// # source samples < 95% of image
	return (nhot < (long)smi.xres*smi.yres*95L/100L);
}

#define PHDMAXPT	(PHDFRES/5)     // based on number of radii, below

// Compute point spread function
bool
FlareFilter::ComputePSF()
{
	const int	brdr = 3;
	double		Amat[PHDMAXPT][PHDFORDER+1];
	double		bvec[PHDMAXPT];
	int		r2min, r2max, br, r2;
	int		x, y, xo, yo, x1, y1;
	float		v, vmin, vminlast;
	ImgStruct	fli;
	int		npts, i, k;
	const float *	ap;
	const float *	wp;
					// clear previous solution
	psfp.coef[psfp.degree=0] = 0;
					// allocate low-res flare image
	fli = smi;
	fli.img = NULL;
	if (!PnewImage(&fli, .0))
		return false;
					// find min. for fit at each radius
	ABitMap2	pMin(fli.xres, fli.yres);
	vminlast = 1.f;
	for (r2min = PHDFR2MIN; (r2max = r2min*3/2) <= PHDFR2MAX; r2min = r2max) {
		br = (int)(sqrt((double)r2max) + .99999);
		PclearImage(&fli, NULL);
					// draw annulus around each hot pixel
		for (x = y = 0; hot.Find(&x, &y); x++) {
			
			v = *(const float *)PpixPtr(&smi,x,y);
			for (yo = -br; yo <= br; yo++) {
			    y1 = y + yo;
			    if (y1 < 0) continue;
			    if (y1 >= fli.yres) break;
			    for (xo = -br; xo <= br; xo++) {
				x1 = x + xo;
				if (x1 < 0) continue;
				if (x1 >= fli.xres) break;
				r2 = xo*xo + yo*yo;
				if ((r2min <= r2) & (r2 < r2max))
					*(float *)PpixPtr(&fli,x1,y1) += v;
			    }
			}
		}
		vmin = 2.f;		// find minimum for PSF inside borders
		for (y = fli.yres-brdr; y-- > brdr; ) {
			ap = (const float *)ProwPtr(&fli,y) + brdr;
			wp = (const float *)ProwPtr(&smi,y) + brdr;
			for (x = brdr; x < fli.xres-brdr; x++, ap++, wp++) {
				if (*ap <= 1e-7f)
					continue;
				v = *wp / *ap;
				if (v >= vmin)
					continue;
				if (hot.Check(x, y))
					continue;
				if (pMin.Check(x, y))
					continue;
				x1 = x; y1 = y;
				vmin = v;
			}
		}
#ifdef phdf_debug
fprintf(stderr, "vmin=%.2e at (%d,%d) for annulus radius %d%s\n",
	vmin, x1, y1, (r2min+r2max)/2, (vmin>=vminlast)?" (rejected)":"");
#endif
		if (vmin >= vminlast)	// enforce monotonic decreasing PSF
			continue;
					// record point
		pMin.Set(x1, y1);
		vminlast = vmin;
	}
#ifdef phdf_debug
WriteBitMap2(hot, "/tmp/hotpix.bmp");
WriteBitMap2(pMin, "/tmp/minpix.bmp");
#endif
	const int	nc = ImgPixelLen[csi.csp->format];
	npts = 0;			// compute overdetermined matrix
	for (x1 = y1 = 0; pMin.Find(&x1, &y1); x1++) {
		DASSERT(npts < PHDMAXPT);
		bvec[npts] = LocalMin(x1, y1, &k);
#ifdef phdf_debug
fprintf(stderr, "local min=%.2e at (%d,%d), chan %d\n",
		bvec[npts], x1, y1, k);
#endif
		for (i = 0; i <= PHDFORDER; i++)
			Amat[npts][i] = 0;
		for (y = csi.yres; y--; ) {
		    ap = (const float *)ProwPtr(&csi,y) + k;
		    for (x = 0; x < csi.xres; x++, ap += nc) {
			double	d = (double)((x1-x)*(x1-x) + (y1-y)*(y1-y));
			if (d < PHDFR2MIN) d = PHDFR2MIN;
			else if (d > PHDFR2MAX) d = PHDFR2MAX;
			d = 1./sqrt(d);
			double	t = *ap;
			for (i = 0; i <= PHDFORDER; i++, t *= d)
				Amat[npts][i] += t;
		    }
		}
		++npts;
	}
	if (npts <= PHDFORDER) {
		DMESG(DMCtrace, "Too few sample points");
		return false;
	}
					// find least-squares polynomial fit
	if (!GJsolveLeastSq(PHDFORDER+1, npts,
				(double *)Amat, (double *)bvec, psfp.coef)) {
		DMESG(DMCtrace, "Cannot fit PSF");
		return false;
	}
	psfp.degree = PHDFORDER;
	return true;
}

// Apply point spread function to remove flare
bool
FlareFilter::ApplyPSF()
{
	const float	xrf = (float)smi.xres/(float)im->xres;
	const float	yrf = (float)smi.yres/(float)im->yres;
	ImgStruct	fli, ffunc;
	int		y, x0, y0;
	int		x, i;
	float		*fp, *fp0;
					// compute PSF template w/ quad symmetry
	ffunc.csp = &ICS_Y;
	ffunc.xres = 2*smi.xres + 1; ffunc.yres = 2*smi.yres + 1;
	ffunc.img = NULL;
	if (!PnewImage(&ffunc, .0))
		return false;
	for (y = 0; y <= smi.yres; y++) {
		const int	y2 = (y-smi.yres)*(y-smi.yres);
		fp = (float *)ProwPtr(&ffunc,y);
		for (x = 0; x <= smi.xres; x++)
			*fp++ = PSF(y2 + (x-smi.xres)*(x-smi.xres));
		--fp;
		for (x = 1; x <= smi.xres; x++)
			fp[x] = fp[-x];
	}
	for (y = 1; y <= smi.yres; y++) {
		fp0 = (float *)ProwPtr(&ffunc,smi.yres-y);
		fp = (float *)ProwPtr(&ffunc,smi.yres+y);
		for (x = ffunc.xres; x--; )
			*fp++ = *fp0++;
	}
#if 0
					// equalization at center
	psfsum = 0;
	for (y = smi.yres; y--; ) {
		fp = (float *)ProwPtr(&ffunc,y);
		for (x = smi.xres; x--; )
			psfsum += *fp++;
	}
	if (psfsum > .25)
		return false;
	*(float *)PpixPtr(&ffunc,smi.xres,smi.yres) -= psfsum;
#endif
					// compute low-res flare image
	fli = csi;
	fli.img = NULL;
	if (!PsetImage(&fli, Pblack)) {
		PfreeImage(&ffunc);
		return false;
	}
	const int	nc = ImgPixelLen[fli.csp->format];
	const int	rowlen = fli.rowsize/sizeof(float);
	for (x0 = y0 = 0; hot.Find(&x0, &y0); x0++) {
		const float * const	cv = (const float *)PpixPtr(&csi,x0,y0);
		for (y = fli.yres; y--; ) {
			fp0 = (float *)PpixPtr(&ffunc,csi.xres-x0,csi.yres+y-y0);
			fp = (float *)ProwPtr(&fli,y);
			for (x = fli.xres; x--; fp0++, fp += nc)
				for (i = nc; i--; )
					fp[i] += *fp0 * cv[i];
		}
	    }
	PfreeImage(&ffunc);
#ifdef phdf_debug
{
ImgWriteBuf iwb;
PsetWriteBuf(&iwb, &fli);
IWInterfaceTIFF.WriteImage("/tmp/flarec.tif", &iwb);
}
#endif
	float	frac = PHDFRMV;		// prevent underflow
	for (y = 0; y < fli.yres; y++) {
		fp0 = (float *)ProwPtr(&csi,y);
		fp = (float *)ProwPtr(&fli,y);
		for (x = 0; x < fli.xres; x++, fp0 += nc, fp += nc)
			for (i = nc; i--; )
				if (fp[i]*frac > fp0[i])
					frac = 0.999f * fp0[i] / fp[i];
	}
					// subtract interpolated flare from original
	for (y = 0; y < im->yres; y++) {
		float	yf = yrf*((float)y+.5f) - .5f;
		y0 = (int)yf;
		if (y0 >= fli.yres-1) y0 = fli.yres-2;
		yf -= (float)y0;
		fp = (float *)ProwPtr(im,y);
		for (x = 0; x < im->xres; x++, fp += nc) {
			float	xf = xrf*((float)x+.5f) - .5f;
			x0 = (int)xf;
			if (x0 >= fli.xres-1) x0 = fli.xres-2;
			xf -= (float)x0;
			fp0 = (float *)ProwPtr(&fli,y0) + x0*nc;
			for (i = nc; i--; ) {
				fp[i] -= frac*((1.f-xf)*(1.f-yf)*fp0[i] +
						xf*(1.f-yf)*fp0[nc + i] +
						(1.f-xf)*yf*fp0[rowlen + i] +
						xf*yf*fp0[nc + rowlen + i]);
				if (fp[i] < 0)
					fp[i] = 0;
			}
		}
	}
					// clean up
	PfreeImage(&fli);
	return true;
}

// Make necessary calls to remove lens flare
bool
FlareFilter::RemoveFlare(ImgStruct *ims, int (*pbar)(const char *, int))
{
	if (ims == NULL || (ims->csp == NULL) | (ims->img == NULL))
		return false;
	if (ims->csp->dtype != IDTfloat) {
		DMESG(DMCparameter, "PHDremoveFlare called with LDR image");
		return false;
	}
	im = ims;
	DMESG(DMCtrace, "Begin flare removal...");
	if (pbar != NULL && !(*pbar)("Finding flare sources", 0)) {
		DMESG(DMCtrace, "User cancelled RemoveFlare");
		return false;
	}
	if (!Resample())
		goto failure;
	DMESG(DMCtrace, "Image resampled");
	if (pbar != NULL && !(*pbar)("Finding flare sources", 20)) {
		DMESG(DMCtrace, "User cancelled RemoveFlare");
		return false;
	}
	if (!FindSources())
		goto failure;
	sprintf(dmessage_buf, "Sources found (%ld/%ld points)",
			(long)hot.SumTotal(), (long)hot.Width()*hot.Height());
	DMESG(DMCtrace, dmessage_buf);
	if (pbar != NULL && !(*pbar)("Computing flare spread", 30)) {
		DMESG(DMCtrace, "User cancelled RemoveFlare");
		return false;
	}
	if (!ComputePSF())
		goto failure;
	strcpy(dmessage_buf, "PSF computed (");
	for (int i = psfp.degree; i > 0; i--)
		sprintf(strchr(dmessage_buf,'\0'),
				"%.2e*r^-%d + ", psfp.coef[i], i);
	sprintf(strchr(dmessage_buf,'\0'), "%.2e)", psfp.coef[0]);
	DMESG(DMCtrace, dmessage_buf);
	if (pbar != NULL && !(*pbar)("Removing lens flare", 60)) {
		DMESG(DMCtrace, "User cancelled RemoveFlare");
		return false;
	}
	if (!ApplyPSF())
		goto failure;
	DMESG(DMCtrace, "Applied PSF to remove flare");
	if (pbar != NULL)
		(*pbar)("Flare removed", 100);
	DMESG(DMCtrace, "Flare removed");
	return true;
failure:
	if (pbar != NULL)
		(*pbar)("Flare removal failed", 100);
	DMESG(DMCwarning, "Flare removal failed");
	return false;
}

// Compute and remove lens flare from high dynamic range image
bool
PHDremoveFlare(ImgStruct *ims, int (*pbar)(const char *, int))
{
	FlareFilter	myFilt;

	return myFilt.RemoveFlare(ims, pbar);
}
