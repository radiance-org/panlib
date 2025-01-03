/*
 *  pfeatures.h
 *  pan
 *
 *  Depends on "panimage.h"
 *
 *  Definitions for image feature-matching.
 *
 *  Created by Greg Ward on 2/24/15.
 *  Copyright 2015 Anyhere Software. All rights reserved.
 *
 */

#ifndef _PFEATURES_H_
#define _PFEATURES_H_

#include "abitmap.h"

#define	BRIEF_BITS	384		// bits per BRIEF descriptor
#define BRIEF_DIA	64		// feature patch diameter
#define BRIEF_DTHRESH	(BRIEF_BITS/5)	// default threshold for match

// An image feature descriptor and location
struct ImgFeature {
	int		x, y;		// descriptor position in image
	ABitMap		descr;		// BRIEF descriptor
			ImgFeature() {
				x = y = -1;
			}
			// Take data from another feature (destructive)
	bool		Take(ImgFeature *srcp) {
				x = srcp->x; y = srcp->y;
				srcp->x = srcp->y = -1;
				return descr.Take(&srcp->descr);
			}
			// Compute distance^2 from another feature location
	int		Dist2(const ImgFeature &f1) const {
				return (x-f1.x)*(x-f1.x) + (y-f1.y)*(y-f1.y);
			}
			// Compute feature distance
	int		Fdist(const ImgFeature &f2) const;
};

// Normalize image to emphasize edges, returning float gray map
bool			PnormEdges(PanImage *imp, const PanImage &im, float blur = 10.f);

// Compute feature (Harris corner) bitmap
extern bool		PdetectFeatures(ABitMap2 *fmap, const PanImage &im, float blur = 0);

// A list of image features
class PanFeatureList {
private:
	int		xres, yres;	// input image resolution
	int		nfeats;		// number of features found
	ImgFeature *	feat;		// feature array
public:
			PanFeatureList() {
				xres=yres=0;
				nfeats=0; feat=0;
			};
			PanFeatureList(const PanImage &im, float blur = 0) {
				xres=yres=0;
				nfeats=0; feat=0;
				FindFeatures(im, blur);
			}
			PanFeatureList(const PanFeatureList &orig) {
				nfeats=0; feat=0;
				*this = orig;
			}
			~PanFeatureList() {
				delete [] feat;
			}
			// Assign features indicated by the given bitmap
	int		SetFeatures(const PanImage &im, ABitMap2 fmap, float blur = 0);
			// Identify features in the given image, return count
	int		FindFeatures(const PanImage &im, float blur = 0) {
				ABitMap2	featMap;
				if (!PdetectFeatures(&featMap, im, blur)) return -1;
				return SetFeatures(im, featMap, blur);
			}
			// Remove possibly confounding features from our list (within rad)
	int		CullSimilar(const double rad = 0, const int thresh = BRIEF_DTHRESH);
			// Clear our feature list
	void		Clear() {
				delete [] feat;
				nfeats=0; feat=0;
				xres=yres=0;
			};
			// Return feature caunt
	int		Count() const {
				return nfeats;
			}
			// Return width
	int		Width() const {
				return xres;
			}
			// Return height
	int		Height() const {
				return yres;
			}
			// Get the indexed feature
	const ImgFeature *
			GetFeature(int n) const {
				if ((n < 0) | (n >= nfeats)) return 0;
				return &feat[n];
			}
	const ImgFeature &
			operator[](int n) const {
				static const ImgFeature	dummy;
				if ((n < 0) | (n >= nfeats)) return dummy;
				return feat[n];
			}
			// Find best-matching n features (within threshold)
	int		BestFeatures(int fndx[], int n, const ImgFeature &mf,
					const int thresh = BRIEF_BITS) const;
	int		BestFeatures(const ImgFeature *farr[], int n,
					const ImgFeature &mf,
					const int thresh = BRIEF_BITS) const {
				if (!farr | (n <= 0)) return 0;
				int *	fndx = new int [n];
				const int	nm = BestFeatures(fndx, n, mf, thresh);
				while (n > nm) farr[--n] = 0;
				while (n--) farr[n] = &feat[fndx[n]];
				delete [] fndx;
				return nm;
			}
			// Find best-matching feature (within threshold)
	const ImgFeature *
			MatchFeature(const ImgFeature &mf,
					const int thresh = BRIEF_DTHRESH) const {
				int	n;
				if (BestFeatures(&n, 1, mf, thresh))
					return &feat[n];
				return 0;
			}
			// Find nearest n features (within threshold)
	int		NearFeatures(int fndx[], int n, const ImgFeature &mf,
					const int thresh = BRIEF_DTHRESH) const;
	int		NearFeatures(const ImgFeature *farr[], int n,
					const ImgFeature &mf,
					const int thresh = BRIEF_DTHRESH) const {
				if (!farr | (n <= 0)) return 0;
				int *	fndx = new int [n];
				const int	nm = NearFeatures(fndx, n, mf, thresh);
				while (n > nm) farr[--n] = 0;
				while (n--) farr[n] = &feat[fndx[n]];
				delete [] fndx;
				return nm;
			}
			// Find nearest feature (within threshold)
	const ImgFeature *
			NearestFeature(const ImgFeature &mf,
					const int thresh = BRIEF_DTHRESH) const {
				int	n;
				if (NearFeatures(&n, 1, mf, thresh))
					return &feat[n];
				return 0;
			}
			// Copy operator
	PanFeatureList &operator=(const PanFeatureList &orig);
};

// Compare two image features and return Hamming distance
inline int
ImgFeature::Fdist(const ImgFeature &f1) const
{
	static int	nthread = 0;	// cheap semaphore
	static ABitMap	statBM;

	if (!descr.Length() | !f1.descr.Length())
		return BRIEF_BITS+1;

	if (++nthread == 1) {		// avoids memory allocation
		statBM = descr;
		statBM ^= f1.descr;
		int	tot = statBM.SumTotal();
		--nthread;
		return tot;
	}
	--nthread;			// multi-thread needs automatic storage
	ABitMap	bxor(descr);
	bxor ^= f1.descr;
	return bxor.SumTotal();
}

// Type definition for callback to trace which features were used
typedef void	PfeatureUsedMethod(const ImgFeature &dstF, const ImgFeature &srcF,
					bool used, void *udp);

// Compute warp grid (source positions) to register feature sets
extern bool		PregisterWarp(float warpGrid[][2], int whres, int wvres,
					const PanFeatureList &dstFL,
					const PanFeatureList &srcFL,
					const int thresh = BRIEF_DTHRESH,
					PfeatureUsedMethod *cb=0, void *udp=0);

#endif	// ! _PFEATURES_H_
