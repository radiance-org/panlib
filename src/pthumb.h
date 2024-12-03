/*
 *  pthumb.h
 *  panlib
 *
 *  Depends on "dbase.h" and "imgreader.h"
 *
 *  Thumbnail cache manager class
 *
 *  Created by gward on Thu Jun 07 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#ifndef _PTHUMB_H_
#define	_PTHUMB_H_

#ifndef _DBASE_H_
#include "dbase.h"
#endif
#ifndef _IMGREADER_H_
#include "imgreader.h"
#endif
#include <fstream>
using std::fstream;

#ifndef PT_NSHIFT
#define PT_NSHIFT	4				// max. denom shift + 1
#endif

extern const int	PBigThumbWidth;			// large thumbnail width
extern int		PThumbNRows;			// # contact rows
extern int		PThumbNCols;			// # contact columns

// Thumbnail database field enumeration
enum PTDBFieldID {
			PTDBFnbytes,			// orig. image size
			PTDBFhashval,			// orig. hash value
			PTDBFimageid,			// contact image ID
			PTDBFimagebytes,		// subregion byte usage
			PTDBFrectangle,			// subregion boundaries
			PTDBFlastuse,			// time at last use
			PTDBFend			// terminator
};

// Thumbnail status
enum PTStatus {PTSbad=0, PTSgood, PTSnothm, PTSmissing};

// Database field information object for Pancine thumbnailer
class PTDBFieldInfo : public DBFieldInfo {
public:
			PTDBFieldInfo();
};
extern const PTDBFieldInfo	PTDBFInfo;		// only need one

extern const char	PancineThumbnailerName[];

// Pancine thumbnail cache database header class
class PTDBHeader : public DBHeader {
public:
			// Set all values -- be sure to clear first!
	virtual void	SetDefaults() {
				DBHeader::SetDefaults();
				strcpy(soft, PancineThumbnailerName);
				*dynamic_cast<DBFieldInfo *>(this) =
					*dynamic_cast<const DBFieldInfo *>(&PTDBFInfo);
			}
			// Check that our DB header values are legal
	virtual bool	CheckValues() const;
};

// Class to hold contact image in progress
class PThumbImage {
private:
	ImgColorSpace	tnCS;				// thumbnail color space
	int		ncols, nrows;			// thumbnail array size
	int		colwidth, rowheight;		// thumbnail dimensions
	uby8 *		img;				// RGB buffer
	int		nextSlot;			// next available slot
public:
			PThumbImage() {
				PcopyCS(&tnCS, &ICS_sRGB);
				ncols = nrows = 1;
				colwidth = rowheight = 0;
				img = NULL;
				nextSlot = 0;
			}
			~PThumbImage() {
				delete [] img;
			}
	int		GetNThumbs() {
				return nextSlot;
			}
	bool		IsEmpty() const {
				return (nextSlot == 0);
			}
	bool		IsFull() const {
				return (nextSlot >= ncols*nrows);
			}
	void		SetThumbnail(ImgStruct *tn,
					const int32 rect[4]) {
				tn->csp = &tnCS;
				tn->xres = rect[2] - rect[0];
				tn->yres = rect[3] - rect[1];
				tn->rowsize = ncols*colwidth*3;
				tn->mbase = NULL;
				tn->img = img + (size_t)rect[1]*tn->rowsize + rect[0]*3;
			}
	void		ClearImage();
	int32		WriteImage(const char *fname) const;
	PTStatus	AddThumb(ImgStruct *tn, int32 rect[4], ImgReader *ir);
	bool		GetThumb(ImgStruct *tn, const int32 rect[4]);
};

class PThumbnailer;					// forward decl.

// Class to cache contact images in memory
class PThumbList : CacheObject {
private:
	int32		imgID;				// contact image ID
	const char *	cdir;				// image directory
	int		denom;				// subsampling rate
	int		oWidth, oHeight;		// original image size
	virtual bool	RestoreMemory();
	virtual void	FreeMemory();
public:
	PThumbList *	next;				// next in list
			// call for first list element
			PThumbList(int iid, const char *cd, int dnm) {
				imgID = iid; cdir = cd;
				denom = dnm; next = NULL;
			}
			// subsequent call to grow list (pln != NULL)
			PThumbList(int iid, PThumbList *pln) {
				imgID = iid; cdir = pln->cdir;
				denom = pln->denom; next = pln;
			}
			~PThumbList() {
				delete next;
				FreeMemory();
			}
			// Get name of subclass
	virtual const char *
			ClassName() const {
				return "PThumbList";
			}
	bool		Swapped() const {
				return (objCache == NULL);
			}
	int32		GetID() const {
				return imgID;
			}
	bool		GetThumb(ImgStruct *tn, const int32 rect[4]);
};

// Pancine thumbnail cache statistics
struct PThumbnailerStats {
	const char *	dName;				// directory name
	long		kAlloc;				// cache Kbytes maximum
	long		kUsed;				// cache Kbytes used
	int		nThumbs;			// number of thumbnails
};

// Pancine thumbnail cache manager class
class PThumbnailer {
protected:
	char		cacheDir[512];			// cache directory name
	long		cacheMax;			// cache limit (Kbytes)
	fstream		dbfs;				// thumbnail DB stream
	PTDBHeader	tnHead;				// thumbnail DB header
	DBAccess	tdb;				// thumbnail DB access
	long		cacheSiz;			// cache size (Kbytes)
	int		nThumbs;			// number of records
	int32		nextImgID;			// next contact ID
	PThumbList *	images[PT_NSHIFT];		// contact image lists
	PThumbImage	nextImage;			// big image in progress
	static int	DenomShift(const ImgStruct *tn, int ow = 0);
	bool		SyncCache();
	long		PurgeCache();
	bool		FindThumbnail(DBRecord *tr, const DBRecord &pr);
	PTStatus	NewThumbnail(ImgStruct *tn, const DBRecord &pr,
					ImgReader *ir);
	bool		LoadThumbnail(ImgStruct *tn, const DBRecord &pr);
	void		ApplyCrop(ImgStruct *tn, const DBRecord &pr) const;
	static void	ApplyOrient(ImgStruct *tn, const DBRecord &pr);
	void		FixThumb(ImgStruct *tn, const DBRecord &pr) const {
				if (doCrop) ApplyCrop(tn, pr);
				if (doOrient) ApplyOrient(tn, pr);
			}
public:
	bool		doCrop;				// apply cropping?
	bool		doOrient;			// reorient thumbs?
	uby8		bRGB[3];			// crop background RGB
			PThumbnailer(const char *dname = NULL) {
				doCrop = true; doOrient = false;
				bRGB[0] = bRGB[1] = bRGB[2] = 200;
				cacheMax = -1; cacheDir[0] = '\0';
				nThumbs = 0;
				for (int i = PT_NSHIFT; i--; ) images[i] = NULL;
				SetCacheDirectory(dname);
			}
			~PThumbnailer() {
				SyncCache();
				for (int i = PT_NSHIFT; i--; ) delete images[i];
			}
			// Set new cache directory
	bool		SetCacheDirectory(const char *dname);
			// Set cache limit to specified # Kbytes (-1==inf)
	void		SetCacheSize(long csiz = -1) {
				bool	purg = (csiz >= 0 &&
					(csiz < cacheMax) | (cacheMax < 0));
				cacheMax = csiz;
				if (purg) { PurgeCache(); tdb.Sync(); }
			}
			// Get cache directory statistics
	bool		GetCacheStats(PThumbnailerStats *ts) const {
				if (!cacheDir[0]) return false;
				ts->dName = cacheDir;
				ts->kAlloc = cacheMax;
				ts->kUsed = cacheSiz;
				ts->nThumbs = nThumbs;
				return true;
			}
			// Get final thumbnail for image, PTSbad on failure
	PTStatus	GetThumbnail(ImgStruct *tn, const DBRecord &pr,
					ImgReader *ir = NULL);
			// Get quick thumbnail for image, true if final
	bool		QuickThumb(ImgStruct *tn, const DBRecord &pr,
					ImgReader *ir = NULL);
			// Delete thumbnail from cache (& fellow contacts)
	void		DeleteThumbnail(const DBRecord &pr);
};

// Figure out JPEG denominator shift for requested thumbnail size
inline int
PThumbnailer::DenomShift(const ImgStruct *tn, int ow)
{
	const int	siz = (tn->xres > tn->yres) ? tn->xres : tn->yres;
	int	i = PT_NSHIFT;
	if (ow <= 0) ow = PBigThumbWidth;
	while (--i)
		if (ow >> i >= siz)
			return i;
	return 0;
}

// Check to see if thumbnail rectangle fits on contact sheet
inline bool
PTcheckRect(const int32 rect[4], int xres, int yres) {
	return (rect[0] >= 0) & (rect[0] < rect[2]) & (rect[2] <= xres) &
		(rect[1] >= 0) & (rect[1] < rect[3]) & (rect[3] <= yres);
}

// Create a small, generic thumbnail image for the given record
extern bool	PgenericThumbnail(ImgStruct *tn, const DBRecord &pr,
					PTStatus status = PTSgood);

#endif	// ! _PTHUMB_H_
