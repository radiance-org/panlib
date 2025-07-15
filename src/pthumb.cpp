/*
 *  pthumb.cpp
 *  panlib
 *
 *  Pancine thumbnail cache manager
 *
 *  Created by gward on Thu Jun 07 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#include <errno.h>
#include <math.h>
#include "pancine.h"
#include "pthumb.h"
#include "imgwriter.h"

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

/*
 *  Thumbnails are cached in JPEG contact sheets to reduce disk
 *  fragmentation and improve compression performance.  The stored
 *  thumbnails are sized so they may be viewed at a higher
 *  resolution when desired.  When the standard resolution is
 *  requested (PBigThumbWidth/PThumbDen), the contact image is
 *  downsampled as it is loaded by the JPEG reader.  The cache
 *  directory holds a database with a standard name and as many
 *  contact sheets as will fit within the allotted maximum.
 *  The memory cache manager (cache.h) is used to control when
 *  thumbnails are swapped in and out of memory while the
 *  thumbnail manager is running.  This should be transparent
 *  to the controlling application.  Total memory use is controlled
 *  by the memory cache manager, and disk usage is set using the
 *  thumbnail manager on a per directory basis.
 */

#ifndef PTM_NRECS
#define	PTM_NRECS	1024			// thumbnail DB read length
#endif
#ifndef PTM_MINPURGE
#define PTM_MINPURGE	3			// minimum % to purge at once
#endif
#ifndef PTM_QUALITY
#define PTM_QUALITY	75			// thumbnail JPEG quality
#endif
#define PTM_MAXCONTACT	16			// max. PThumbNRows*PThumbNCols

#define	PTnow()		((int32)((time(0)/60)&0x7fffffff))

const char	PancineThumbnailerName[] = "Pancine Thumbnailer";

const char	PancineThumbnailDBFile[] = "tncache.adb";

const int	PBigThumbWidth = 320;		// large thumbnail width
int		PThumbNRows = 2;		// thumb rows per image
int		PThumbNCols = 3;		// thumb columns per image

const PTDBFieldInfo	PTDBFInfo;		// Thumbnailer field information

// Compute kilobytes from bytes
inline static int
Kbytes(int32 bytes)
{
	return (bytes + 1023) >> 10;
}

// Create a small, generic thumbnail image for the given record
bool
PgenericThumbnail(ImgStruct *tn, const DBRecord &, PTStatus status)
{
	extern const int	PThumbIconWidth;
	extern const uby8	PThumbIconData[];
	ImgStruct *		tndest = tn;
	ImgStruct		tnnew;
						// check arguments
	if (!PmatchColorSpace(tn->csp, &ICS_sRGB, PICMptype))
		return false;
	if ((tn->xres < PThumbIconWidth) | (tn->yres < PThumbIconWidth)) {
		tnnew.csp = tn->csp;
		tnnew.img = NULL;
		tnnew.mbase = NULL;
		tn = &tnnew;
	}
	tn->xres = tn->yres = PThumbIconWidth;
	if (!PnewImage(tn, .0))
		return false;
						// copy gray image -> RGB
	const uby8 *	sp = PThumbIconData;
	for (int y = 0; y < tn->yres; y++) {
		uby8 *	dp = ProwPtr(tn, y);
		int		i;
		if (status == PTSgood)
			for (i = tn->xres; i--; ) {
				*dp++ = *sp;
				*dp++ = *sp;
				*dp++ = *sp++;
			}
		else if (status == PTSnothm)
			for (i = tn->xres; i--; ) {
				*dp++ = *sp >> 1;
				*dp++ = *sp;
				*dp++ = *sp++ >> 1;
			}
		else /* PTSmissing or PTSbad */
			for (i = tn->xres; i--; ) {
				*dp++ = *sp;
				*dp++ = *sp >> 1;
				*dp++ = *sp++ >> 2;
			}
	}
	if (status == PTSbad)			// add shotgun blast
		for (int n = PThumbIconWidth*PThumbIconWidth/50; n--; ) {
			double	r = 1./RAND_MAX*rand(); r *= r;
			double	p = 2.*M_PI/RAND_MAX*rand();
			int	x = int(.5*PThumbIconWidth*(1. + r*cos(p)));
			int	y = int(.5*PThumbIconWidth*(1. + r*sin(p)));
			uby8 *	dp = PpixPtr(tn, x, y);
			*dp++ = 0xff;
			*dp++ = 0xff;
			*dp = 0xff;
		}
	if (tndest != tn)
		return PrenderImageI(tndest, tn);
	return true;
}

// Initialize thumbnailer field object (order must match PTDBFieldID)
PTDBFieldInfo::PTDBFieldInfo()
{
	static const char	eemsg[] = "Enumeration error in PTDBFieldInfo";

	if (Field("NBytes", DBDTint, "Original image size (bytes)",
			DBFFrequired|DBFFindex) != PTDBFnbytes)
		DMESG(DMCassert, eemsg);
	if (Field("HashVal", DBDTint, "Original hash value",
			DBFFrequired|DBFFindex) != PTDBFhashval)
		DMESG(DMCassert, eemsg);
	if (Field("ImageID", DBDTint, "Contact image ID",
			DBFFrequired) != PTDBFimageid)
		DMESG(DMCassert, eemsg);
	if (Field("ImageBytes", DBDTint, "Contact image subsize (bytes)",
			DBFFrequired) != PTDBFimagebytes)
		DMESG(DMCassert, eemsg);
	if (Field("Rectangle", DBDTint, "Contact image subregion",
			DBFFarray|DBFFrequired) != PTDBFrectangle)
		DMESG(DMCassert, eemsg);
	if (Field("LastUse", DBDTint, "Time at last usage",
			DBFFrequired) != PTDBFlastuse)
		DMESG(DMCassert, eemsg);
}

// Check thumbnail header values
bool
PTDBHeader::CheckValues() const
{
	if (strcmp(soft, PancineThumbnailerName) != 0)
		return false;
	int	fieldMap[DB_MAXFIELD];
	if (MapIDs(fieldMap, &PTDBFInfo))
		return false;			// no field adjustments allowed
	return DBHeader::CheckValues();
}

// Thumbnail contact image file name generator
inline static char *
PcontactName(char *fname, const char *dir, int32 iid)
{
	sprintf(fname, "%s%cC%09d.jpg", dir, DIRSEP, (int)iid);
	return fname;
}

// Tally cache size for block of records and get next image ID
static int
countCache(const DBRecordList &rl, void *cdp)
{
	int32	*	nIDp = (int32 *)cdp;
	int32		blockBytes = 0;
	int		i = rl.GetSize();
	
	while (i--) {
		int32	iv;
		if (rl.Get(i)[PTDBFimageid].Get(&iv) && iv > *nIDp)
			*nIDp = iv;
		blockBytes += rl.Get(i)[PTDBFimagebytes].GetInt();
	}
	return Kbytes(blockBytes);
}

// Set new cache directory for thumbnail manager
bool
PThumbnailer::SetCacheDirectory(const char *dname)
{
	if (dname != NULL && !strcmp(dname, cacheDir))
		return true;
	if (!SyncCache())			// write out in-core structures
		DMESG(DMCdata, "Thumbnail cache sync failed");
	nextImage.ClearImage();
	if (!tdb.Close())			// close DB
		DMESG(DMCdata, "Cannot close thumbnail database");
	tnHead.Attach(NULL);
	dbfs.close();
	cacheDir[0] = '\0';
	nThumbs = 0;
	if (dname == NULL || !*dname)
		return false;
	char	dbname[1024];
	sprintf(dbname, "%s%c%s", dname, DIRSEP, PancineThumbnailDBFile);
	dbfs.clear();
	dbfs.open(dbname, ios::in|ios::out|ios::binary);
	if (!dbfs.is_open()) {			// try adding trunc flag
		dbfs.clear();
		dbfs.open(dbname, ios::in|ios::out|ios::trunc|ios::binary);
	}
	if (!dbfs.is_open()) {
		DMESGF(DMCresource, "Cannot open thumbnail database %s", dbname);
		return false;
	}
	if (!tdb.Init(&tnHead, &dbfs)) {
		DMESGF(DMCdata, "Corrupt thumbnail database %s", dbname);
		dbfs.close();
		return false;
	}
	nextImgID = 0;				// get cache size and last ID
	cacheSiz = tdb.ForEachBlock(&countCache, &nextImgID);
	++nextImgID;
	nThumbs = tdb.TotalRecords();
	strcpy(cacheDir, dname);		// cache database open
	return true;
}

// Write out contact image in progress and check cache limit
bool
PThumbnailer::SyncCache()
{
	if (!cacheDir[0])
		return true;
	if (!nextImage.IsEmpty()) {	// write out contact sheet
		char	imgfile[1024];
		PcontactName(imgfile, cacheDir, nextImgID);
		int32	siz = nextImage.WriteImage(imgfile);
		if (siz <= 0) {
			DMESGF(DMCdata, "Error writing contact image '%s'",
					imgfile);
			return false;
		}
					// update thumbnail records with size
		DBSearch	dbs(&tdb);
		DBRecord	query(&tnHead);
		int		nt = nextImage.GetNThumbs();
		DBRecord *	rlist = new DBRecord [nt];
		if (!query.SetField(PTDBFimageid, nextImgID) ||
				!dbs.AddSelector(query))
			DMESG(DMCassert, "AddSelector failed");
		nt = dbs.FindRecords(rlist, nt);
		DTEST(nt < nextImage.GetNThumbs(), DMCwarning,
				"Missing thumbnail records in SyncCache");
		tdb.DeleteRecords(rlist, nt);
		for (int i = nt; i--; )
			rlist[i].SetField(PTDBFimagebytes,
				(int32)(i ? siz/nt : siz/nt + (siz%nt)));
		tdb.AddRecords(rlist, nt);
		delete [] rlist;
		cacheSiz += Kbytes(siz);
	}
	PurgeCache();			// reassert cache limit
	return tdb.Sync();		// sync database
}

// Priority struct needed by PurgeCache
struct PTNPurgeObj {
	int32		id;		// contact image id if it exists
	int32		lu;		// last use for this image
	int32		siz;		// size of image (bytes)
};

// Holder for deletion list
struct PTNDelList {
	ABitMap		delMap;		// deletion bitmap
	DBRecordList	toDel;		// records to delete
			PTNDelList(uint32 dlen = 0) : delMap(dlen) {}
};

// Sort LRU first, then largest to smallest for images with same lu value
// Relegate unassigned records (id == 0) to end of list
static int
PTNpurgeCmp(const void *p1, const void *p2)
{
	const PTNPurgeObj *	po1 = (const PTNPurgeObj *)p1;
	const PTNPurgeObj *	po2 = (const PTNPurgeObj *)p2;
	if (po1->id <= 0) return (po2->id > 0);
	if (po2->id <= 0) return -1;
	if (po1->lu < po2->lu) return -1;
	if (po1->lu > po2->lu) return 1;
	if (po1->siz < po2->siz) return 1;
	if (po1->siz > po2->siz) return -1;
	return 0;
}

// Accumulate sizes and last use measures for block
static int
PTNaccumCache(const DBRecordList &rl, void *lp)
{
	PTNPurgeObj *	parr = (PTNPurgeObj *)lp;
	int		psiz = (parr++)->id;
	int		i = rl.GetSize();
	
	while (i--) {
		int32	iid, lu;
		if (!rl.Get(i)[PTDBFimageid].Get(&iid) || iid >= psiz)
			continue;
		parr[iid].id = iid;
		parr[iid].siz += rl.Get(i)[PTDBFimagebytes].GetInt();
		if (rl.Get(i)[PTDBFlastuse].Get(&lu) && lu > parr[iid].lu)
			parr[iid].lu = lu;
	}
	return 1;
}

// Gather records for deletion
static int
PTNgatherDoomed(const DBRecordList &rl, void *lp)
{
	PTNDelList *	dlp = (PTNDelList *)lp;
	int		i = rl.GetSize();
	
	while (i--)
		if (dlp->delMap.Check(rl.Get(i)[PTDBFimageid].GetInt()))
			dlp->toDel.NewRecord() = rl.Get(i);
	return 1;
}

// Purge thumbnail cache directory
long
PThumbnailer::PurgeCache()
{
	if (!cacheDir[0] | (cacheMax < 0) | (cacheSiz <= cacheMax))
		return 0;

	PTNPurgeObj *	parr;
						// get last use values & sizes
	parr = (PTNPurgeObj *)calloc(nextImgID+1, sizeof(PTNPurgeObj));
	if (parr == NULL) {
		DMESG(DMCmemory, "Cannot allocate priority array in PurgeCache");
		return 0;
	}
	parr[0].id = nextImgID;			// first record marks size
	tdb.ForEachBlock(&PTNaccumCache, parr);
	++parr;
						// sort with LRU first
	qsort(parr, nextImgID, sizeof(PTNPurgeObj), PTNpurgeCmp);
						// choose files to purge
	long	minPurge = cacheSiz*PTM_MINPURGE/100 + 1024;
	long	nKtopurge = cacheSiz - cacheMax;
	long	nKpurged = 0;
	int	n;
	if (!cacheMax)				// check special cases
		nKtopurge += minPurge;
	else if (nKtopurge < minPurge)
		nKtopurge = minPurge;
	for (n = 0; (n < nextImgID) & (nKpurged < nKtopurge); n++) {
		if (!parr[n].id)
			break;
		nKpurged += Kbytes(parr[n].siz);
	}
	PTNDelList	dList(nextImgID);	// delete files and records
	while (n-- > 0) {
		char	imgfile[1024];		// remove file
		PcontactName(imgfile, cacheDir, parr[n].id);
		if (remove(imgfile) < 0) {
			if (errno != ENOENT) {
				DMESGF(DMCresource,
					"Unlink '%s' failed", imgfile);
				free(parr);
				return 0;
			}
			DMESGF(DMCwarning, "Sheet '%s' missing", imgfile);
		}
		dList.delMap.Set(parr[n].id);	// flag contact ID removed
		cacheSiz -= Kbytes(parr[n].siz);
	}
	--parr;
	free(parr);				// done with list
						// delete corresponding records
	tdb.ForEachBlock(&PTNgatherDoomed, &dList);
	tdb.DeleteRecordList(dList.toDel);
	if ((nThumbs = tdb.TotalRecords()) <= 0) {
		cacheSiz = 0;			// reset if cleared all
		nextImage.ClearImage();
		nextImgID = 1;
	} else if (cacheSiz < 0)
		cacheSiz = 0;
	sprintf(dmessage_buf, "%s %ld Kbytes from thumbnail cache",
			cacheSiz ? "Purged" : "Emptied", nKpurged);
	DMESG(DMCinfo, dmessage_buf);
	return nKpurged;			// all done
}

// Find a thumbnail that matches the given Pancine image record
bool
PThumbnailer::FindThumbnail(DBRecord *tr, const DBRecord &pr)
{
	if (tr != NULL)
		tr->Clear();
	if (nThumbs <= 0)
		return false;
	DBSearch	dbs(&tdb);
	DBRecord	query(&tnHead);
					// match nbytes & hashval
	if (!query.SetField(PTDBFnbytes, PDBgetField(pr, PDBFnbytes)) ||
		    !query.SetField(PTDBFhashval, PDBgetField(pr, PDBFhashval)))
		return false;
	if (!dbs.AddSelector(query))
		DMESG(DMCassert, "AddSelector failed");
	if (tr == NULL)
		tr = &query;
	return dbs.FindRecords(tr);	// search for thumbnail record
}

// Get final thumbnail for the given Pancine image record
PTStatus
PThumbnailer::GetThumbnail(ImgStruct *tn, const DBRecord &pr, ImgReader *ir)
{
	if (!cacheDir[0]) {		// no cache?
		PgenericThumbnail(tn, pr, PTSnothm);
		return PTSnothm;
	}
					// check for existing thumbnail
	if (LoadThumbnail(tn, pr))
		return PTSgood;
					// generate new thumbnail
	bool	myreader = (ir == NULL);
	if (myreader)
		ir = PopenImageD(pr, true);
	if (ir == NULL) {
		PgenericThumbnail(tn, pr, PTSmissing);
		return PTSmissing;
	}
	PTStatus	st = NewThumbnail(tn, pr, ir);
	if (myreader)
		IRclose(ir);
	if (st != PTSgood)
		PgenericThumbnail(tn, pr, st);
	return st;
}

// Get quick thumbnail for the given Pancine image record, true if final
bool
PThumbnailer::QuickThumb(ImgStruct *tn, const DBRecord &pr, ImgReader *ir)
{
					// check for existing thumbnail
	if (LoadThumbnail(tn, pr))
		return true;
					// try image reader
	bool	myreader = (ir == NULL);
	if (myreader)
		ir = PopenImageD(pr, true);
	bool	ok = (ir != NULL);
	if (ok) {
					// may work, may not
		ImgReadErr	re = IRgetThumbnail(ir, tn);
		if (myreader)
			IRclose(ir);
		if (re == IREnone) {
			FixThumb(tn, pr);
			return false;	// not final
		}
	}
					// use an icon as our last resort
	PgenericThumbnail(tn, pr, (ok ? PTSgood : PTSmissing));
	return false;
}

// Apply crop border to rendered thumbnail
void
PThumbnailer::ApplyCrop(ImgStruct *tn, const DBRecord &pr) const
{
	int32			cr[4], oxr, oyr;
	if (PDBgetField(pr,PDBFcrop).Get(cr, 4) != 4 ||
			!PDBgetField(pr,PDBFxsize).Get(&oxr) ||
			!PDBgetField(pr,PDBFysize).Get(&oyr))
		return;
	uby8		myRGB[3] = {bRGB[0], bRGB[1], bRGB[2]};
	ImgRect		imreg;
	ImgRect		clreg;
	imreg.xleft = cr[0]*tn->xres/oxr;
	imreg.ytop = cr[1]*tn->yres/oyr;
	imreg.xright = cr[2]*tn->xres/oxr;
	imreg.ybottom = cr[3]*tn->yres/oyr;
	clreg.ytop = 0; clreg.ybottom = imreg.ytop;
	clreg.xleft = 0; clreg.xright = tn->xres;
	PclearRect(tn, &clreg, &myRGB);
	clreg.ytop = imreg.ybottom; clreg.ybottom = tn->yres;
	PclearRect(tn, &clreg, &myRGB);
	clreg.ytop = imreg.ytop; clreg.ybottom = imreg.ybottom;
	clreg.xleft = 0; clreg.xright = imreg.xleft;
	PclearRect(tn, &clreg, &myRGB);
	clreg.xleft = imreg.xright; clreg.xright = tn->xres;
	PclearRect(tn, &clreg, &myRGB);
}

// Apply orientation to rendered thumbnail
void
PThumbnailer::ApplyOrient(ImgStruct *tn, const DBRecord &pr)
{
	POrient		ior(&pr);
	ior.SetOrigSize(tn->xres, tn->yres);
	ior.OrientImage(tn);
}

// Delete thumbnail from cache (necessarily removes other thumbs in contact)
void
PThumbnailer::DeleteThumbnail(const DBRecord &pr)
{
	DBRecord	tr;
	if (!FindThumbnail(&tr, pr))
		return;
	int32	iid;
	if (!tr[PTDBFimageid].Get(&iid) || (iid <= 0) | (iid > nextImgID))
		return;
					// find all thumbs in contact image
	DBRecord	rlist[PTM_MAXCONTACT];
	DBSearch	dbs(&tdb);
	DBRecord	query(&tnHead);
	if (!query.SetField(PTDBFimageid, iid) || !dbs.AddSelector(query))
		DMESG(DMCassert, "AddSelector failed in DeleteThumbnail");
	int		nr = dbs.FindRecords(rlist, PTM_MAXCONTACT);
					// delete them all
	tdb.DeleteRecords(rlist, nr);
	if (iid == nextImgID) {
		nextImage.ClearImage();	// restart current contact image
		return;
	}
	char	imgfile[1024];		// else remove contact image file
	PcontactName(imgfile, cacheDir, iid);
	if (remove(imgfile) < 0) {
		DMESGF(DMCresource, "Cannot remove '%s'", imgfile);
		return;
	}
	int32	totsiz = 0;
	for (int i = nr; i--; )
		totsiz += rlist[i][PTDBFimagebytes].GetInt();
	cacheSiz -= Kbytes(totsiz);
}

// Load the thumbnail for the given Pancine image record from memory or disk
bool
PThumbnailer::LoadThumbnail(ImgStruct *tn, const DBRecord &pr)
{
	if (!cacheDir[0])
		return false;
	DBRecord	tr;		// find thumbnail
	if (!FindThumbnail(&tr, pr))
		return false;
					// get contact image ID
	int32	iid;
	int32	rect[4];
	if (!tr[PTDBFimageid].Get(&iid) || iid <= 0)
		return false;
	if (tr[PTDBFrectangle].Get(rect, 4) != 4)
		return false;
					// update last use
	if (!tdb.DeleteRecord(tr) ||
			!tr.SetField(PTDBFlastuse, PTnow()) ||
			!tdb.AddRecord(tr))
		DMESG(DMCdata, "Cannot alter thumbnail last use value");

					// check if it's "in progress"
	if (iid == nextImgID) {
		if (!nextImage.GetThumb(tn, rect))
			return false;
		FixThumb(tn, pr);
		return true;
	}
					// else search our list
	const int	li = DenomShift(tn);
	int		nempties = 0;
	PThumbList *	tnl;
	for (tnl = images[li]; tnl != NULL; tnl = tnl->next) {
		if (tnl->GetID() == iid)
			break;
		nempties += tnl->Swapped();
	}
	if (tnl == NULL) {		// need to load image
		if (images[li] == NULL)
			tnl = new PThumbList(iid, cacheDir, 1<<li);
		else
			tnl = new PThumbList(iid, images[li]);
		images[li] = tnl;
		if (nempties > 50) {	// do a little housekeeping
			PThumbList *	tnn = tnl->next;
			while (tnn != NULL) {
				if (tnn->Swapped()) {
					tnl->next = tnn->next;
					tnn->next = NULL;
					delete tnn;
					tnn = tnl;
				}
				tnn = (tnl = tnn)->next;
			}
			tnl = images[li];
			DMESGF(DMCtrace, "Cleared %d swapped contacts", nempties);
		}
	}
	if (tnl->GetThumb(tn, rect)) {	// load thumbnail
		FixThumb(tn, pr);
		return true;
	}
					// contact sheet vanished...
	DMESGF(DMCwarning, "Missing contact sheet (%d)", (int)iid);
	DBSearch	dbs(&tdb);
	DBRecord	query(&tnHead);
	if (!query.SetField(PTDBFimageid, iid) || !dbs.AddSelector(query))
		DMESG(DMCassert, "AddSelector failed");
	DBRecord	rlist[PTM_MAXCONTACT];
	int		nr = dbs.FindRecords(rlist, PTM_MAXCONTACT);
	int32		totsiz = 0;
					// delete all references
	nThumbs -= tdb.DeleteRecords(rlist, nr);
	while (nr--)			// correct cache size
		totsiz += rlist[nr][PTDBFimagebytes].GetInt();
	cacheSiz -= Kbytes(totsiz);
	return false;
}

// Create new thumbnail and add it to our database (private)
PTStatus
PThumbnailer::NewThumbnail(ImgStruct *tn, const DBRecord &pr, ImgReader *ir)
{
	if (nextImage.IsFull()) {
		if (!SyncCache()) {	// write out full contact sheet
			DMESG(DMCdata, "Cannot sync thumbnail cache");
			return PTSbad;
		}
		nextImage.ClearImage();
		++nextImgID;		// and we're on to the next one
	}
	int32		rect[4];	// create new thumbnail
	PTStatus	st = nextImage.AddThumb(tn, rect, ir);
	if (st != PTSgood)
		return st;
					// create new thumbnail record
	DBRecord	newrec(&tnHead);
	if (!newrec.SetField(PTDBFnbytes, PDBgetField(pr, PDBFnbytes)) ||
			!newrec.SetField(PTDBFhashval,
					PDBgetField(pr, PDBFhashval)) ||
			!newrec.SetField(PTDBFimageid, nextImgID) ||
			!newrec.SetField(PTDBFrectangle, rect, 4) ||
			!newrec.SetField(PTDBFlastuse, PTnow()) ||
			!tdb.AddRecord(newrec)) {
		DMESG(DMCdata, "Error adding new thumbnail record");
		return PTSbad;
	}
	++nThumbs;			// increment thumbnail total
	FixThumb(tn, pr);
	return PTSgood;
}

// Clear contact image and prepare for output
void
PThumbImage::ClearImage()
{
	int	oldsiz = 3 * ncols*colwidth * nrows*rowheight;

	ncols = PThumbNCols;		nrows = PThumbNRows;
	colwidth = PBigThumbWidth;	rowheight = PBigThumbWidth;

	int	newsiz = 3 * ncols*colwidth * nrows*rowheight;

	if (newsiz != oldsiz) {
		delete [] img;
		img = new uby8 [newsiz];
	}
	
	memset(img, 0x7f, newsiz);

	nextSlot = 0;
}

// Check chromaticity coordinates to be sure they're not too extreme
static bool
chromasOK(const float chroma[4][2])
{
	for (int i = 4; i--; )
		if ((chroma[i][0] < .01f) | (chroma[i][0] > .8f) |
				(chroma[i][1] < .01f) | (chroma[i][1] > .7f))
			return false;
	return true;
}

// Add thumbnail from reader to our contact sheet image in progress
PTStatus
PThumbImage::AddThumb(ImgStruct *tn, int32 rect[4], ImgReader *ir)
{
	if (IsFull())
		return PTSbad;		// caller is messing with me!
	if (img == NULL)
		ClearImage();		// allocate & initialize
					// assign rectangle
	rect[0] = (nextSlot % ncols) * colwidth;
	rect[1] = nextSlot/ncols * rowheight;
	rect[2] = rect[0] + colwidth;
	rect[3] = rect[1] + rowheight;
	ImgStruct	thm;		// render RGB thumbnail
	SetThumbnail(&thm, rect);
	if (tnCS.format == ir->cs.format && chromasOK(ir->cs.chroma))
		memcpy(tnCS.chroma, ir->cs.chroma, sizeof(tnCS.chroma));
					// writes directly on contact image
	if (!PrenderImageR(&thm, ir, true))
		return ir->errCode ? PTSbad : PTSnothm;
					// correct thumbnail dimensions
	rect[2] = rect[0] + thm.xres;
	rect[3] = rect[1] + thm.yres;
					// increment slot counter
	nextSlot++;
					// render copy for caller
	return (PTStatus)PrenderImageI(tn, &thm);
}

// Get thumbnail region (x0 y0 x1 y1) from contact image in progress
bool
PThumbImage::GetThumb(ImgStruct *tn, const int32 rect[4])
{
	if (!PTcheckRect(rect, ncols*colwidth, nrows*rowheight))
		return false;
	
	ImgStruct	thm;
	SetThumbnail(&thm, rect);
	return PrenderImageI(tn, &thm);
}

// Write out a contact sheet and report the number of bytes written
int32
PThumbImage::WriteImage(const char *fname) const
{
	if (IsEmpty())
		return 0;
	ImgWriteBuf	wbuf;
	wbuf.xres = ncols*colwidth;
	wbuf.yres = nrows*rowheight;
	wbuf.rowsize = wbuf.xres*3;
	wbuf.pixAspect = 1;
	wbuf.csp = &ICS_sRGB;
	wbuf.info.flags = IIFquality;
	wbuf.info.quality = PTM_QUALITY;
	wbuf.img = img;
	return (*IWInterfaceJPEG.WriteImage)(fname, &wbuf);
}

// Get thumbnail region (x0 y0 x1 y1) from contact image
bool
PThumbList::GetThumb(ImgStruct *tn, const int32 rect[4])
{
	if (GetCacheObject() == NULL)	// make sure image is loaded
		return false;
	if (!PTcheckRect(rect, oWidth, oHeight)) {
		ReleaseObject(objCache);
		return false;
	}
	ImgStruct	thm;		// render copy of thumbnail
	int		cwidth = oWidth / denom;
	int		rx = rect[0] / denom;
	int		ry = rect[1] / denom;
	thm.csp = &ICS_sRGB;
	thm.xres = (rect[2] - rect[0]) / denom;
	if (thm.xres <= 0) thm.xres = 1;
	thm.yres = (rect[3] - rect[1]) / denom;
	if (thm.yres <= 0) thm.yres = 1;
	thm.rowsize = cwidth*3;
	thm.mbase = NULL;
	thm.img = (uby8 *)objCache + (size_t)ry*thm.rowsize + rx*3;
	bool	ok = PrenderImageI(tn, &thm);
	ReleaseObject(objCache);	// release object memory
	return ok;
}

// Restore contact image from disk cache
bool
PThumbList::RestoreMemory()
{
	if (IsResident())
		return true;
	char	imgfile[1024];
	PcontactName(imgfile, cdir, imgID);
	ImgReader *	irdr = PopenImageF(imgfile, true, &IRInterfaceJPEG);
	if (irdr == NULL)
		return false;		// doesn't exist or not JPEG
	if ((irdr->pixAspect < .995f) | (irdr->pixAspect > 1.005f))
		return false;		// should never happen!!
	oWidth = irdr->xres;
	oHeight = irdr->yres;
	ImgStruct	csheet;
	csheet.csp = &ICS_sRGB;
	csheet.xres = oWidth / denom;
	csheet.yres = oHeight / denom;
	csheet.rowsize = csheet.xres*3;
	objSize = (size_t)csheet.yres*csheet.rowsize;
	objCache = Cmalloc(objSize);
	if (objCache == NULL) {
		objSize = 0;
		return false;
	}
	csheet.img = (uby8 *)objCache;
	csheet.mbase = NULL;
					// read (and resample) image
	if (!PrenderImageR(&csheet, irdr, true)) {
		IRclose(irdr);
		Cfree(objCache); objCache = NULL;
		objSize = 0;
		return false;
	}
	IRclose(irdr);			// close reader and return success
	return true;
}

// Free memory associated with contact image
void
PThumbList::FreeMemory()
{
	if (!IsPurgeable())
		return;
	Cfree(objCache);
	objCache = NULL;
}
