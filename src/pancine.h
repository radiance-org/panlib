/*
 *  pancine.h
 *  panlib
 *
 *  Includes and classes for Pancine imaging library
 *
 *  Created by gward on Wed Jun 06 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#ifndef _PANCINE_H_
#define _PANCINE_H_

#include <iostream>
#include <fstream>
using namespace std;
#include <stdlib.h>
#include "system.h"
#include "rtio.h"
#include "dmessage.h"
#include "pimage.h"
#include "dbase.h"

#define	P_NHASHBITS	31			// # bits in file hash value

#define	P_FIRSTYR	1980			// earliest year we handle

// Pancine database field enumeration
enum PDBFieldID {
			PDBFfile,			// file name
			PDBFlocation,			// file directory
			PDBFnbytes,			// file length
			PDBFhashval,			// file hash value
			PDBFformat,			// file format
			PDBFmodtime,			// file modified time
			PDBFcapdate,			// capture date
			PDBFgmtdate,			// precise GMT date/time
			PDBFlatitude,			// degrees North latitude
			PDBFlongitude,			// degrees East longitude
			PDBFaltitude,			// meters above sea level
			PDBFbearing,			// sighting deg. E of N
			PDBFsubject,			// subject
			PDBFtitle,			// title
			PDBFowner,			// owner
			PDBFhistory,			// history
			PDBFalbum,			// album(s)
			PDBFposition,			// album position(s)
			PDBFcaption,			// caption(s)
			PDBFkeyword,			// search keyword(s)
			PDBFcomment,			// comment
			PDBFsource,			// image source class
			PDBFmake,			// device manufacturer
			PDBFmodel,			// device model
			PDBFversion,			// software revision
			PDBFquality,			// compression quality
			PDBFcrop,			// crop x0 y0 x1 y1
			PDBForient,			// image rotation
			PDBFflip,			// post-rotation flip
			PDBFredeye,			// red-eye removal boxes
			PDBFspotexp,			// spot-metered boxes
			PDBFxsize,			// original width
			PDBFysize,			// original height
			PDBFxdensity,			// horizontal ppi
			PDBFydensity,			// vertical ppi
			PDBFxangle,			// horizontal degrees
			PDBFyangle,			// vertical degrees
			PDBFview,			// Radiance view spec.
			PDBFstonits,			// sample-to-Nits factor
			PDBFexptime,			// exposure seconds
			PDBFaperture,			// f-stop
			PDBFasa,			// ASA
			PDBFflash,			// flash setting
			PDBFwhitebal,			// illuminant
			PDBFfocus,			// focus distance
			PDBFend				// terminator
};

// Pancine database field information object
class PDBFieldInfo : public DBFieldInfo {
public:
			PDBFieldInfo();
};
extern const PDBFieldInfo	PDBFInfo;		// only need one

// Pancine principal database search fields
class PDBSearchFieldInfo : public DBFieldInfo {
public:
			PDBSearchFieldInfo();
};
extern const PDBFieldInfo	PDBSFInfo;		// only need one

extern const char	PancineSoftwareName[];
extern const char	PDBUniqBoolName[];
extern const char	PDBLastSubjName[];
extern const char	PDBLastAlbmName[];
extern const char	PDBLastOwnrName[];
extern const char	PDBLastKeywName[];

					// for tracking recent settings...
#define kNRecentSubjects	7
#define kNRecentAlbums		7
#define kNRecentOwners		5
#define kNRecentKeywords	9
#define kValueLen		64	// upper average value length

// Pancine database header class
class PDBHeader : public DBHeader {
private:
	char		lastSubj[kNRecentSubjects*kValueLen];
	char		lastAlbm[kNRecentAlbums*kValueLen];
	char		lastOwnr[kNRecentOwners*kValueLen];
	char		lastKeyw[kNRecentKeywords*kValueLen];
	int		uniq;		// only allow one record per image?
public:
			PDBHeader() {
				uniq = 0;
				hvars = new DBHeaderLink(PDBUniqBoolName,
							&uniq, 1, hvars);
				hvars = new DBHeaderLink(PDBLastSubjName,
						lastSubj, sizeof(lastSubj), hvars);
				hvars = new DBHeaderLink(PDBLastAlbmName,
						lastAlbm, sizeof(lastAlbm), hvars);
				hvars = new DBHeaderLink(PDBLastOwnrName,
						lastOwnr, sizeof(lastOwnr), hvars);
				hvars = new DBHeaderLink(PDBLastKeywName,
						lastKeyw, sizeof(lastKeyw), hvars);
			}
			// Get value of unique record boolean
	bool		Unique() const {
				return (uniq != 0);
			}
			// Set value of unique record boolean
	void		SetUnique(bool val) {
				hvars->SetVal(&uniq, (int)val);
			}
			// Set all values -- be sure to clear first!
	virtual void	SetDefaults() {
				DBHeader::SetDefaults();
				strcpy(soft, PancineSoftwareName);
				*dynamic_cast<DBFieldInfo *>(this) =
					*dynamic_cast<const DBFieldInfo *>(&PDBFInfo);
				hvars->FindVar(PDBLastSubjName)->nset = 0;
				hvars->FindVar(PDBLastAlbmName)->nset = 0;
				hvars->FindVar(PDBLastOwnrName)->nset = 0;
				hvars->FindVar(PDBLastKeywName)->nset = 0;
				SetUnique(true);
			}
			// Get nth most recent value
	const char *	GetRecent(const char *nm, int n = 0) const {
				const DBHeaderLink *	hl = hvars->GetVar(nm);
				if (hl == NULL) return NULL;
				return hl->GetSVal(n);
			}
			// Get number of recent values
	int		GetNRecent(const char *nm) const {
				const DBHeaderLink *	hl = hvars->GetVar(nm);
				if (hl == NULL || hl->vtype != DBDTstring)
					return -1;
				return hl->nset;
			}
			// Find index of a recent value, or -1 if no match
	int		FindRecent(const char *nm, const char *val) const {
				const DBHeaderLink *	hl = hvars->GetVar(nm);
				if (hl == NULL) return -1;
				return hl->FindSVal(val);
			}
			// Insert recent value at the specified point in list
	bool		InsertRecent(const char *nm, const char *val, int n = 0) {
				DBHeaderLink *	hl = hvars->FindVar(nm);
				if (hl == NULL) return false;
				return hl->InsertSVal(val, n);
			}
			// Add recent value to the end of our list
	bool		AddRecent(const char *nm, const char *val) {
				return InsertRecent(nm, val, -1);
			}
			// Delete indicated recent value (-1 == last)
	void		DeleteRecent(const char *nm, int n = -1) {
				DBHeaderLink *	hl = hvars->FindVar(nm);
				if (hl == NULL) return;
				hl->DeleteSVal(n);
			}
			// Check that our DB header values are legal
	virtual bool	CheckValues() const {
				if (strcmp(soft, PancineSoftwareName) != 0)
					return false;
				return DBHeader::CheckValues();
			}
			// Compute hash value for current settings
	virtual int32	Hash() const {
				return ( DBHeader::Hash() ^
					(int32)memhash(lastSubj,
						(char *)&uniq - lastSubj, 31) ^
					Unique()<<1 );
			}
};

// Class to track file migrations
class PFileMigration {
	const char *		dold;		// old directory
	const char *		dnew;		// new directory
	int			doldlen;	// strlen(dold)
	int			dnewlen;	// strlen(dnew)
	PFileMigration *	next;		// next in list
	static bool		GetChange(char *d1, char *d2,
						const char *p1, const char *p2);
public:
				PFileMigration(const char *op, const char *np,
						PFileMigration *nfm = NULL);
				~PFileMigration() {
					delete next;
					strunlink(dold); strunlink(dnew);
				}
				// Add new migration path to end of list
	PFileMigration *	AddMigration(const char *op, const char *np);
				// Promote the given member to this position
	PFileMigration *	Promote(const PFileMigration *fmp);
				// Return pointer to next migration path
	const PFileMigration *	GetNext() const {
					return next;
				}
				// Get old directory path
	const char *		GetOldDir() const {
					return dold;
				}
				// Get new directory path
	const char *		GetNewDir() const {
					return dnew;
				}
				// Match the given path to next in our list
	const PFileMigration *	Match(const char *op) const;
				// Match a path pair
	const PFileMigration *	Match(const char *op, const char *np) const;
				// Migrate path using this object (must match)
	bool			Migrate(char *np, int len, const char *op) const;
};

// The following is a really good candidate for someone's preferences file...
extern PFileMigration *	PFMigrations;		// our file migration list

struct DBChangeList;

// Abstract base class to manage snapshot of interface
class PSnapshot {
public:
	virtual		~PSnapshot() {}
			// Restore saved state
	virtual bool	Restore(const DBChangeList *cl, bool before) = 0;
};

// Change list for Undo and Redo functions
struct DBChangeList {
	char		name[32];	// name of change
	DBRecordList	rdel, radd;	// records to delete and add
	PSnapshot *	snap;		// interface snapshot
	DBChangeList *	prev;		// previous change
			DBChangeList(const DBFieldInfo *fi = NULL,
					const char *nm = NULL,
					PSnapshot *ss = NULL) {
				if (fi != NULL) {
					rdel.Init(fi); radd.Init(fi);
				}
				snap = NULL;
				SetName(nm);
				SetSnapshot(ss);
				prev = NULL;
			}
			~DBChangeList() {
				delete snap;
				delete prev;
			}
			// Previous in this named change
	const DBChangeList *
			More() const {
				if (prev == NULL) return NULL;
				if (prev->name[0]) return NULL;
				return prev;
			}
			// Name this change
	void		SetName(const char *nm) {
				if (nm == NULL) nm = "";
				strlcpy(name, nm, sizeof(name));
			}
			// Associate "before" snapshot with this change
	void		SetSnapshot(PSnapshot *ss) {
				delete snap; snap = ss;
			}
			// Update record list according to changes
	bool		ReflectChange(DBRecordList *rl, bool rev = false) const;
			// Get record count for this change
	int		GetSize() const {
				int				n = 0;
				const DBChangeList *	cl;
				for (cl = this; cl != NULL; cl = cl->More())
					n += cl->rdel.GetSize() +
						cl->radd.GetSize();
				return n;
			}
			// Restore snapshot of interface
	bool		Restore(bool before) const {
				if (snap == NULL) return true;
				return snap->Restore(this, before);
			}
};

// Database access manager with undo facility
class PDBAccess : public DBAccess {
private:
	DBChangeList *	undoList;	// list of changes we made
	DBChangeList *	pending;	// change we're about to make
	DBChangeList *	redoList;	// list of changes we unmade
	long		nunsynced;	// number of unsynced record changes
public:
			PDBAccess(DBHeader *hdr=NULL,
					iostream *dbstrm=NULL, bool ro=false) :
						DBAccess(hdr, dbstrm, ro) {
				undoList = NULL;
				pending = NULL;
				redoList = NULL;
				nunsynced = 0;
			}
	virtual		~PDBAccess() {
				ClearChanges();
			}
			// Get name of subclass
	virtual const char *
			ClassName() const {
				return "PDBAccess";
			}
			// Apply pending change and name new one
	bool		Checkpoint(const char *name = NULL,
					PSnapshot *snap = NULL);
			// Undo last change and push onto redo stack
	const DBChangeList *
			Undo();
			// Redo change we undid
	const DBChangeList *
			Redo();
			// Get name of last named change
	const char *	UndoName() const {
				if (pending != NULL && pending->name[0])
					return pending->name;
				const DBChangeList *	cl;
				for (cl = undoList; cl != NULL; cl = cl->prev)
					if (cl->name[0]) return cl->name;
				return NULL;
			}
			// Get name of last undo we can redo
	const char *	RedoName() const {
				const DBChangeList *	cl;
				for (cl = redoList; cl != NULL; cl = cl->prev)
					if (cl->name[0]) return cl->name;
				return NULL;
			}
			// Retrieve our undo change list
	const DBChangeList *
			UndoList() const {
				return undoList;
			}
			// Retrieve our redo change list
	const DBChangeList *
			RedoList() const {
				return redoList;
			}
			// Clear undo and redo lists (also syncs for safety)
	bool		ClearChanges() {
				bool	ok = Sync();
				delete undoList; undoList = NULL;
				delete redoList; redoList = NULL;
				return ok;
			}
			// Execute callback over all database records
	virtual int     ForEachBlock(int (*f)(const DBRecordList &, void *),
					void *cd = NULL,
					const DBFieldSort *ord = NULL) {
				if (!Checkpoint()) return 0;
				return DBAccess::ForEachBlock(f, cd, ord);
			}
			// Pretend to add records to database
	virtual int	AddRecords(const DBRecord *rlist, int n);
			// Pretend to delete records from database
	virtual int	DeleteRecords(const DBRecord *rlist, int n);
			// Get n records starting at rn, return #found
	virtual int	GetRecords(DBRecord *rlist, int n, int rn = 0) {
				if (!Checkpoint()) return 0;
				return DBAccess::GetRecords(rlist, n, rn);
			}
			// Find n records matching query starting from *rnp
	virtual int	FindRecords(const DBQuery *dq, DBRecord *rlist, int n,
					int *rnp = NULL) {
				if (!Checkpoint()) return 0;
				return DBAccess::FindRecords(dq, rlist, n, rnp);
			}
			// Return the total records in database
	virtual int	TotalRecords() {
				int	tot = DBAccess::TotalRecords();
				if ((tot < 0) | (pending == NULL)) return tot;
				return tot + pending->radd.GetSize() -
						pending->rdel.GetSize();
			}
			// Write all database records to a stream
	virtual void	Write(ostream *os, const char tabch = '\t') {
				Checkpoint(); DBAccess::Write(os, tabch);
			}
			// Get record count since last sync
	long		GetNUnsynced() const {
				long	n = nunsynced;
				if (pending != NULL) n += pending->GetSize();
				return n;
			}
			// Synchronize in-memory data to disk
	virtual bool	Sync() {
				bool	ok = (Checkpoint() && DBAccess::Sync());
				nunsynced = 0;
				return ok;
			}
			// Close database object, return true if OK
	virtual bool	Close() {
				bool	ok = ClearChanges();
				return ok & DBAccess::Close();
			}
};

typedef char		PDate[20];		// "YYYY:MM:DD HH:MM:SS"

// Set up query to match the given record
bool		PCsetQuery(DBQuery *dq, const DBRecord &dr, bool fileOnly);

// Translate field id for specific record or record list
#define PDBfieldID(r,fid)	((r).GetFieldInfo()->XlateID(&PDBFInfo,fid))

// Index operator for translating and retrieving a Pancine field
inline const DBField &
PDBgetField(const DBRecord &r, PDBFieldID fid)
{
	return r[PDBfieldID(r,fid)];
}

// Set field value (macro to avoid typing issues)
#define PDBsetField(rp,id,v)	(rp)->SetField(PDBfieldID(*(rp),id), v)

inline int
PDBsetText(DBRecord *rp, PDBFieldID fid, const char *txt, int nl = '\n')
{
	return rp->SetText(PDBfieldID(*rp,fid), txt, nl);
}

// Set field array
#define PDBarrayField(rp,id,a,n) (rp)->SetField(PDBfieldID(*(rp),id), a, n)

// Clear field
#define PDBclearField(rp,id)	(rp)->ClearField(PDBfieldID(*(rp),id))

// Pancine cached image handler -- READ WARNING for CacheObject in "cache.h"
class PCacheImage : CacheObject {
	ImgReader *	irdr;		// image reader
	bool		rdrMine;	// do I own the reader?
	ImgStruct	isrc;		// image loading structure
	ImgColorSpace	myCS;		// local color space holder
	virtual bool	RestoreMemory();
	virtual void	FreeMemory();
public:
	DBRecord	ircd;		// image database record
			PCacheImage() {
				irdr = NULL;
			}
			PCacheImage(const char *fpath) {
				irdr = NULL;
				SetImage(fpath, false);
			}
			PCacheImage(const DBRecord &pr) {
				irdr = NULL;
				SetImage(pr, false);
			}
			PCacheImage(ImgReader *ir) {
				irdr = NULL;
				SetImage(ir);
			}
	virtual		~PCacheImage() {
				CloseImage();
			}
			// Get name of subclass
	virtual const char *
			ClassName() const {
				return "PCacheImage";
			}
			// Someone wants to keep holder
	int		HolderRetain() {
				return CacheObject::HolderRetain();
			}
			// Replacement for delete operator
	virtual void	HolderRelease() {
				CacheObject::HolderRelease();
			}
			// How many people are retaining us?
	int		HolderCount() const {
				return CacheObject::HolderCount();
			}
	void		CloseImage();
	bool		SetImage(const char *fpath, bool quiet = false);
	bool		SetImage(const DBRecord &pr, bool quiet = false);
	bool		SetImage(ImgReader *ir);
	bool		Ready() const {
				return (irdr != NULL);
			}
	ImgReader *	GetReader() const {
				return irdr;
			}
			// These calls work in concert with PfreeImage()
	bool		GetSubimage(ImgStruct *ims, const ImgRect *r);
	bool		GetImage(ImgStruct *ims) {
				return GetSubimage(ims, NULL);
			}
	bool		GetCropped(ImgStruct *ims) {
				int32	rect[4];
				if (PDBgetField(ircd,PDBFcrop).Get(rect,4) != 4)
					return GetImage(ims);
				ImgRect	rct;
				if (irdr->pixAspect < 0.995f) {
					rct.xleft = irdr->pixAspect*rect[0] + .5f;
					rct.xright = irdr->pixAspect*rect[2] + .5f;
				} else {
					rct.xleft = rect[0]; rct.xright = rect[2];
				}
				if (irdr->pixAspect > 1.005f) {
					rct.ytop = rect[1]/irdr->pixAspect + .5f;
					rct.ybottom = rect[3]/irdr->pixAspect + .5f;
				} else {
					rct.ytop = rect[1]; rct.ybottom = rect[3];
				}
				return GetSubimage(ims, &rct);
			}
};

// Struct to facilitate image reorientation (rotations are mod 90 degrees)
struct POrient {
	int32		xres, yres;		// size in original orientation 
	bool		hflip, vflip;		// are horiz. and/or vert. flipped?
	bool		xyswap;			// x and y coords then swapped?
			POrient(const DBRecord *ir = NULL) {
				SetRecord(ir);
			}
			// Set transform according to record values
	void		SetRecord(const DBRecord *ir);
			// Set original image size
	void		SetOrigSize(int xr, int yr) {
				xres = xr; yres = yr;
			}
			// Rotate final image specified degrees clockwise
	void		RotateCW(int cwdeg = 90) {
				while (cwdeg < 0) cwdeg += 360;
				while (cwdeg >= 360) cwdeg -= 360;
				if ((cwdeg <= 45) | (cwdeg >= 305)) return;
				if (cwdeg <= 135) {		// 90 deg CW
					VFlip(); xyswap = !xyswap;
				} else if (cwdeg <= 225) {	// 180 deg
					hflip = !hflip; vflip = !vflip;
				} else {			// 90 deg CCW
					HFlip(); xyswap = !xyswap;
				}
			}
			// Flip final image horizontally
	void		HFlip() {
				if (xyswap) vflip = !vflip;
				else hflip = !hflip;
			}
			// Flip final image vertically
	void		VFlip() {
				if (xyswap) hflip = !hflip;
				else vflip = !vflip;
			}
			// Get final image width
	int32		GetWidth() const {
				if (xyswap) return yres;
				return xres;
			}
			// Get final image height
	int32		GetHeight() const {
				if (xyswap) return xres;
				return yres;
			}
			// Is this mapping a no-op?
	bool		IsNoop() const {
				return !(hflip | vflip | xyswap);
			}
			// Map original to corrected pixel position
	void		MapToCorrect(int *xcp, int *ycp, int xo, int yo) const {
				if (hflip) xo = xres-1 - xo;
				if (vflip) yo = yres-1 - yo;
				if (xyswap) { *xcp = yo; *ycp = xo; }
				else { *xcp = xo; *ycp = yo; }
			}
			// Map corrected to original pixel position
	void		MapToOrig(int *xop, int *yop,
					int xc, int yc) const {
				if (xyswap) { *xop = yc; *yop = xc; }
				else { *xop = xc; *yop = yc; }
				if (hflip) *xop = xres-1 - *xop;
				if (vflip) *yop = yres-1 - *yop;
			}
			// Get inverse mapping
	void		GetInverse(POrient *inv) const {
				if (inv == NULL) return;
				if ((inv->xyswap = xyswap)) {
					inv->xres = yres; inv->yres = xres;
					inv->hflip = vflip; inv->vflip = hflip;
				} else {
					inv->xres = xres; inv->yres = yres;
					inv->hflip = hflip; inv->vflip = vflip;
				}
				return;
			}
			// Reorient an image (can copy in place)
	bool		OrientImage(ImgStruct *dst, const ImgStruct *src = NULL) const;
			// Compare two orient objects for equality
	bool		operator==(const POrient &that) const {
				if (xres != that.xres) return false;
				if (yres != that.yres) return false;
				if (hflip != that.hflip) return false;
				if (vflip != that.vflip) return false;
				return (xyswap == that.xyswap);
			}
};

inline bool
operator!=(const POrient &por1, const POrient &por2)
{
	return !(por1 == por2);
}

// Return values for PlocateFile()
enum PLocateFileRes {
	PLFRfail = -1,			// file not found
	PLFRchanged = 0,		// file contents have changed
	PLFRok = 1,			// file unmoved and unchanged
	PLFRputback,			// file found in previous location
	PLFRmigrated,			// file migrated to another folder
	PLFRrenamed			// file was renamed
};

// Locate Pancine image file and confirm size and hash value
extern PLocateFileRes	PlocateFile(char *pnbuf, int blen, const DBRecord &dr);

// Same as above but call query function to update location with user's help
extern PLocateFileRes	PlocateFile(char *pnbuf, int blen, DBRecord *dr,
				bool (*queryf)(char *, int, DBRecord *));

// Get current file path without any checks at all (return filename pointer)
extern char *		PcurrentFile(char *pnbuf, int blen, const DBRecord &dr);

// Get hash value and length for file
extern bool		PhashFile(int32 *hval, int32 *flen, const char *fname);

// Open an image from the given Pancine image record
// If quiet is set to true, PopenImageD won't report any errors
extern ImgReader *	PopenImageD(const DBRecord &pr, bool quiet = false);

// Fill record fields based on information from file
extern void		PgrokFile(DBRecord *dr, const char *fpath);

// Fill record fields based on info from open image
extern void		PgrokImage(DBRecord *dr, ImgReader *ir);

// Set image information based on database record and return pixel aspect
extern float		PsetInfo(ImgInfo *iip, const DBRecord &pr);

// Check to see if two image orientations match
extern bool		PmatchOrient(const DBRecord &dr1, const DBRecord &dr2);

// Convert "YYYY:MM:DD HH:MM:SS" to relative seconds
extern unsigned long	PsecsFromDate(const PDate dt);

// Reverse conversion
extern char *		PdateFromSecs(PDate dt, unsigned long sec);

// Custom field formatting functions
extern char *		PDBformatFolder(char *, int, const DBField &);
extern char *		PDBformatDate(char *, int, const DBField &);
extern char *		PDBformatExposure(char *, int, const DBField &);
extern char *		PDBformatAperture(char *, int, const DBField &);
extern char *		PDBformatFlash(char *, int, const DBField &);
extern char *		PDBformatWhiteBal(char *, int, const DBField &);
extern char *		PDBformatSource(char *, int, const DBField &);

#endif	// ! _PANCINE_H_
