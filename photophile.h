/*
 *  photophile.h
 *  panlib
 *
 *  Include after "pancine.h"
 *
 *  Classes and declarations for Pancine Photosphere browser
 *
 *  Created by gward on Thu Jun 14 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#ifndef _PHOTOPHILE_H_
#define	_PHOTOPHILE_H_

#ifndef _PANCINE_H_
#include "pancine.h"
#endif

// Class for linked events in timeline
struct PTimeLink {
private:
	PTimeLink *	next;			// next event
	PTimeLink *	prev;			// previous event
	const char *	subject;		// subject of this event
	PDate		begDate, endDate;	// date strings
	unsigned long	begTime, endTime;	// beginning and ending times
	bool		Set(const char *subj, const PDate dt);
	bool		InRange(unsigned long eTime) const {
				return (eTime >= begTime - slopTime &&
						eTime <= endTime + slopTime);
			}
	const PTimeLink *
			MatchEvent(const char *subj, unsigned long eTime) const;
	PTimeLink *	InsertEvent(const char *subj, unsigned long eTime, const PDate dt);
	bool		MergeNext();
	bool		MergePrevious() {
				if (prev == NULL) return false;
				return prev->MergeNext();
			}
	void		TimeOut() {
				if (prev != NULL) prev->next = next;
				if (next != NULL) next->prev = prev;
				prev = next = NULL;
			}
public:
	static unsigned long
			slopTime;		// seconds slop time
			PTimeLink() {
				subject = NULL;
				begTime = endTime = 0L;
				prev = next = NULL;
			}
			PTimeLink(const char *subj, const PDate dt) {
				subject = NULL;
				prev = next = NULL;
				Set(subj, dt);
			}
			PTimeLink(const PTimeLink &orig) {
				subject = NULL;
				prev = next = NULL;
				*this = orig;
			}
			~PTimeLink() {
				KillTime();
			}
			// Find event matching the one given
	const PTimeLink *
			MatchEvent(const char *subj, const PDate dt) const;
			// Insert event in our timeline
	PTimeLink *	InsertEvent(const char *subj, const PDate dt);
			// Add to list with what's happening on a given date
	void		GetHaps(PStringSet *ss, const PDate dt) const;
			// Get event subject
	const char *	GetSubject() const {
				return subject;
			}
			// Get starting date string
	const char *	GetBeginDate() const {
				if (!begTime) return NULL;
				return begDate;
			}
			// Get final date string
	const char *	GetEndDate() const {
				if (!endTime) return NULL;
				return endDate;
			}
			// Return pointer to previous event
	const PTimeLink *
			PrevEvent() const {
				return prev;
			}
			// Return pointer to next event
	const PTimeLink *
			NextEvent() const {
				return next;
			}
			// Return pointer to first event
	const PTimeLink *
			FirstEvent() const {
				if (prev != NULL) return prev->FirstEvent();
				return (begTime ? this : (PTimeLink *)NULL);
			}
			// Return pointer to last event
	const PTimeLink *
			LastEvent() const {
				if (next != NULL) return next->LastEvent();
				return (begTime ? this : (PTimeLink *)NULL);
			}
			// Delete all events from our timeline
	void		KillTime();
			// Copy operator
	PTimeLink &	operator=(const PTimeLink &src);
};

// Class for proxy access to doubly-linked event list
class PTimeline {
private:
	PTimeLink *	timel;			// pointer to timeline
public:
			PTimeline() {
				timel = new PTimeLink;
			}
			PTimeline(const char *subj, const PDate dt) {
				timel = new PTimeLink(subj, dt);
			}
			PTimeline(const PTimeline &orig) {
				timel = new PTimeLink(*orig.timel);
			}
			~PTimeline() {
				delete timel;
			}
			// Return the current event
	const PTimeLink *
			CurrentEvent() const {
				return timel;
			}
			// Find matching event in our timeline
	const PTimeLink *
			MatchEvent(const char *subj, const PDate dt) const {
				return timel->MatchEvent(subj, dt);
			}
			// Add to list with what's happening on a given date
	void		GetHaps(PStringSet *ss, const PDate dt) const {
				return timel->GetHaps(ss, dt);
			}
			// Insert event in our timeline
	PTimeLink *	InsertEvent(const char *subj, const PDate dt) {
				PTimeLink *	newtl = timel->InsertEvent(subj, dt);
				if (newtl == NULL) return NULL;
				return timel = newtl;
			}
			// Return pointer to first event
	const PTimeLink *
			FirstEvent() const {
				return timel->FirstEvent();
			}
			// Return pointer to last event
	const PTimeLink *
			LastEvent() const {
				return timel->LastEvent();
			}
			// Clear  our timeline
	void		KillTime() {
				timel->KillTime();
			}
	PTimeline &	operator=(const PTimeline &src) {
				if (&src == this) return *this;
				*timel = *src.timel;
				return *this;
			}
};

#define PC_SOURCECHAR	'*'			// non-directory start char
#define PC_DBSOURCE	"*DB"			// records from database
#define PC_EXTSOURCE	"*EX"			// records from external source

// Holder for context state
struct PContextState {
protected:
	mutable bool	uniq;			// only one record per file?
	DBRecordList	rcommon;		// DB overlap with rdisp
	DBRecordList	rdisp_link;		// DB links back to rdisp
	DBFieldInfo *	rdisp_info;		// allocated info object
public:
	PStringSet	subjSet;		// subject set
	PTimeline	timeline;		// subjects sorted into events
	PStringSet	albmSet;		// album set
	PStringSet	keywSet;		// keyword set
	PStringSet	ownrSet;		// owner set
	char		rsrc[1024];		// current record source/directory
	DBRecordList	rdisp;			// current display records
	ABitMap		rsel;			// current selection
	char		album[256];		// selected album
			PContextState() {
				rsrc[0] = '\0';
				uniq = false;
				rdisp_info = NULL;
				album[0] = '\0';
			}
			~PContextState() {
				delete rdisp_info;
			}
			// Return true if images are unique
	bool		Unique() const {
				return uniq;
			}
			// Select all display records
	void		SelectAll() {
				rsel.NewBitMap(rdisp.GetSize(), true);
			}
			// Deselect all
	void		SelectNone() {
				rsel.NewBitMap(0);
			}
			// Invert record selection
	void		SelectInverse() {
				if (!rsel.Length()) SelectAll();
				else rsel.Invert();
			}
			// Select (or deselect) the specified record
	void		Select(int i, bool switchon = true) {
				if (switchon & !rsel.Length())
					rsel.NewBitMap(rdisp.GetSize());
				if ((i >= 0) & ((uint32)i < rsel.Length()))
					rsel.Set(i, switchon);
			}
			// Is the given record selected?
	bool		IsSelected(int i) const {
				return rsel.Check(i);
			}
			// Get the next selected record
	bool		NextSelected(int *ip) const {
				if (ip == NULL) return false;
				if (*ip < 0) *ip = 0;
				uint32	pos = rsel.Find(*ip);
				if (pos == ABMend) {
					*ip = rdisp.GetSize();
					return false;
				}
				*ip = pos;
				return true;
			}
			// How many records are selected?
	int		TotalSelected() const {
				return (int)rsel.SumTotal();
			}
			// Get corresponding original DB record
	const DBRecord *
			GetDBOrig(int i) const;
			// Get common record bitmap
	void		GetCommon(ABitMap *cbm) const;
			// How many records in common with DB?
	int		TotalCommon() const {
				return rcommon.GetSize();
			}
			// Find matching display record
	int		FindI(const DBRecord &dr) {
				DBQuery	query;
				if (!PCsetQuery(&query, dr, uniq)) return -1;
				return rdisp.FindFirst(&query);
			}
			// Find linked display record
	int		FindL(const DBRecord &dr) const {
				int	i = rdisp.GetSize();
				while (i--) if (rdisp.Get(i).IsLinked(&dr)) return i;
				return -1;
			}
			// Sort display records, preserving selection
	void		SortDisplay();
			// Add a list of records to our catalog information
	void		Catalog(const DBRecord *rlist, int nr);
			// Catalog a single record
	void		Catalog(const DBRecord &rec) {
				Catalog(&rec, 1);
			}
			// Clear context records
	void		ClearRecords();
			// Clear using the given field info
	void		ClearRecords(const DBFieldInfo *fi, bool mkCopy = false) {
				ClearRecords();
				if (mkCopy & (fi != NULL))
					fi = rdisp_info = new DBFieldInfo(*fi);
				rdisp.Init(fi);
			}
			// Clear context state
	void		ClearSets();
			// Clear everything
	void		ClearState() {
				ClearRecords();
				ClearSets();
			}
			// Utility to transfer record fields
	static bool	TransferFields(DBRecord *rdest, const DBRecord &rsrc);
			// Repair broken links in rdisp_link
	int		FixLinks(bool force = false);
			// Reflect the given change in our context
	void		ReflectChange(const DBChangeList *cl, bool rev = false);
			// Memberwise copy operator (maintains links)
	PContextState &	operator=(const PContextState &src);
};

// Snapshot of context state for undo command
class PCSnapshot : public PSnapshot {
	PContextState *	pcs;			// pointer to context state
	PContextState	s;			// saved context state
	bool		stateSaved;		// state has been set
public:
			PCSnapshot(PContextState *cont = NULL, bool sav = true) {
				SetState(cont, sav);
			}
			// Set the state to the given context
	void		SetState(PContextState *cont, bool sav = true) {
				s.ClearState(); stateSaved = false;
				if ((pcs = cont) != NULL && sav) {
					s = *cont; stateSaved = true;
				}
			}
			// Clear saved state, but not pointer to context
	void		ClearState() {
				s.ClearState(); stateSaved = false;
			}
			// Return saved context state, or NULL if unset
	const PContextState *
			GetState() const {
				if (stateSaved) return &s;
				return NULL;
			}
			// Restore the context state
	virtual bool	Restore(const DBChangeList *cl, bool before) {
				if ((cl == NULL) | (pcs == NULL)) return false;
				if (stateSaved) {
					*pcs = s;
					if (!before) pcs->ReflectChange(cl);
				} else pcs->ReflectChange(cl, before);
				return true;
			}
};	

#ifndef PC_SYNCINTVL
#define PC_SYNCINTVL	32			// record updates before sync
#endif
#ifndef PC_SYNCOPS
#define PC_SYNCOPS	8			// operations before sync
#endif

// Shadow database header class
class SDBHeader : public DBHeader {
public:
			// Set all values -- be sure to clear first!
	virtual void	SetDefaults() {
				DBHeader::SetDefaults();
				strcpy(soft, PancineSoftwareName);
				*dynamic_cast<DBFieldInfo *>(this) =
					*dynamic_cast<const DBFieldInfo *>(&PDBFInfo);
				details[PDBFsubject].flags &= ~DBFFrequired;
			}
};

// Struct to hold Pancine Photosphere browser context
struct PContext : PContextState {
protected:
	char		lock_path[1024];	// DB lock file name
	char		slock_path[1024];	// shadow DB lock file name
	fstream		dbfs;			// DB stream
	fstream		sdbfs;			// shadow DB stream
	SDBHeader	shd;			// shadow DB header
	DBAccess	sdb;			// shadow DB
	static bool	WriteLock(char *lpath, const char *dbpath);
public:
	char		dbname[64];		// Pancine database name
	PDBHeader	phd;			// DB header
	PDBAccess	pdb;			// DB accessor
	DBSearch	dbs;			// current DB search
	int		sync_intvl;		// record changes before sync
	int		sync_nops;		// operations before sync
	int		sync_opcnt;		// operation countdown
			PContext() {
				slock_path[0] = '\0';
				lock_path[0] = '\0';
				dbname[0] = '\0';
				sync_intvl = PC_SYNCINTVL;
				sync_opcnt = sync_nops = PC_SYNCOPS;
			}
			~PContext() {
				WriteLock(slock_path, NULL);
				WriteLock(lock_path, NULL);
			}
			// Set shadow database file
	bool		ShadowDB(const char *dbpath);
			// Open new database and initialize
	bool		OpenDB(const char *dbpath, bool ro = false);
			// Is database hooked up and ready to go?
	bool		Ready() const {
				return (dbname[0] && pdb.Ready());
			}
			// Is database not writeable?
	bool		ReadOnly() const {
				if (!Ready()) return true;
				return pdb.ReadOnly();
			}
			// Are DB images unique?
	bool		Unique() const {
				if (Ready()) uniq = phd.Unique();
				return uniq;
			}
			// Set DB header for unique images
	void		SetUnique(bool val) {
				uniq = val;
				if (ReadOnly()) return;
				pdb.ClearChanges();
				phd.SetUnique(uniq);
			}
			// Clear display records
	void		Clear(const DBFieldInfo *fi = NULL) {
				if (fi != NULL) ClearRecords(fi, true);
				else if (Ready()) ClearRecords(&phd, false);
				else ClearRecords();
			}
			// Return editable field information, or NULL if none
	DBFieldInfo *	EditableFieldInfo() {
				if (rdisp_info != NULL) return rdisp_info;
				if (rdisp.GetFieldInfo() == &phd) return &phd;
				return NULL;
			}
			// Use dbs to search database for display records
	int		Research() {
				if (!Ready()) return 0;
				Clear();
				dbs.FindAll(&rdisp);
				strcpy(rsrc, PC_DBSOURCE);
				Reconcile();
				return rdisp.GetSize();
			}
			// Compute rcommon and transfer field values to rdisp
	int		Reconcile();
			// Reconcile a single record with DB
	bool		Reconcile(DBRecord *drp);
			// Recatalog the entire DB
	void		Recatalog();
			// Undo previous change to DB
	bool		Undo() {
				if (ReadOnly()) return false;
				const DBChangeList *	cl = pdb.Undo();
				if (cl == NULL) return false;
				if (!cl->rdel.GetSize())	// undoing Add?
					AddShadowList(cl->radd);
				return true;
			}
			// Redo last undone change to DB
	bool		Redo() {
				if (ReadOnly()) return false;
				return (pdb.Redo() != NULL);
			}
			// Clear undo/redo lists
	bool		ClearUndo() {
				if (ReadOnly()) return false;
				sync_opcnt = sync_nops;
				return pdb.ClearChanges();
			}
			// Checkpoint the database (and sync if time)
	bool		Checkpoint(const char *name = NULL,
					PSnapshot *snap = NULL);
			// Synchronize the database
	bool		SyncDB() {
				sync_opcnt = sync_nops;
				if (ReadOnly()) return true;
				return pdb.Sync();
			}
			// Add current (selected) records to database
	int		AddRecords();
			// Delete current (selected) records from database
	int		DeleteRecords(bool fromDB = true, bool fromSDB = true);
			// Write change in display records to database
	int		RecordChange(const char *nm = NULL,
					ABitMap *bm = NULL,
					PCSnapshot *ss = NULL);
			// Add shadow records
	int		AddShadows(const DBRecord *rlist, int n);
	int		AddShadowList(const DBRecordList &rl) {
				return AddShadows(rl.GetArray(), rl.GetSize());
			}
			// Add record to shadow catalog
	bool		AddShadow(const DBRecord &dr) {
				return (bool)AddShadows(&dr, 1);
			}
			// Delete shadow records
	int		DeleteShadows(const DBRecord *rlist, int n);
	int		DeleteShadowList(const DBRecordList &rl) {
				return DeleteShadows(rl.GetArray(), rl.GetSize());
			}
			// Delete record from shadow catalog
	bool		DeleteShadow(const DBRecord &dr) {
				return DeleteShadows(&dr, 1);
			}
			// Close database and free resources
	bool		CloseDB();
};

#endif	// ! _PHOTOPHILE_H_
