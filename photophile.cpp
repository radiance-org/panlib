/*
 *  photophile.cpp
 *  panlib
 *
 *  Pancine Photosphere browser support routines
 *
 *  Created by gward on Thu Jun 14 2001.
 *  Copyright (c) 2003 Anyhere Software. All rights reserved.
 *
 */

#include "pancine.h"
#include "pstrings.h"
#include "photophile.h"
#include <fcntl.h>
#include <signal.h>
#if defined(_WIN32) || defined(_WIN64)
#include <process.h>				// contains getpid() declaration
#endif

#ifndef P_MAXARR
#define P_MAXARR		128		// largest array we handle
#endif
#ifndef P_MAXSTATE
#define P_MAXSTATE		10240		// max. records to save in state
#endif
#ifndef PTL_SLOP
#define	PTL_SLOP		(2L*24L*3600L)	// default slop time (seconds)
#endif

unsigned long	PTimeLink::slopTime = PTL_SLOP;	// seconds slop time

/*
 *  The PTimeLink class keeps events sorted by date, merging
 *  events with the same subject that are within slopTime of
 *  one another.  If a different subject comes up in the middle
 *  of another event, that event is split so as to maintain
 *  order in the timeline.  Other than that, events may be
 *  extended to overlap each other, just so long as one event
 *  does not end up inside another.  The list is kept sorted
 *  by event start time.
 */

// Set the current event to a single time
bool
PTimeLink::Set(const char *subj, const PDate dt)
{
	if ((subj == NULL) | (dt == NULL))
		return false;
	unsigned long	sec = PsecsFromDate(dt);
	if (!sec)
		return false;
	begTime = endTime = sec;
	if (subject != NULL)
		strunlink(subject);
	subject = strlink(subj);
	strcpy(begDate, dt);
	strcpy(endDate, dt);
	return true;
}

// Search timeline for a match to the given subject+time (public interface)
const PTimeLink *
PTimeLink::MatchEvent(const char *subj, const PDate dt) const
{
	if ((subj == NULL) | (dt == NULL))
		return NULL;
	if (prev != NULL)
		return prev->MatchEvent(subj, dt);
	return MatchEvent(subj, PsecsFromDate(dt));
}

// Search timeline for a match (private interface)
const PTimeLink *
PTimeLink::MatchEvent(const char *subj, unsigned long eTime) const
{
	if (InRange(eTime) && (subj == NULL || !istrcmp(subj, subject)))
		return this;
	if ((eTime <= begTime) | (next == NULL))
		return NULL;
	return next->MatchEvent(subj, eTime);
}

// Add to list with what's happening on a given date
void
PTimeLink::GetHaps(PStringSet *ss, const PDate dt) const
{
	if ((ss == NULL) | (dt == NULL))
		return;
	if (prev != NULL) {
		prev->GetHaps(ss, dt);
		return;
	}
	int	sec = PsecsFromDate(dt);
	if (!sec)
		return;
	const PTimeLink *	ep = MatchEvent(NULL, sec);
	while (ep != NULL && ep->begTime <= sec + ep->slopTime) {
		ss->Insert(ep->subject);
		ep = ep->next;
	}
}

// Insert event in our timeline (public interface)
PTimeLink *
PTimeLink::InsertEvent(const char *subj, const PDate dt)
{
	if ((subj == NULL) | (dt == NULL))
		return NULL;
	if (prev != NULL)			// back up to start
		return prev->InsertEvent(subj, dt);
	if (!begTime) {				// initial setting
		if (Set(subj, dt))
			return this;
		return NULL;
	}
	unsigned long	sec = PsecsFromDate(dt);
	if (!sec)
		return NULL;
	return InsertEvent(subj, sec, dt);
}

// Insert event into our timeline (private interface)
PTimeLink *
PTimeLink::InsertEvent(const char *subj, unsigned long eTime, const PDate dt)
{
						// check for merge
	if (InRange(eTime) && !istrcmp(subj, subject)) {
		if (eTime < begTime) {		// extend backwards
			begTime = eTime;
			strcpy(begDate, dt);
			MergePrevious();
		} else if (eTime > endTime) {	// extend forwards
			endTime = eTime;
			strcpy(endDate, dt);
			MergeNext();
		}
		return this;
	}
	if (eTime <= begTime) {			// need to insert before
		PTimeLink *	enew = new PTimeLink(subj, dt);
		if ((enew->prev = prev) != NULL)
			prev->next = enew;
		enew->next = this;
		return prev = enew;
	}
	if (next == NULL) {			// new end of the line
		next = new PTimeLink(subj, dt);
		next->prev = this;
		return next;
	}
						// still looking...
	return next->InsertEvent(subj, eTime, dt);
}

// Merge next event if appropriate
bool
PTimeLink::MergeNext()
{
	PTimeLink *	nr = next;
	if (nr == NULL)
		return false;
	if (nr->subject != subject)
		return false;			// different subject
	DASSERT(begTime != 0);
	DASSERT(nr->begTime != 0);
	if (endTime + slopTime < nr->begTime)
		return false;			// no overlap

	if (endTime < nr->endTime) {		// else merge
		strcpy(endDate, nr->endDate);
		endTime = nr->endTime;
	}
	if (begTime > nr->begTime) {		// should never happen
		strcpy(begDate, nr->begDate);
		begTime = nr->begTime;
	}
	next = nr->next;			// delete next
	if (next != NULL)
		next->prev = this;
	nr->next = nr->prev = NULL;
	delete nr;
	return true;				// we did it
}

// Delete all events from our timeline
void
PTimeLink::KillTime()
{
	if (prev != NULL) {
		prev->next = NULL;
		delete prev;
		prev = NULL;
	}
	if (next != NULL) {
		next->prev = NULL;
		delete next;
		next = NULL;
	}
	strunlink(subject);
	subject = NULL;
	begTime = endTime = 0L;
}

// Copy operator
PTimeLink &
PTimeLink::operator=(const PTimeLink &src)
{
	if (this == &src)
		return *this;
	KillTime();
	subject = strlink(src.subject);
	memcpy(this, &src, sizeof(PTimeLink));
	if (src.next != NULL) {
		PTimeLink	chain;
		memcpy(&chain, src.next, sizeof(PTimeLink));
		chain.prev = NULL;
		next = new PTimeLink(chain);
		next->prev = this;
		memset(&chain, 0, sizeof(PTimeLink));
	}
	if (src.prev != NULL) {
		PTimeLink	chain;
		memcpy(&chain, src.prev, sizeof(PTimeLink));
		chain.next = NULL;
		prev = new PTimeLink(chain);
		prev->next = this;
		memset(&chain, 0, sizeof(PTimeLink));
	}
	return *this;
}

// Get corresponding original DB record
const DBRecord *
PContextState::GetDBOrig(int i) const
{
	if ((i < 0) | (i >= rdisp.GetSize()))
		return NULL;
	const DBRecord &	dr = rdisp.Get(i);
	if (!dr.IsLinked())
		return NULL;
						// check obvious place first
	if (i < rdisp_link.GetSize() && dr.IsLinked(&rdisp_link.Get(i)))
		return &rcommon.Get(i);
						// try searching
	for (i = rdisp_link.GetSize(); i--; )
		if (dr.IsLinked(&rdisp_link.Get(i)))
			return &rcommon.Get(i);
	return NULL;
}

// Get common record bitmap
void
PContextState::GetCommon(ABitMap *cbm) const
{
	if (cbm == NULL)
		return;
	cbm->NewBitMap(rdisp.GetSize(), true);
	if (!strcmp(rsrc, PC_DBSOURCE))
		return;				// all common
	for (int i = rdisp.GetSize(); i--; )
		if (GetDBOrig(i) == NULL)
			cbm->Reset(i);
}

// Clear context records
void
PContextState::ClearRecords()
{
	SelectNone();
	rdisp.Init(&PDBFInfo);
	rcommon.Init();
	rdisp_link.Init();
	delete rdisp_info; rdisp_info = NULL;
}

// Clear context state
void
PContextState::ClearSets()
{
	subjSet.Clear();
	timeline.KillTime();
	ownrSet.Clear();
	albmSet.Clear();
	keywSet.Clear();
}

// Reflect the given change in our context
void
PContextState::ReflectChange(const DBChangeList *cl, bool rev)
{
	if (cl == NULL)
		return;
	SelectNone();
	cl->ReflectChange(&rdisp, rev);
	cl->ReflectChange(&rcommon, rev);
	FixLinks(true);
	if (rev)
		return;
				// update catalog info.
	do
		Catalog(cl->radd.GetArray(), cl->radd.GetSize());
	while ((cl = cl->More()) != NULL);
}

// Catalog a list of records (note new subjects, collections, etc.)
void
PContextState::Catalog(const DBRecord *rlist, int nr)
{
	if (rlist == NULL)
		return;
	const char *	sl[P_MAXARR];
	int		n;
	for ( ; nr-- > 0; rlist++) {
		if (!rlist->GetNAlloc())
			continue;
		if (PDBgetField(*rlist,PDBFsubject).Get(sl)) {
			subjSet.Insert(sl[0]);
			const char *	dt;
			if (PDBgetField(*rlist,PDBFcapdate).Get(&dt))
				timeline.InsertEvent(sl[0], dt);
		}
		if (PDBgetField(*rlist,PDBFowner).Get(sl))
			ownrSet.Insert(sl[0]);
		n = PDBgetField(*rlist,PDBFalbum).Get(sl, P_MAXARR);
		if (n > 0)
			albmSet.Add(sl, n);
		n = PDBgetField(*rlist,PDBFkeyword).Get(sl, P_MAXARR);
		if (n > 0)
			keywSet.Add(sl, n);
	}
}

// Repair broken links in rdisp_link
int
PContextState::FixLinks(bool force)
{
	int	nfixed = 0;
	int	nbroken = 0;
	int	rc = rcommon.GetSize();

	if (force || rdisp_link.GetSize() != rc) {
		rdisp_link.Init(NULL, nbroken = rc);
		rc = 0;
	}
	while (rc--) {			// try working on broken links only
		if (!rdisp_link.Get(rc).GetNAssigned()) {
			++nbroken;	// no value set!
			continue;
		}
		if (rdisp_link.Get(rc).IsLinked())
			continue;	// assume link is OK
		int	i = rdisp.GetSize();
		while (i--)
			if (rdisp_link.Get(rc) == rdisp.Get(i))
				break;
		if (i < 0) {
			++nbroken;	// missing link
			continue;
		}
		rdisp_link[rc].Link(&rdisp[i]);
		++nfixed;
	}
	if (!nbroken)			// got them all?
		return nfixed;
					// else try inexact search
	if (!force)
		DMESG(DMCwarning, "Resorting to record search in FixLinks");
	for (rc = rcommon.GetSize(); rc--; ) {
		if (rdisp_link.Get(rc).IsLinked())
			continue;
		DBQuery	query;
		query.matchEmpty = true;
		if (!PCsetQuery(&query, rcommon.Get(rc), uniq))
			DMESG(DMCassert, "PCsetQuery failed in FixLinks");
		int	i = rdisp.FindFirst(&query);
		if ((i < 0) & !uniq) {
			PCsetQuery(&query, rcommon.Get(rc), true);
			i = rdisp.FindFirst(&query);
		}
		if (i < 0) {
			DMESG(DMCwarning, "Lost DB link in FixLinks");
			rcommon[rc].Init();
			rdisp_link[rc].Init();
			continue;
		}
					// reestablish link
		rdisp_link[rc].Link(&rdisp[i]);
		++nfixed;
	}
	rcommon.Compact();
	rdisp_link.Compact();
	return nfixed;
}

// Utility to transfer fields from catalog record
bool
PContextState::TransferFields(DBRecord *rdest, const DBRecord &rsrc)
{
	if (rdest == NULL || !rsrc.GetNAssigned())
		return false;
	bool		didSome = false;
	int		sidmap[DB_MAXFIELD];
	rsrc.GetFieldInfo()->MapIDs(sidmap, rdest->GetFieldInfo());
	for (int i = rdest->GetNFields(); i-- > 0; ) {
		const DBField *	srcf = rsrc.GetField(sidmap[i]);
		if (srcf == NULL)
			continue;
		if (rdest->GetField(i) == NULL)
			didSome |= rdest->SetField(i, *srcf);
		else if (rdest->GetFieldDetails(i)->flags & DBFFarray) {
			DBField	fv;	// append certain array fields
			switch (PDBFInfo.XlateID(rdest->GetFieldInfo(),i)) {
			case PDBFlocation:
			case PDBFkeyword:
			case PDBFcomment:
			case PDBFhistory:
				fv = (*rdest)[i];
				fv |= *srcf;
				if (fv.GetNV() > (*rdest)[i].GetNV())
					didSome |= rdest->SetField(i, fv);
			}
		} else			// special case overrides
			switch (PDBFInfo.XlateID(rdest->GetFieldInfo(),i)) {
			case PDBFcapdate:
			case PDBFlatitude:
			case PDBFlongitude:
			case PDBFaltitude:
			case PDBFgmtdate:
			case PDBForient:
			case PDBFflip:
				if (*srcf != (*rdest)[i])
					didSome |= rdest->SetField(i, *srcf);
			}
	}
	return didSome;
}

// Local record pointer comparison function for qsort()
static int
rptr_cmp(const void *rpp1, const void *rpp2)
{
	return (**(DBRecord * const *)rpp1).Compare(**(DBRecord * const *)rpp2);
}

// Sort display records, preserving selection
void
PContextState::SortDisplay()
{
	if (rdisp.GetSize() <= 0 || rdisp.sord.GetSort(0) < 0)
		return;
	const int	nsel = TotalSelected();
	const bool	savSel = (2*nsel <= rdisp.GetSize());
	DBRecordList	savLink(NULL, savSel ? nsel : rdisp.GetSize()-nsel);
	int		n;
						// simple method for short list
	if (savLink.GetSize()*double(rdisp.GetSize()) < 1e4) {
		n = 0;				// remember (de)selected records
		for (uint32 ui = 0; rsel.Find(&ui, savSel); ui++)
			savLink[n++].Link(&rdisp[ui]);
		rdisp.Sort();			// sort display records */
		if (savSel)			// reset selection
			SelectNone();
		else
			SelectAll();
		while (n--)			// (de)select linked records
			Select(FindL(savLink.Get(n)), savSel);
		return;
	}
						// sort long list by pointers
	if (!savLink.LinkAll(&rdisp)) {
		SelectNone();			// punt (loses selection)
		rdisp.Sort();
		return;
	}
	const ABitMap	rbef = rsel;
	DBRecord **	sortp = new DBRecord * [rdisp.GetSize()];
	for (n = rdisp.GetSize(); n--; ) {
		savLink[n].sortord = &rdisp.sord;
		sortp[n] = savLink.Array() + n;
	}
	qsort(sortp, rdisp.GetSize(), sizeof(*sortp), rptr_cmp);
						// link back in order
	for (n = rdisp.GetSize(); n--; )
		rdisp[n].Link(sortp[n]);
						// reshuffle selection, too
	if (savSel) {
		SelectNone();
		for (n = rdisp.GetSize(); n--; )
			if (rbef.Check(sortp[n] - savLink.GetArray()))
				Select(n, true);
	} else {
		SelectAll();
		for (n = rdisp.GetSize(); n--; )
			if (!rbef.Check(sortp[n] - savLink.GetArray()))
				Select(n, false);
	}
	delete [] sortp;			// that's it
}

// Memberwise copy operator (maintains links)
PContextState &
PContextState::operator=(const PContextState &src)
{
	if (this == &src)
		return *this;
	delete rdisp_info; rdisp_info = NULL;
	if (src.rdisp_info != NULL)
		rdisp_info = new DBFieldInfo(*src.rdisp_info);
	subjSet = src.subjSet;
	timeline = src.timeline;
	ownrSet = src.ownrSet;
	albmSet = src.albmSet;
	keywSet = src.keywSet;
	memcpy(rsrc, src.rsrc, sizeof(rsrc));
	memcpy(album, src.album, sizeof(album));
	uniq = src.uniq;
	rdisp = src.rdisp;
	rcommon = src.rcommon;
	rsel = src.rsel;
	int	i, rc;			// duplicate link relationships
	rdisp_link.Init(NULL, rc = rcommon.GetSize());
	while (rc--) {
		i = rdisp.GetSize();
		while (i--)
			if (src.rdisp.Get(i).IsLinked(&src.rdisp_link.Get(rc)))
				break;
		if (i < 0)
			DMESG(DMCassert, "Missing DB link in context copy");
		
		if (!rdisp_link[rc].Link(&rdisp[i]))
			DMESG(DMCassert, "Link failed in context copy");
	}
	return *this;
}

int
PContext::AddShadows(const DBRecord *rlist, int n)
{
	if (!sdb.Ready())
		return 0;
	DeleteShadows(rlist, n);
					// make sure fields match
	DBRecordList	rlistn(sdb.GetHeader(), n);
	for (int i = 0; i < n; i++)
		rlistn[i] = rlist[i];
	return sdb.AddRecordList(rlistn);
}

// Delete matching records from shadow catalog
int
PContext::DeleteShadows(const DBRecord *rlist, int n)
{
	if (!sdb.Ready() || sdb.TotalRecords() <= 0)
		return 0;
	DBSearch	sdbs(&sdb);
	DBQuery		query;
	DBRecordList	matching;
					// set up shadow DB query
	while (n-- > 0) {
					// XXX should just check location/file
		if (!PCsetQuery(&query, rlist[n], true))
			continue;
		if (!sdbs.AddSelector(query.rlo, query.rhi))
			return 0;
	}
	if (sdbs.FindAll(&matching) <= 0)
		return 0;
	return sdb.DeleteRecordList(matching);
}

// Checkpoint the database (and sync if time)
bool
PContext::Checkpoint(const char *name, PSnapshot *snap)
{
	if (ReadOnly()) {
		delete snap;
		return false;
	}
					// time to sync DB?
	if (!--sync_opcnt || (sync_intvl > 0 && pdb.GetNUnsynced() >= sync_intvl))
		if (!SyncDB()) {
			delete snap;
			return false;
		}
	if (snap != NULL) {		// check for too much saved state
		int			nsaved = 0;
		const PContextState *	cs = (*(PCSnapshot *)snap).GetState();
		const DBChangeList *	cl = pdb.UndoList();
		if (cs != NULL)
			nsaved = cs->rdisp.GetSize() + cs->TotalCommon();
		while (cl != NULL) {
			if (cl->snap != NULL &&
					(cs = (*(PCSnapshot *)cl->snap).GetState())
						!= NULL) {
				nsaved += cs->rdisp.GetSize();
				nsaved += cs->TotalCommon();
				if (nsaved > P_MAXSTATE)
					break;
			}
			cl = cl->prev;
		}
		while (cl != NULL) {	// clear the remainder of undo list
			if (cl->snap != NULL)
				(*(PCSnapshot *)cl->snap).ClearState();
			cl = cl->prev;
		}
	}
	return pdb.Checkpoint(name, snap);
}

// Find common records and reconcile field values between DB and display
int
PContext::Reconcile()
{
	rcommon.Init(&phd);		// clear previous common records
	rdisp_link.Init();
	if (!rdisp.GetSize())
		return 0;
	int	i;
	if (!strcmp(rsrc, PC_DBSOURCE)) {
		rcommon = rdisp;	// all records are common!
		for (i = rcommon.GetSize(); i--; )
			rcommon[i].SetReadOnly();
		rdisp_link.LinkAll(&rdisp);
		return rcommon.GetSize();
	}
	dbs.Clear();			// set up DB search for commonalities
	DBSearch	sdbs(&sdb);
	const bool      doQuery = (Ready() && pdb.TotalRecords() > 0);
	const bool	doShadow = (sdb.Ready() && sdb.TotalRecords() > 0);
	const bool	defUniq = (!doQuery | uniq);
	ABitMap		hashMap(primeAbove(rdisp.GetSize()*10));
	DBQuery		query;
	for (i = 0; i < rdisp.GetSize(); i++) {
		DBRecord &	dr = rdisp[i];
		const char *	fn;
		if (!PDBgetField(dr,PDBFfile).Get(&fn)) {
			DMESG(DMCwarning, "Ignoring record w/o file name");
			dr.Init();
			continue;
		}
		if (!PCsetQuery(&query, dr, defUniq)) {
			char	fname[1024];
			int32	hval, nbyt;
			if (PcurrentFile(fname, sizeof(fname), dr) == NULL ||
					!PhashFile(&hval, &nbyt, fname)) {
				DMESGF(DMCwarning, "Cannot hash '%s'", fname);
				dr.Init();
				continue;
			}
			if (!PDBsetField(&dr, PDBFhashval, hval) ||
					!PDBsetField(&dr, PDBFnbytes, nbyt) ||
					!PCsetQuery(&query, dr, defUniq)) {
				DMESGF(DMCparameter, "Bad record '%s'", fn);
				dr.Init();
				continue;
			}
		}
		bool	maybe = true;   // check if we have it already
		long	hsh = query.HashIndex();
		if (hsh >= 0)		// use hash map to accelerate search
			maybe = !hashMap.TestAndSet(hsh % hashMap.Length());
		int	j = (maybe ? rdisp.FindFirst(&query) : i);
		if (j < 0)
			DMESG(DMCassert, "Search botch in Reconcile");
		if (j < i) {
			DMESGF(DMCwarning, "Duplicate image '%s' ignored", fn);
			TransferFields(&dr, rdisp.Get(j));
			rdisp[j].Take(&dr);
			continue;
		}
					// add it to our query list(s)
		if (doQuery && !dbs.AddSelector(query.rlo, query.rhi))
			DMESG(DMCassert, "AddSelector failed in Reconcile");
		if (doShadow) {
			if (defUniq || PCsetQuery(&query, dr, true))
				sdbs.AddSelector(query.rlo, query.rhi);
		}
	}
	hashMap.NewBitMap(0);
	rdisp.Compact();		// deletes duplicate images

	if (!(doQuery | doShadow))	// shall we query DB?
		return 0;
					// perform search and transfer fields
	if (doQuery) {
		i = dbs.FindAll(&rcommon);
		dbs.Clear();
		rdisp_link.Init(NULL, i);
	}
	DBRecordList	rshadow;
	if (doShadow) {
		rshadow.Init(&shd);
		sdbs.FindAll(&rshadow);
	}
	for (i = rdisp.GetSize(); i--; ) {
		DBRecord &	dr = rdisp[i];
		int		j;
					// identify DB record if one
		PCsetQuery(&query, dr, defUniq);
		if ((j = rcommon.FindFirst(&query)) < 0) {
			if (!doShadow)
				continue;
			if (!defUniq && !PCsetQuery(&query, dr, true))
				continue;
					// check shadow database
			if ((j = rshadow.FindFirst(&query)) >= 0) {
				TransferFields(&dr, rshadow.Get(j));
				rshadow[j].Init();
			}
			continue;
		}
					// establish link
		rdisp_link[j].Link(&dr);
					// fill in unassigned fields
		TransferFields(&dr, rcommon.Get(j));
					// eliminate duplicate matches
		while ((j = rcommon.FindNext()) >= 0)
			rcommon[j].Init();
	}
	rcommon.Compact();		// deletes duplicate matches
	rdisp_link.Compact();
	DASSERT(rcommon.GetSize() == rdisp_link.GetSize());
	if (doShadow) {			// delete unneeded shadow matches
		rshadow.Compact();
		sdb.DeleteRecordList(rshadow);
	}
					// rcommon records are read-only
	for (i = rcommon.GetSize(); i--; )
		rcommon[i].SetReadOnly();
	return rcommon.GetSize();
}

// Reconcile a single record with DB (does not affect rcommon)
bool
PContext::Reconcile(DBRecord *drp)
{
	if (drp == NULL || drp->GetFieldInfo() == NULL)
		return false;
	const bool	defUniq = (!Ready() | uniq);
	DBQuery		query;
	DBRecord	dbr;
	if (!PCsetQuery(&query, *drp, defUniq))
		return false;
	if (!Ready() || !pdb.FindRecords(&query, &dbr, 1)) {
		if (!sdb.Ready())
			return false;
		if (!defUniq && !PCsetQuery(&query, *drp, true))
			return false;
					// check shadow DB
		if (!sdb.FindRecords(&query, &dbr, 1))
			return false;
	}
	return TransferFields(drp, dbr);
}

// Write change in display records to database
int
PContext::RecordChange(const char *nm, ABitMap *bm, PCSnapshot *ss)
{
	const int	nrec = rdisp.GetSize();
	if (!nrec) {
		delete ss;
		return 0;
	}
	ABitMap		chMap;
	if (bm == NULL) {
		chMap.NewBitMap(nrec, true);
		bm = &chMap;
	} else if (bm->Find() == ABMend) {
		delete ss;
		return 0;
	}
	const bool	chngDB = !ReadOnly();
	const bool	onlyDB = (chngDB && !strcmp(rsrc, PC_DBSOURCE));
	const bool	chngShadow = sdb.Ready();
	int *		repl = new int [nrec];
	DBRecordList	inShadow(&shd);
	uint32		i;
	for (i = nrec; i--; ) {		// have they actually changed?
		repl[i] = -1;
		if (!bm->Check(i))
			continue;
		const DBRecord *	crp = GetDBOrig(i);
		if (crp == NULL) {
			if (onlyDB)
				DMESG(DMCwarning,
					"Record not found in RecordChange");
			else if (chngShadow)
				inShadow.NewRecord() = rdisp[i];
			bm->Reset(i);	// DB doesn't change
		} else if (chngDB && rdisp.Get(i) != *crp)
			repl[i] = crp - rcommon.GetArray();
		else
			bm->Reset(i);	// no change
	}
	if (chngShadow)			// update shadow records
		AddShadowList(inShadow);
	int	nChange = bm->SumTotal();
	if (!nChange) {
		delete [] repl;
		delete ss;
		return 0;		// no change to DB
	}
					// create checkpoint
	if (!Checkpoint(nm, ss)) {
		delete [] repl;
		DMESG(DMCdata, "Database update error in RecordChange");
		return -1;
	}
					// update records
	DBRecordList	toDel(NULL, nChange);
	DBRecordList	toAdd(NULL, nChange);
	int		j;
	for (i = j = 0; bm->Find(&i); i++, j++) {
		toDel[j].Take(&rcommon[repl[i]]);
		rcommon[repl[i]] = rdisp.Get(i);
		rcommon[repl[i]].SetReadOnly();
		toAdd[j].Link(&rdisp[i]);
	}
	pdb.DeleteRecordList(toDel);
	pdb.AddRecordList(toAdd);
	Catalog(toAdd.GetArray(), toAdd.GetSize());
	delete [] repl;
	return nChange;
}

// Delete selected display records from display (and DB if fromDB is set)
int
PContext::DeleteRecords(bool fromDB, bool fromSDB)
{
	int	nDel = TotalSelected();
	if (!nDel)
		return 0;		// none selected means delete none!
	DBRecordList	delList(NULL, nDel);
	PCSnapshot *	ss = NULL;
	DBRecordList	toDel;
	fromDB &= !ReadOnly();		// identify DB records
	if (fromDB)
		ss = new PCSnapshot(this);
	const bool	onlyDB = (fromDB && !strcmp(rsrc, PC_DBSOURCE));
	nDel = 0;
	for (int i = 0; NextSelected(&i); i++) {
		const DBRecord *	crp = GetDBOrig(i);
		delList[nDel++].Take(&rdisp[i]);
		if (crp == NULL) {
			if (onlyDB)
				DMESG(DMCwarning,
					"Record not found in DeleteRecords");
			continue;
		}
		int	rc = crp - rcommon.GetArray();
		toDel.NewRecord().Take(&rcommon[rc]);
		rdisp_link[rc].Init();
	}
	SelectNone();
	rdisp.Compact();
	rcommon.Compact();
	rdisp_link.Compact();
	fromSDB &= sdb.Ready();
	if (fromSDB)			// delete records from shadow DB
		DeleteShadowList(delList);
	fromDB &= (toDel.GetSize() > 0);
	if (!fromDB) {
		delete ss;
		return 0;		// don't delete from DB
	}
					// create checkpoint
	if (!Checkpoint("Delete", ss)) {
		DMESG(DMCdata, "Database update error in DeleteRecords");
		return -1;
	}
					// delete DB records
	pdb.DeleteRecordList(toDel);
	return toDel.GetSize();
}

// Add current (selected) display records to database
int
PContext::AddRecords()
{
	if (!strcmp(rsrc, PC_DBSOURCE))
		return 0;
	if (ReadOnly())
		return 0;
	const int	nrec = rdisp.GetSize();
	if (!nrec)
		return 0;
	ABitMap		addMap = rsel;
	int		nAdd = addMap.SumTotal();
	uint32		i;
	if (!nAdd)			// none selected means do all
		addMap.NewBitMap(nAdd=nrec, true);
					// make sure we have required fields
	for (i = 0; addMap.Find(&i); i++) {
		const DBRecord &	dr = rdisp.Get(i);
		for (int j = dr.GetNFields(); j--; )
			if (dr.GetFieldDetails(j)->flags & DBFFrequired
					&& dr.GetField(j) == NULL) {
				DMESG(DMCwarning,
					"Incomplete record not added to DB");
				addMap.Reset(i);
				nAdd--;
			}
	}
	if (!nAdd)
		return 0;
					// identify replacement records
	DBRecordList	reprecs;
	int *		repl = new int [nrec];
	for (i = nrec; i--; ) {
		repl[i] = -1;
		if (!addMap.Check(i))
			continue;
		const DBRecord *	crp = GetDBOrig(i);
		if (crp == NULL)
			continue;
					// avoid updating identical records
		if (rdisp.Get(i) == *crp) {
			addMap.Reset(i);
			nAdd--;
			continue;
		}
		repl[i] = crp - rcommon.GetArray();
		reprecs.NewRecord().Link(&rcommon[repl[i]]);
	}
	if (nAdd <= 0) {
		delete [] repl;
		return 0;		// nothing new to add
	}
					// checkpoint database
	if (!Checkpoint("Add", new PCSnapshot(this))) {
		DMESG(DMCdata, "Database update error in AddRecords");
		delete [] repl;
		return -1;
	}
					// first, delete ones we'll replace
	pdb.DeleteRecordList(reprecs);
					// now, add new records
	DBRecordList	addList(NULL, nAdd);
	int		j;
	for (i = j = 0; addMap.Find(&i); i++, j++) {
		addList[j].Link(&rdisp[i]);
		if (repl[i] >= 0) {
			rcommon[repl[i]] = rdisp.Get(i);
			rcommon[repl[i]].SetReadOnly();
			// already linked
		} else {
			rcommon.NewRecord() = rdisp.Get(i);
			rcommon[rcommon.GetSize()-1].SetReadOnly();
			rdisp_link.NewRecord().Link(&rdisp[i]);
		}
	}
	DeleteShadowList(addList);	// remove from shadow DB
					// add to catalog DB
	pdb.AddRecordList(addList);
	Catalog(addList.GetArray(), nAdd);
					// clean up and return number added
	delete [] repl;
	return nAdd;
}

// Recatalog callback
static int
cat_call(const DBRecordList &drl, void *myObject)
{
	(*(PContext *)myObject).Catalog(drl.GetArray(), drl.GetSize());
	return drl.GetSize();
}

// Recatalog what's in our database
void
PContext::Recatalog()
{
	ClearSets();
	if (!Ready())
		return;
	pdb.ForEachBlock(cat_call, this);
}

// Get/release write lock on a catalog file
bool
PContext::WriteLock(char *lpath, const char *dbpath)
{
	static const char	msgfmt[] = "Catalog write lock held by process %d";
	static char		ourhost[128] = "@";
	
	if (lpath[0]) {				// release previous
		remove(lpath);
		lpath[0] = '\0';
	}
	if (dbpath == NULL || !*dbpath)
		return false;			// called for release
	char		message[256];
	int		msglen;
	const char *	ext = PgetSuffix(dbpath);
	if (ext == NULL)
		ext = dbpath + strlen(dbpath);
	else
		--ext;
	strcpy(lpath, dbpath);			// name lock file
	strcpy(lpath + (ext-dbpath), ".lock");
	if (ourhost[0]=='@' &&			// get host name
			gethostname(ourhost, sizeof(ourhost)) < 0)
		ourhost[0] = '\0';
						// check/set lock
	int	fd = open(lpath, O_WRONLY|O_CREAT|O_EXCL, 0644);
	if (fd < 0) {				// pre-existing lock?
		int	pid = -1;
		fd = open(lpath, O_RDONLY);	// read old lock file
		if (fd >= 0) {
			msglen = read(fd, message, sizeof(message));
			close(fd);
			msglen -= (msglen > 0);	// ignore eol
			if (msglen >= 0) {	// get process ID
				message[msglen] = '\0';
				if (sscanf(message, msgfmt, &pid) != 1)
					pid = 0;	// bogus lock file
			}
			fd = -1;
		}
		if (pid >= 0) {
			if (pid && (ext = strstr(message, " on ")) != NULL)
				msglen = ext - message + 4;
						// replace dead lock
			if ((!pid || (!strcmp(message+msglen, ourhost) &&
					kill(pid, 0) < 0)) &&
					remove(lpath) == 0)
				fd = open(lpath, O_WRONLY|O_CREAT|O_EXCL, 0644);
		}
		if (fd < 0) {			// lock active or error
			DMESGF(DMCwarning, "Cannot capture lock file '%s'", lpath);
			lpath[0] = '\0';
			return false;		// return failure
		}
		DMESG(DMCinfo, "Ignoring previous DB lock (dead process)");
	}
	sprintf(message, msgfmt, getpid());	// write lock message
	msglen = strlen(message);
	if (ourhost[0]) {
		sprintf(message+msglen, " on %s", ourhost);
		msglen += 4 + strlen(ourhost);
	}
	message[msglen++] = '\n';		// add eol & close
	write(fd, message, msglen);
	close(fd);
	return true;				// return success
}

// Set shadow database
bool
PContext::ShadowDB(const char *dbpath)
{
	if (sdb.Ready()) {			// close previous shadow
		sdb.Close();
		shd.Attach(NULL);
		sdbfs.close();
		WriteLock(slock_path, NULL);
	}
	if (dbpath == NULL || !*dbpath)
		return false;			// no shadow
	if (!WriteLock(slock_path, dbpath))
		return false;			// cannot get lock
	sdbfs.clear();
	sdbfs.open(dbpath, ios::in|ios::out|ios::binary);
	if (!sdbfs.is_open()) {
		sdbfs.clear();			// try truncating it
		sdbfs.open(dbpath, ios::in|ios::out|ios::trunc|ios::binary);
		if (!sdbfs.is_open()) {
			DMESGF(DMCresource, "Cannot open shadow database '%s'",
					dbpath);
			WriteLock(slock_path, NULL);
			return false;
		}
	}
	if (!shd.Attach(&sdbfs)) {
		DMESGF(DMCdata, "Header error in shadow database '%s'", dbpath);
		shd.Attach(NULL);
		sdbfs.close();
		WriteLock(slock_path, NULL);
		return false;
	}
	if (!sdb.Init(&shd)) {
		DMESGF(DMCdata, "Cannot initialize shadow database '%s'", dbpath);
		sdb.Init(NULL);
		shd.Attach(NULL);
		sdbfs.close();
		WriteLock(slock_path, NULL);
		return false;
	}
	return true;				// shadow DB ready
}

// Open new database and initialize
bool
PContext::OpenDB(const char *dbpath, bool ro)
{
	const char *	ext;
						// reset operations count
	sync_opcnt = sync_nops;
						// check for DB already open
	if ((dbpath != NULL) & (lock_path[0] != '\0') & !ro &&
			(ext = PgetSuffix(dbpath)) != NULL &&
			!strncmp(dbpath, lock_path, ext-dbpath)) {
		ClearUndo();			// syncs catalog to disk
		Recatalog();			// same effect as reopen
		return true;
	}
	CloseDB();				// close previous DB
	if (dbpath == NULL || !*dbpath)
		return false;
						// obtain write lock
	if (!ro && WriteLock(lock_path, dbpath)) {
		dbfs.clear();			// open for appending
		dbfs.open(dbpath, ios::in|ios::out|ios::binary);
		if (!dbfs.is_open()) {
			dbfs.clear();		// new file needs truncate?
			dbfs.open(dbpath, ios::in|ios::out|ios::trunc|ios::binary);
			if (!dbfs.is_open())	// permission problem?
				WriteLock(lock_path, NULL);
		}
	}
	if (!dbfs.is_open()) {			// otherwise open read-only
		dbfs.clear();
		dbfs.open(dbpath, ios::in|ios::binary);
		if (!dbfs.is_open()) {
			DMESGF(DMCresource, "Cannot open database '%s'", dbpath);
			return false;
		}
		if (!ro) {			// was hoping to write it
			DMESGF(DMCwarning,
				"Database '%s' opened read-only", dbpath);
			ro = true;
		}
	}
	if (!phd.Attach(&dbfs, ro)) {
		DMESGF(DMCdata, "Bad header in database '%s'", dbpath);
		phd.Attach(NULL);
		dbfs.close();
		WriteLock(lock_path, NULL);
		return false;
	}
	for (int fid = 0; fid < PDBFend; fid++)	// update DB fields
		if (phd.XlateID(&PDBFInfo,fid) < 0)
			phd.Field(PDBFInfo.GetName(fid), PDBFInfo.GetDetails(fid));
	bool	needsync = !ro && phd.HasChanged();
	phd.CopyFormat(&PDBFInfo);		// custom formatting functions
	phd.SetChangedFlag(needsync);		// avoids spurious header sync
	uniq = phd.Unique();			// copy unique boolean
	if (!pdb.Init(&phd)) {
		DMESGF(DMCdata, "Cannot initialize database '%s'", dbpath);
		pdb.Init(NULL);
		phd.Attach(NULL);
		dbfs.close();
		WriteLock(lock_path, NULL);
		return false;
	}
	dbs.SetDB(&pdb);			// set search to us
	strcpy(dbname, PgetFilename(dbpath));	// record DB name
	ext = PgetSuffix(dbname);
	if (ext != NULL)
		dbname[ext-dbname-1] = '\0';
	Recatalog();				// record catalog information
	if (!pdb.Ready()) {			// catalog scan failed?
		ClearSets();
		dbname[0] = '\0';
		dbs.SetDB(NULL);
		pdb.Init(NULL);
		phd.Attach(NULL);
		dbfs.close();
		WriteLock(lock_path, NULL);
		return false;
	}
	Reconcile();				// recompute rcommon
	return true;
}

// Close out database and free resources
bool
PContext::CloseDB()
{
	if (!dbname[0])				// database even open?
		return true;
	if (rdisp.GetSize() > 0) {		// preserve display records
		if (rdisp.GetFieldInfo() == &phd) {
			DBRecordList	savlst;	// header is going away...
			savlst.LinkAll(&rdisp);
			Clear(&phd);		// copies header info
			rdisp = savlst;		// copy back records
			if (rsrc[0] == PC_SOURCECHAR)
				strcpy(rsrc, PC_EXTSOURCE);
		} else {
			rcommon.Init();		// else clear common records
			rdisp_link.Init();
		}
	}
	ClearSets();				// clear sets
	dbs.SetDB(NULL);			// close database
	bool	ok = pdb.Close();
	phd.Attach(NULL);
	dbfs.close();
	WriteLock(lock_path, NULL);		// release lock
	if (!ok)
		DMESGF(DMCdata, "Error closing database '%s'", dbname);
	dbname[0] = '\0';
	return ok;
}
