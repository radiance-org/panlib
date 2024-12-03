/*
 *  dbaccess.cpp
 *  panlib
 *
 *  Database access class implementations.
 *
 *  Created by gward on Mon May 07 2001.
 *  Modified by plonghurst on Thursday Feb 15 2007.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#include "system.h"
#include "dbase.h"
#include "dmessage.h"

// pl: Necesarry under Windows
using std::ios;

#ifndef DB_MAXARR
#define DB_MAXARR	512			// biggest array we handle
#endif
						// field memory block size
#define DB_FBSIZ	(DB_FBLEN*sizeof(DBFieldS))
						// query hash set length
const uint32		DB_HASHSETLEN = primeAbove(32000);

#define DB_RBMAGIC	0x1a72			// database block magic #
#define DB_RBMAGIC_SWAP	0x721a			// swapped block magic #

// Structure to store fixed-size block of data in database
// This struct gets written out directly to disk with possible byte-swap
struct DBRecordBlock {
	int16		magic;			// magic number
	int16		nfields;		// record width
	int32		nrecords;		// records used
	int32		nrawbytes;		// raw bytes used
	DBFieldS	fieldarr[DB_FBLEN];	// field array
	const char *	GetEndData() const {
				return (const char *)(fieldarr+DB_FBLEN);
			}
	char *		EndData() {
				return (char *)(fieldarr+DB_FBLEN);
			}
	void		Init() {
				magic = DB_RBMAGIC;
				nfields = 0;
				nrecords = 0;
				nrawbytes = 0;
				memset(fieldarr, 0, DB_FBSIZ);
			}
	bool		BlockOK() const;
	bool		IsSwapped() const {
				return (magic == DB_RBMAGIC_SWAP);
			}
	int		NBytesFree() const {
				return DB_FBSIZ - nrawbytes -
					nrecords*nfields*sizeof(DBFieldS);
			}
	int		NRecordsFree() const {
				if (nfields <= 0) return DB_FBLEN/DB_MAXFIELD;
				if (nrecords <= 0) return DB_FBLEN/nfields;
				return NBytesFree()/(nfields*sizeof(DBFieldS) +
						nrawbytes/nrecords);
			}
	void *		Balloc(int nbytes) {
				if (nbytes & 03) nbytes += 4 - (nbytes&03);
				if (nbytes > NBytesFree())
					return NULL;
				nrawbytes += nbytes;
				return (void *)(EndData() - nrawbytes);
			}
	const DBFieldS *
			GetRecordB(int rn) const {
				if ((rn < 0) | (rn >= nrecords))
					return NULL;
				return fieldarr + nfields*rn;
			}
	DBFieldS *	NewRecordB() {
				if (NBytesFree() < nfields*(int)sizeof(DBFieldS))
					return NULL;
				DBFieldS * fp = fieldarr + nfields*nrecords++;
				memset(fp, 0, nfields*sizeof(DBFieldS));
				return fp;
			}
	const char *	StrAlloc(const char *sa[], int ns);
	void		DeleteRecords(const ABitMap &toDel);
	void		Compact();
	bool		Consistent(const DBFieldInfo *fi = NULL) const;
	bool		SwapBlock();
};

// Basic consistency check for record block
inline bool
DBRecordBlock::BlockOK() const
{
	if (magic != DB_RBMAGIC)
		return false;
	if ((nfields < 0) | (nfields > DB_MAXFIELD)
			| (nrecords < 0) | (nrawbytes < 0))
		return false;
	if ((nfields == 0) & (nrecords > 0))
		return false;
	if (nrecords*nfields*sizeof(DBFieldS) + nrawbytes > DB_FBSIZ)
		return false;
	return true;
}

// Swap bytes in a short, int, or long word
static inline void
swapBytes(char *bytes, int n)
{
	char	t;
	switch (n) {
	case 2:
		t = bytes[0]; bytes[0] = bytes[1]; bytes[1] = t;
		break;
	case 4:
		t = bytes[0]; bytes[0] = bytes[3]; bytes[3] = t;
		t = bytes[1]; bytes[1] = bytes[2]; bytes[2] = t;
		break;
	case 8:
		t = bytes[0]; bytes[0] = bytes[7]; bytes[7] = t;
		t = bytes[1]; bytes[1] = bytes[6]; bytes[6] = t;
		t = bytes[2]; bytes[2] = bytes[5]; bytes[5] = t;
		t = bytes[3]; bytes[3] = bytes[4]; bytes[4] = t;
		break;
	default:
		DMESG(DMCparameter, "Illegal data size in swapBytes");
	}
}

#define SwapBytes(w)	swapBytes((char *)&(w), sizeof(w))

// Attempt to merge strings with block memory or allocate space and copy
const char *
DBRecordBlock::StrAlloc(const char *sa[], int ns)
{
	if ((ns <= 0) | (sa == NULL))
		return NULL;
	int		i;
	int		len = 0;
	for (i = ns; i--; )			// compute total length
		if (sa[i] == NULL) {
			DMESG(DMCwarning, "NULL string in StrAlloc");
			sa[i] = ""; len++;
		} else
			len += strlen(sa[i]) + 1;
						// search for match
	const int	firstc = sa[0][0];
	char *		cp;
	char *		beg = EndData() - len;
	if (len & 0x3)
		beg -= 4 - (len&0x3);
	for ( ; beg >= EndData() - nrawbytes; beg -= 4) {
		if (*beg != firstc)
			continue;
		cp = beg;			// possible beginning
		for (i = 0; i < ns; i++) {
			if (istrcmp(cp, sa[i]))
				break;
			while (*cp++)
				;
		}
		if (i == ns)			// matched them all?
			return beg;		// then return pointer
	}
	beg = (char *)Balloc(len);		// else allocate
	if (beg == NULL)
		return NULL;			// no room at the inn
	cp = beg;
	for (i = 0; i < ns; i++) {		// OK, copy sequence
		strcpy(cp, sa[i]);
		while (*cp++)
			;
	}
	return beg;				// return location
}

// Delete indicated records from this block -- does not compact array data
void
DBRecordBlock::DeleteRecords(const ABitMap &toDel)
{
	if (toDel.Length() != nrecords) {	// check if valid
		DMESG(DMCparameter, "Bad call to DBRecordBlock::DeleteRecords");
		return;
	}
						// close gaps
	const void *		mp;
	DBFieldS *		fd = fieldarr;
	const DBFieldS *	fs = fieldarr;
	int			i, n;
	for (i = 0; i < toDel.Length(); i++) {
		if (toDel.Check(i)) {		// this one is going
			fs += nfields;
			--nrecords;
			continue;
		}
		if (fd == fs) {			// no gap yet
			fs = fd += nfields;
			continue;
		}
		for (n = nfields; n--; )	// else move record
			if ((mp = fs->GetPointer()) != NULL) {
				fd->SetOffsetArray((DBDataType)fs->dt, fs->nv, mp);
				fd++; fs++;
			} else
				*fd++ = *fs++;	// single value
	}
	if (nrecords <= 0)
		Init();				// emptied block
}

// Compact array data within this record block
void
DBRecordBlock::Compact()
{
	if (!BlockOK())				// basic check first
		return;
	if (!nrecords) {			// reset empty block
		if (nfields)
			Init();
		return;
	}
	if (nrawbytes <= 4)
		return;
	/*
	 *  Compute a bitmap of 32-bit data words currently in use.
	 *  Word 0 is at top of block, and increases downward.
	 */
	const int	nwords = (nrawbytes + 3) >> 2;
	ABitMap		usemap(nwords);
	const char *	ptr;
	int		i;
	for (i = nfields*nrecords; i--; ) {	// look at each field
		int	siz;
		ptr = (const char *)fieldarr[i].GetPointer();
		if (ptr == NULL)
			continue;
		siz = GetEndData() - ptr;
		if ((siz <= 0) | (siz > nrawbytes)) {
			DMESG(DMCwarning, "Bad offset pointer in Compact");
			continue;
		}
		siz = fieldarr[i].nv;
		if (siz <= 0) {
			DMESG(DMCwarning, "Illegal size in Compact");
			continue;
		}
		switch (fieldarr[i].dt) {
		case DBDTint:			// integer array size
			siz *= sizeof(int32);
			break;
		case DBDTfloat:			// float array size
			siz *= sizeof(float);
			break;
		case DBDTstring: {		// string array size
			const char *	cp = ptr;
			while (siz--)		// find sequence end
				while (*cp++)
					;
			siz = cp - ptr;
			} break;
		default:			// trouble
			DMESG(DMCparameter, "Missing type in Compact");
			return;
		}
		siz = (siz + 3) >> 2;		// round up to #words
		int	wrd = (GetEndData() - ptr + 3) >> 2;
		while (siz--)			// set "in use" bits
			usemap.Set(--wrd);
	}
	/*
	 *  Make a "save sequence" from the bitmap of used words:
	 *  a sequence of (#words_to_crush, #words_to_save),
	 *  starting from the top of our raw data and working down.
	 *  A crush pair of (0,0) terminates the list.
	 */
	unsigned short	saveseq[DB_FBSIZ/4/2 + 1][2];
	bool		saving = false;
	int		j = 0;
	saveseq[0][false] = saveseq[0][true] = 0;
	for (i = 0; i < nwords; i++) {
		if (saving) {
			if (!usemap.Check(i)) {
				++j;
				saveseq[j][false] = saveseq[j][true] = 0;
				saving = false;
			}
		} else
			saving = usemap.Check(i);
		saveseq[j][saving]++;
	}
	++j;					// terminate list
	saveseq[j][false] = saveseq[j][true] = 0;
	/*
	 * Perform actual compaction on data.
	 */
	char *		dstp = EndData();
	const char *	srcp = dstp;
	if (!saveseq[0][false]) {
		if (!saveseq[1][false])
			return;			// nothing to compact
		srcp -= saveseq[0][true] << 2;	// avoids no-op copy
		dstp -= saveseq[0][true] << 2;
		j = 1;				// skip to first crush
	} else
		j = 0;
	do {					// crush unused words
		srcp -= saveseq[j][false] << 2;
		for (i = saveseq[j][true] << 2; i--; )
			*--dstp = *--srcp;
	} while (saveseq[++j][false]);
	nrawbytes = EndData() - dstp;		// correct space use
	/*
	 * Adjust data offset pointers.
	 */
	for (i = nfields*nrecords; i--; ) {	// correct each array field
		ptr = (const char *)fieldarr[i].GetPointer();
		if (ptr == NULL)
			continue;		// must not be array
		srcp = dstp = EndData();
		j = 0;
		do {				// else find position
			srcp -= (saveseq[j][false] + saveseq[j][true]) << 2;
			dstp -= saveseq[j][true] << 2;
			if (ptr >= srcp) {	// found it, so correct it
				dstp += ptr - srcp;
				if (dstp != ptr &&
						!fieldarr[i].SetOffsetArray(
							(DBDataType)fieldarr[i].dt,
							fieldarr[i].nv, dstp))
					DMESG(DMCassert,
					"SetOffsetArray() failed in Compact()");
				break;
			}
		} while (saveseq[++j][false]);
	}
}

// Check disk image record block for consistency
bool
DBRecordBlock::Consistent(const DBFieldInfo *fi) const
{
	int			n;
	int			i;
	const DBFieldS *	fp;
						// basic checks
	if (!BlockOK())
		return false;
	if (nrecords == 0)			// OK if empty
		return true;
						// type checks
	if (fi != NULL) {
		if (nfields > fi->GetNFields())
			return false;
		for (fp = fieldarr, n = nrecords; n--; )
			for (i = 0; i < nfields; i++, fp++) {
				if (!fp->nv) continue;
				const DBFDetails *	dp = fi->GetDetails(i);
				if (dp == NULL)
					return false;
				if (dp->dtype == DBDTunknown) {
					switch (fp->dt) {
					case DBDTint:
					case DBDTfloat:
					case DBDTstring:
						break;
					default:
						return false;
					}
					if (fp->nv > 1 && !(fp->fl & DBFFarray))
						return false;
					continue;
				}
				if (fp->dt != (int)dp->dtype)
					return false;
				if (fp->nv > 1 && !(dp->flags & DBFFarray))
					return false;
			}
	} else {
		for (fp = fieldarr, i = nrecords*nfields; i--; fp++) {
			if (!fp->nv) continue;
			switch (fp->dt) {
			case DBDTint:
			case DBDTfloat:
			case DBDTstring:
				break;
			default:
				return false;
			}
			if (fp->nv > 1 && !(fp->fl & DBFFarray))
				return false;
		}
	}
						// array bounds checks
	const char * const	maxptr = GetEndData();
	const char * const	minptr = maxptr - nrawbytes;
	for (fp = fieldarr, i = nrecords*nfields; i--; fp++) {
		if (fp->nv > 1 || (fp->nv == 1) & (fp->dt == (int)DBDTstring)) {
			if ((fp->fl & (DBFFoffset|DBFFownmem)) != DBFFoffset)
				return false;
			if ((fp->v.offset & 03) != 0)
				return false;	// 32-bit alignment required
			const char *	ptr = (const char *)fp + fp->v.offset;
			int		siz = 0;
			if (ptr < minptr)
				return false;
			switch (fp->dt) {
			case DBDTint:
				siz = fp->nv*sizeof(int32);
				break;
			case DBDTfloat:
				siz = fp->nv*sizeof(float);
				break;
			case DBDTstring:
				if (fp->nv > 1 && !(fp->fl & DBFFstrseq))
					return false;
				siz = fp->nv;	// minimum
				break;
			}
			if (ptr + siz > maxptr)
				return false;
		}
	}
						// string terminator checks
	for (fp = fieldarr, i = nrecords*nfields; i--; fp++)
		if (fp->nv > 0 && fp->dt == (int)DBDTstring) {
			const char *	cp = (const char *)fp +
							fp->v.offset;
			for (n = fp->nv; n > 0; n -= !*cp++)
				if (cp >= maxptr)
					return false;
		}
	// We're not checking for array and string overlaps at this point
	return true;				// seems OK
}

// Swap bytes in a database record block, returning true on success
bool
DBRecordBlock::SwapBlock()
{
	const bool		fromSwapped = IsSwapped();
	int			i;
	const DBFieldS *	fp;
	
	if (fromSwapped) {			// pre-swap header & counts
		SwapBytes(magic);
		SwapBytes(nfields);
		SwapBytes(nrecords);
		SwapBytes(nrawbytes);
		for (fp = fieldarr, i = nrecords*nfields; i-- > 0; fp++)
			SwapBytes(fp->nv);
	}
	if (!BlockOK()) {			// perform basic checks
		DMESG(DMCparameter, "Basic block check failed in SwapBlock");
		return false;
	}
						// valid pointer range
	const char *	maxptr = GetEndData();
	const char *	minptr = maxptr - nrawbytes;
						// fix field values
	for (fp = fieldarr, i = nrecords*nfields; i-- > 0; fp++) {
		if (!fp->nv)
			continue;
		if (fp->dt == DBDTstring) {	// string (sequence)
			SwapBytes(fp->v.offset);
			continue;
		}
		if (fp->nv == 1)		// single value
			switch (fp->dt) {
			case DBDTint:
				SwapBytes(fp->v.mint);
				continue;
			case DBDTfloat:
				SwapBytes(fp->v.mfloat);
				continue;
			default:
				DMESG(DMCparameter,
					"Unknown data type in SwapBlock");
				return false;
			}
						// array data
		if (fromSwapped)
			SwapBytes(fp->v.offset);
		if (fp->dt == DBDTint) {
			int32 *	ip = (int32 *)((char *)fp + fp->v.offset);
			if (((char *)ip < minptr) | ((char *)(ip+fp->nv) > maxptr)) {
				DMESG(DMCparameter, "Bad offset in SwapBlock");
				return false;
			}
			for (int n = fp->nv; n--; ip++)
				SwapBytes(*ip);
		} else if (fp->dt == DBDTfloat) {
			float *	rp = (float *)((char *)fp + fp->v.offset);
			if (((char *)rp < minptr) | ((char *)(rp+fp->nv) > maxptr)) {
				DMESG(DMCparameter, "Bad offset in SwapBlock");
				return false;
			}
			for (int n = fp->nv; n--; rp++)
				SwapBytes(*rp);
		} else {
			DMESG(DMCparameter, "Unknown array type in SwapBlock");
			return false;
		}
		if (!fromSwapped)
			SwapBytes(fp->v.offset);
	}
	if (!fromSwapped) {			// post-swap counts & header
		for (fp = fieldarr, i = nrecords*nfields; i-- > 0; fp++)
			SwapBytes(fp->nv);
		SwapBytes(magic);
		SwapBytes(nfields);
		SwapBytes(nrecords);
		SwapBytes(nrawbytes);
	}
	return true;				// OK as far as we know
}

// Assign an offset array to a field
bool
DBFieldS::SetOffsetArray(DBDataType typ, int nva, const void *mp)
{
	if (nva < 1 || (nva == 1 && typ != DBDTstring) || mp == NULL)
		return false;
	Clear();
	fl = DBFFoffset|DBFFarray;
	nv = nva;
	if ((dt = typ) == DBDTstring)
		if (nva > 1)
			fl |= DBFFstrseq;
		else
			fl &= ~DBFFarray;
	v.offset = (const char *)mp - (char *)this;
	return true;
}

// Get generic pointer to ith data item
const void *
DBFieldS::GetPointer(int i) const
{
	if ((i < 0) | (i >= nv))
		return NULL;
	const char *	cp = (const char *)GetPointer();
	if (cp == NULL) {
		if (nv == 1)
			return (const void *)&v;
		return NULL;
	}
	if (fl & DBFFstrseq) {
		while (i--)
			while (*cp++)
				;
		return (const void *)cp;
	}
	switch (dt) {
	case DBDTint:
		cp += i*sizeof(int32);
		return (const void *)cp;
	case DBDTfloat:
		cp += i*sizeof(float);
		return (const void *)cp;
	case DBDTstring:
		if (nv == 1)
			return (const void *)cp;
		// FALL THROUGH
	default:
		DMESG(DMCparameter, "Botched type in DBFieldS::Get");
		break;
	}
	return NULL;
}

// Field reference operator (does not allocate arrays)
DBField &
DBField::operator=(const DBFieldS *ref)
{
	Clear();				// clear previous value
	if (!ref->nv)
		return *this;
	fl = ref->fl & ~(DBFFownmem|DBFFoffset);
	dt = ref->dt;
	if ((nv = ref->nv) == 1)		// check for inline value
		switch (dt) {
		case DBDTint:
			v.mint = ref->v.mint;
			return *this;
		case DBDTfloat:
			v.mfloat = ref->v.mfloat;
			return *this;
		}
	v.ptr = ref->GetPointer();
	DASSERT(v.ptr != NULL);
	return *this;				// all done
}

// Field copy operator (duplicates memory for arrays)
DBField &
DBField::operator=(const DBFieldS &src)
{
	Clear();				// clear previous value
	if (!src.nv)
		return *this;
	fl = src.fl & ~DBFFoffset;
	dt = src.dt;
	if ((nv = src.nv) == 1)			// check for inline value
		switch (dt) {
		case DBDTint:
			v.mint = src.v.mint;
			return *this;
		case DBDTfloat:
			v.mfloat = src.v.mfloat;
			return *this;
		}
	const void *	ptr = src.GetPointer();
	DASSERT(ptr != NULL);
						// allocate and copy memory
	switch (dt) {
	case DBDTint:
		v.aint = new int32 [nv];
		memcpy(v.aint, ptr, nv*sizeof(int32));
		break;
	case DBDTfloat:
		v.afloat = new float [nv];
		memcpy(v.afloat, ptr, nv*sizeof(float));
		break;
	case DBDTstring:
		if (nv == 1) {
			v.str = strlink((const char *)ptr);
		} else if (fl & DBFFstrseq) {
			const char *	cp = (const char *)ptr;
			v.astr = new const char * [nv];
			for (int i = 0; i < nv; i++) {
				v.astr[i] = strlink(cp);
				while (*cp++)
					;
			}
			fl &= ~DBFFstrseq;
		} else {
			const char * const *	sa =
						(const char * const *)ptr;
			v.astr = new const char * [nv];
			for (int i = nv; i--; )
				v.astr[i] = strlink(sa[i]);
		}
		break;
	default:
		DMESG(DMCparameter, "Botched type in DBFieldS assignment");
		v.ptr = NULL;
		return *this;
	}
	fl |= DBFFownmem;
	return *this;
}

// Close the database
bool
DBAccess::Close()
{
	if (header == NULL)
		return true;
	bool	aok = true;
	if (InUse()) {
		do
			ReleaseObject(objCache);
		while (InUse());
		DMESG(DMCparameter, "Database still in use when closed!");
		aok = false;
	}
	FreeMemory();
	aok &= (state == DBPSclean);
	if (next != NULL) {
		aok &= next->Close();
		delete next; next = NULL;
	} else
		aok &= header->Sync();
	header = NULL;
	state = DBPSclean;
	return aok;
}

// Initialize database access object
bool
DBAccess::Init(DBHeader *hdr, iostream *dbstrm, bool ro)
{
	if (!Close())				// close previous database
		DMESG(DMCdata, "Error closing previous database");
	if (hdr == NULL) {			// no header?
		strpos = 0; rwidth = 1;
		nrused = 0; nravail = 0;
		state = DBPSclean;
		return false;
	}
	if (dbstrm != NULL) {			// from new stream
		strpos = hdr->Attach(dbstrm, ro);
		if (strpos <= 0)
			hdr->Attach(dbstrm = NULL);
	} else {				// use header's stream
		dbstrm = hdr->GetStream();
		ro = hdr->ReadOnly();
		strpos = hdr->GetSize();
		if (strpos <= 0)
			dbstrm = NULL;
	}
	if (dbstrm == NULL)
		return false;
	bihset.NewBitMap(DB_HASHSETLEN);	// create hash index
	header = hdr;
	dbstrm->seekg(0, ios::end);		// get file size
	off_t		dblen = dbstrm->tellg();
	objSize = sizeof(DBRecordBlock);
	if (dblen < strpos + sizeof(DBRecordBlock)) {
		next = NULL;			// create empty block
		VirginBlock();
	} else if (dblen == strpos + sizeof(DBRecordBlock)) {
		next = NULL;			// on the last block
		UnknownBlock();
	} else {				// else link next block
		next = new DBAccess(this, dblen);
		UnknownBlock();
	}
	state = DBPSclean;			// ready to roll
	return true;
}

// Continue database initialization from previous block
DBAccess::DBAccess(const DBAccess *prev, off_t dblen)
{
	bihset.NewBitMap(DB_HASHSETLEN);	// create hash index
	header = prev->header;
	objSize = sizeof(DBRecordBlock);
	strpos = prev->strpos + sizeof(DBRecordBlock);
	if (dblen < strpos + sizeof(DBRecordBlock)) {
		next = NULL;			// create empty block
		VirginBlock();
	} else if (dblen == strpos + sizeof(DBRecordBlock)) {
		next = NULL;			// on the last block
		UnknownBlock();
	} else {				// else link next block
		next = new DBAccess(this, dblen);
		UnknownBlock();
	}
	state = DBPSclean;			// page not resident
}

// Print database statistics
void
DBAccess::PrintStats(ostream *ostr)
{
	if (header == NULL)
		return;
	int32		nfset[DB_MAXFIELD];
	int32		nfbytes[DB_MAXFIELD];
	int32		nrblocks = 0;
	int32		nru_total = 0;
	int32		nra_total = 0;
	int32		free_total = 0;
	int32		nrawb_total = 0;
	int		i;
	DBAccess *	ap;
						// initialize counts
	for (i = header->GetNFields(); i--; )
		nfset[i] = nfbytes[i] = 0;
						// gather our statistics
	for (ap = this; ap != NULL; ap = ap->next) {
		const DBRecordBlock *	rb;
		rb = (const DBRecordBlock *)ap->GetCacheObject();
		if (rb == NULL)
			continue;
		++nrblocks;
		nru_total += ap->NRecords();
		nra_total += ap->NRecordsAvail();
		free_total += rb->NBytesFree();
		nrawb_total += rb->nrawbytes;
		for (int j = ap->brl.GetSize(); j--; )
			for (i = ap->brl.Get(j).GetNAlloc(); i--; ) {
				const DBField &	f = ap->brl.Get(j)[i];
				if (f.GetNV() <= 0)
					continue;
				nfset[i] += f.GetNV();
				if ((f.GetNV() == 1) & (f.GetDT() != DBDTstring))
					continue;
				switch (f.GetDT()) {
				case DBDTint:
					nfbytes[i] += f.GetNV() * sizeof(int32);
					break;
				case DBDTfloat:
					nfbytes[i] += f.GetNV() * sizeof(float);
					break;
				case DBDTstring:
					{
					const char *	stra[DB_MAXARR];
					int		ns = f.Get(stra,DB_MAXARR);
					while (ns--)
						nfbytes[i] += strlen(stra[ns]) + 1;
					}
					break;
				}
			}
		ap->ReleaseObject(rb);
	}
						// print our statistics
	if (ostr == NULL)
		return;
	int32	array_total = 0;
	int32	waste_total = 0;
	for (i = header->GetNFields(); i--; ) {
		array_total += nfbytes[i];
		waste_total += (nru_total - nfset[i]) * 4;
	}
	*ostr << "==== Database Statistics ====\n";
	*ostr << "Header size:\t\t\t" << header->GetSize() << " bytes\n";
	*ostr << "Total record blocks:\t\t" << nrblocks
			<< " (" << nrblocks*sizeof(DBRecordBlock)/1024
			<< " Kbytes)\n";
	*ostr << "Total records used:\t\t" << nru_total << '\n';
	*ostr << "Total records free:\t\t" << nra_total << '\n';
	*ostr << "Unused space:\t\t\t" << free_total/1024 << " Kbytes\n";
	*ostr << "Unset field waste:\t\t" << waste_total/1024 << " Kbytes\n";
	*ostr << "Used array space:\t\t" << nrawb_total/1024 << " Kbytes\n";
	if (array_total > 0)
		*ostr << "Array data sharing:\t\t"
			<< 100.f - 100.f*nrawb_total/array_total << " %\n";
	*ostr << "Fields set:\n";
	for (i = 0; i < header->GetNFields(); i++) {
		if (!nfset[i])
			continue;
		*ostr << '\t' << nfset[i] << ' ' << header->GetName(i)
			<< " values";
		if (nfbytes[i] > 0)
			*ostr << " (" << nfbytes[i] << " array bytes)";
		*ostr << '\n';
	}
	*ostr << "Unused fields:";
	for (i = 0; i < header->GetNFields(); i++)
		if (!nfset[i])
			*ostr << ' ' << header->GetName(i);
	*ostr << '\n';
}

// Rebuild pseudo record list (block must be loaded)
void
DBAccess::ListRecords()
{
	const DBRecordBlock *	rb = (const DBRecordBlock *)objCache;
	DASSERT(rb != NULL);
	brl.Init(NULL, rb->nrecords);
	for (int i = rb->nrecords; i-- > 0; )
		if (!brl[i].Pseudo(header, rb->GetRecordB(i), rb->nfields))
			DMESG(DMCassert, "Pseudo() failed in ListRecords()");
}

// Recompute record index hash set (block must be loaded)
void
DBAccess::IndexRecords()
{
	bihset.ClearBitMap(false);		// clear hash set
	for (int i = brl.GetSize(); i-- > 0; ) {
		long	ih = brl.Get(i).HashIndex();
		if (ih < 0L) {
			bihset.ClearBitMap(true);
			return;			// this is useless
		}
		bihset.Set(ih % DB_HASHSETLEN);
	}
}

// Restore database page from disk
bool
DBAccess::RestoreMemory()
{
	if (IsResident())
		return (state != DBPSerror);
	objCache = Cmalloc(objSize = sizeof(DBRecordBlock));
	if (objCache == NULL)			// out of memory!
		return false;
	DBRecordBlock *	rb = (DBRecordBlock *)objCache;
	bool		need2index = false;
	iostream *	strm = header->GetStream();
	if (nrused == 0) {			// virgin block
		rb->Init();
		return true;
	}
	strm->clear();				// seek to block position
	if (strm->seekg(strpos).bad()) {
		DMESG(DMCsystem, "Cannot seek to block in database");
		goto errclean;
	}
						// load it from disk
	strm->read((char *)objCache, sizeof(DBRecordBlock));
	if ((size_t)strm->gcount() != sizeof(DBRecordBlock)) {
		DMESG(DMCsystem, "Database block read error");
		goto errclean;
	}
	if (rb->IsSwapped() != header->NeedByteSwap())
		DMESG(DMCwarning, "Unexpected byte swap in database block");
	if (rb->IsSwapped() && !rb->SwapBlock()) {
		DMESG(DMCdata, "Could not swap bytes in database block");
		goto errclean;
	}
	if (!rb->Consistent(header)) {
		DMESG(DMCdata, "Corrupted database disk block");
		goto errclean; 
	}
	if (rb->nfields > 0) {
		need2index = (nrused < 0);	// initial load?
		rwidth = rb->nfields;
		nrused = rb->nrecords;
		nravail = rb->NRecordsFree();
		ListRecords();			// rebuild pseudo record list
	} else
		VirginBlock();			// loaded virgin
	state = DBPSclean;
	if (need2index)
		IndexRecords();
	return true;
errclean:
	Cfree(objCache);
	objCache = NULL;
	UnknownBlock();
	state = DBPSerror;			// we're hosed
	return false;
}

// Sync this database block page to disk
bool
DBAccess::SyncPage(bool keep)
{
	if (state == DBPSclean)
		return true;			// already sync'ed
	if (state == DBPSerror)
		return false;			// already screwed
	if (ReadOnly())
		return false;			// no can do

	iostream *	strm = header->GetStream();
	DBRecordBlock *	rb = (DBRecordBlock *)objCache;
	DASSERT(objSize == sizeof(DBRecordBlock));
	DASSERT(rb != NULL);

	if (state >= DBPScompact) {		// compact if requested
		brl.Init();
		rb->Compact();
		if (keep) {			// good time to re-index
			ListRecords();
			IndexRecords();
		}
	}
	state = DBPSerror;			// prepare for catastrophe
	if (!rb->Consistent(header)) {
		DMESG(DMCdata, "Corrupted database memory block");
		goto cleanup;			// block is corrupt
	}
	if (header->NeedByteSwap()) {
		if (keep) {			// swap our own copy
			DBRecordBlock *	rnew =
				(DBRecordBlock *)malloc(sizeof(DBRecordBlock));
			if (rnew == NULL) {
				DMESG(DMCmemory, "malloc() failed in SyncPage");
				goto cleanup;
			}
			memcpy(rnew, rb, sizeof(DBRecordBlock));
			rb = rnew;
		}
		if (!rb->SwapBlock())
			goto cleanup;		// swap failed
	}
	strm->clear();
	if (strm->seekp(strpos).bad() ||
			strm->write((const char *)rb, sizeof(DBRecordBlock)).bad() ||
			strm->flush().bad()) {
		DMESG(DMCsystem, "Database block write error");
		goto cleanup;			// write error!
	}
	state = DBPSclean;			// else we're OK
cleanup:
	if (rb != (DBRecordBlock *)objCache)
		free(rb);			// free swapped copy
	return (state == DBPSclean);		// return success/failure
}

// Free database block (sync'ing to disk first of course)
void
DBAccess::FreeMemory()
{
	if (!IsPurgeable())			// free to free?
		return;
	if ((state != DBPSclean) & (state != DBPSerror))
		if (!SyncPage(false)) {		// sync to disk
			DMESG(DMCdata, "Could not sync database block");
			UnknownBlock();		// we're hosed
		}
	brl.Init();				// clear pseudo record list
	Cfree(objCache);			// free memory
	objCache = NULL;
}

// Get n records starting at rn, return number retrieved
int
DBAccess::GetRecords(DBRecord *rlist, int n, int rn)
{
	if ((n <= 0) | (rlist == NULL) | (rn < 0) | (header == NULL))
		return 0;
	if (rn + n > NRecords()) {		// break on border
		int	nfound = GetRecords(rlist, NRecords() - rn, rn);
		if (next == NULL)
			return nfound;
		return nfound + next->GetRecords(rlist + nfound, n - nfound,
						rn + nfound - NRecords());
	}
	if (GetCacheObject() == NULL)
		return 0;
						// copy each record
	for (int i = 0; i < n; i++)
		rlist[i] = brl.Get(rn+i);
						// release and return
	ReleaseObject(objCache);
	return n;
}

// Execute callback over all database blocks (in order if given)
int
DBAccess::ForEachBlock(int (*f)(const DBRecordList &, void *),
				void *client_data, const DBFieldSort *ord)
{
	if (f == NULL)
		return 0;
	int		rvsum;
	DBAccess *      dba;
	if (ord == NULL || ord->GetSort(0) < 0) {
		rvsum = 0;			// any order will do
		for (dba = this; dba != NULL; dba = dba->next) {
			if (dba->GetCacheObject() == NULL)
				continue;	// should be error?
			int	rv = (*f)(dba->brl, client_data);
			dba->ReleaseObject(dba->objCache);
			if (rv < 0)
				return rv;
			rvsum += rv;
		}
		return rvsum;
	}
	DBRecordList    allRec;			// load entire DB
	for (dba = this; dba != NULL; dba = dba->next) {
		const int       siz0 = allRec.GetSize();
		int		i;
		if (dba->GetCacheObject() == NULL)
			continue;		// should be error?
		allRec.Resize(siz0 + dba->NRecords());
		for (i = dba->NRecords(); i-- > 0; )
			allRec[siz0+i].Link(&dba->brl[i]);
	}
	allRec.sord = *ord;			// sort it
	allRec.Sort();
						// single block callback
	rvsum = (*f)(allRec, client_data);
	allRec.Init();				// release memory
	for (dba = this; dba != NULL; dba = dba->next)
		dba->ReleaseObject(dba->objCache);
	return rvsum;
}

// Local block -> record callback and data type
struct BlockData {
	int			(*rcb)(const DBRecord &, void *);
	void *			cdata;
};

// callback for ForEachRecord()
static int
block2record(const DBRecordList &rlist, void *client_data)
{
	BlockData *	bdp = (BlockData *)client_data;
	int		res = 0;
	int		i;

	for (i = 0; i < rlist.GetSize(); i++) {
		int	rv = (*bdp->rcb)(rlist.Get(i), bdp->cdata);
		if (rv < 0)
			return rv;
		res += rv;
	}
	return res;
}

// Execute callback over all database records (in order if given)
int
DBAccess::ForEachRecord(int (*f)(const DBRecord &, void *),
				void *client_data, const DBFieldSort *ord)
{
	if (f == NULL)
		return 0;
	BlockData	bd;
	bd.rcb = f;
	bd.cdata = client_data;

	return ForEachBlock(&block2record, &bd, ord);
}

// Write header and database records to a stream (outputs last record first)
void
DBAccess::Write(ostream *os, const char tabch)
{
	if ((os == NULL) | (header == NULL))
		return;
	if (next == NULL)			// write header if last block
		header->Write(os, tabch);
	else					// else write next block
		next->Write(os, tabch);
						// check stream status
	if (os->bad())
		return;
						// write this block
	if (GetCacheObject() == NULL)
		return;
						// write records from end
	for (int i = NRecords(); i-- > 0; )
		brl.Get(i).Write(os, tabch);
						// all done
	ReleaseObject(objCache);
}

// client data for write_blk() callback
struct write_blk_s {
	ostream *       os;
	char		tabch;
};

// Write out a database block (DBAccess::WriteSorted callback)
static int
write_blk(const DBRecordList &dl, void *dptr)
{
	write_blk_s *	wrp = (write_blk_s *)dptr;
	int		i;
						// check stream status
	if (wrp->os->bad())
		return 0;
						// write the block
	for (i = 0; i < dl.GetSize(); i++)
		dl.Get(i).Write(wrp->os, wrp->tabch);
	return i;
}

// Write out database records, sorted
void
DBAccess::WriteSorted(ostream *os, const char tabch)
{
						// write out header
	header->Write(os, tabch);
						// sort & write out records
	write_blk_s      wrs;
	wrs.os = os; wrs.tabch = tabch;
	ForEachBlock(&write_blk, &wrs, &header->sord);
}

// Read records from stream, adding to our DB
int
DBAccess::Read(istream *ins, char *tabcp, DBFieldInfo *fi)
{
	if (ins == NULL)
		return 0;
	if (ReadOnly()) {
		DMESG(DMCparameter, "Cannot add records to read-only DB");
		return 0;
	}
	char		tabch = '\0';		// default delimiter (any)
	DBFieldInfo	finf;
	int		nadd;
	if (tabcp == NULL)
		tabcp = &tabch;
	if (fi == NULL)
		fi = &finf;
						// read header if none given
	if (!fi->GetNFields() && (nadd = fi->Read(ins, tabcp)) <= 0)
		return nadd;
	DBAccess *	addto = this;
	DBRecord	rec(fi);		// init record holder
	rec.fieldlock = true;
	nadd = 0;
	while (ins->good()) {			// read & add each record
		int	nf = rec.Read(ins, *tabcp);
		if (nf < 0) {			// EOF or read error
			if (!ins->eof())
				DMESGF(DMCdata, "Read error at record %d", nadd+1);
			break;
		}
		if (!nf)
			continue;		// empty record
		while (addto->NRecordsAvail() <= 0 && addto->next != NULL)
			addto = addto->next;	// optimization
		if (!addto->AddRecord(rec))
			return -1;		// add error
		++nadd;
	}
	return nadd;
}

// Find n records matching query starting from *rnp
int
DBAccess::FindRecords(const DBQuery *dq, DBRecord *rlist, int n, int *rnp)
{
	if ((dq == NULL) | (n <= 0) | (rlist == NULL) | (header == NULL))
		return 0;			// nothing to do
	int	rcount = 0;
	if (rnp == NULL)			// make sure we have *rnp
		rnp = &rcount;
	else if (*rnp < 0)
		*rnp = 0;
						// check this block
	int	nfound = 0;
	bool	checkIt = (*rnp < NRecords());
						// filter on index hash
	if (checkIt & (dq->hashSet != NULL) && dq->hashSet->Length() == DB_HASHSETLEN) {
		ABitMap		matchSet = *dq->hashSet;
		matchSet &= bihset;
		checkIt = (matchSet.Find() != ABMend);
	} else if (checkIt & (dq->next == NULL))
		checkIt = !HashReject(dq->HashIndex());
						// advance through block
	if (checkIt && GetCacheObject() != NULL) {
		while (nfound < n && *rnp < NRecords()) {
						// copy if match
			if (brl.Get(*rnp).MatchQuery(dq))
				rlist[nfound++] = brl.Get(*rnp);
			++(*rnp);		// advance to next record
		}
		ReleaseObject(objCache);	// done with block
	} else if (*rnp < NRecords())
		*rnp = NRecords();
						// check subsequent blocks
	if ((nfound < n) & (next != NULL)) {
		*rnp -= NRecords();
		nfound += next->FindRecords(dq, rlist + nfound, n - nfound, rnp);
		*rnp += NRecords();
	}
	return nfound;				// return number found
}

// Add single record to this block if we can
bool
DBAccess::PutRecord(const DBRecord &rec)
{
	if (!rec.GetNAssigned()) {
		DMESG(DMCwarning, "Pretending to add empty record to database");
		return true;
	}
	if (rec.GetNAlloc() > RecordWidth())
		return false;			// too wide for this block
	if (NRecordsAvail() <= 0)
		return false;			// space is too tight
	if (InUse())
		DMESG(DMCwarning, "Adding record to block in use");
	GetCacheObject();			// else lock down block
	DBRecordBlock *	rb = (DBRecordBlock *)objCache;
	if (rb == NULL)
		return false;
	DASSERT(objSize == sizeof(DBRecordBlock));
	if (rb->nfields <= 0)
		rb->nfields = RecordWidth();
	int		orig_nrawbytes = rb->nrawbytes;
	DBFieldS *	df = rb->NewRecordB();	// allocate block record
	if (df == NULL) {
		ReleaseObject(objCache);
		nravail = 0;
		return false;			// out of space
	}
	DASSERT(rec.GetFieldInfo() != NULL);
	DASSERT(header != NULL);
	long		ih;
	int		idmap[DB_MAXFIELD];
	header->MapIDs(idmap, rec.GetFieldInfo());
	for (int i = 0; i < rec.GetNAlloc(); i++) {
		const DBField &	f = rec[i];
		if (f.GetNV() <= 0)		// check field
			continue;
		const int		j = idmap[i];
		const DBFDetails *	d = header->GetDetails(j);
		if (d == NULL ||
			(d->dtype != DBDTunknown && d->dtype != f.GetDT()) ||
				(f.GetNV() > 1 && !(d->flags & DBFFarray))) {
			DMESGF(DMCdata, "Cannot assign field %s in PutRecord",
					rec.GetFieldName(i));
			continue;
		}
		if (d->dtype == DBDTstring) {	// allocate string(s)
			const char *	str[DB_MAXARR];
			int		ns = f.Get(str, DB_MAXARR);
			if (ns < f.GetNV())
				DMESGF(DMCparameter,
				"Too many strings in field %s in PutRecord",
						rec.GetFieldName(i));
			const char *	sp = rb->StrAlloc(str, ns);
			if (sp == NULL)
				goto failure;
			if (!df[j].SetOffsetArray(DBDTstring, ns, sp))
				goto failure;
			continue;
		}
		if (f.GetNV() == 1) {		// single value
			df[j].fl = 0;
			df[j].nv = 1;
			switch (df[j].dt = f.GetDT()) {
			case DBDTint:
				f.Get(&df[j].v.mint);
				break;
			case DBDTfloat:
				f.Get(&df[j].v.mfloat);
				break;
			default:
				DMESG(DMCassert, "Type botch in PutRecord()");
			}
		} else if (d->dtype == DBDTint) {	// integer array
			int32 *	ip = (int32 *)
					rb->Balloc(f.GetNV()*sizeof(int32));
			if (ip == NULL)
				goto failure;
			if (f.Get(ip, f.GetNV()) != f.GetNV())
				goto failure;
			if (!df[j].SetOffsetArray(DBDTint, f.GetNV(), ip))
				goto failure;
		} else if (d->dtype == DBDTfloat) {	// float array
			float *	fp = (float *)
					rb->Balloc(f.GetNV()*sizeof(float));
			if (fp == NULL)
				goto failure;
			if (f.Get(fp, f.GetNV()) != f.GetNV())
				goto failure;
			if (!df[j].SetOffsetArray(DBDTfloat, f.GetNV(), fp))
				goto failure;
		} else
			DMESG(DMCparameter, "Missing type in PutRecord");
	}
						// end here on success
	if (!brl.NewRecord().Pseudo(header, df, rb->nfields))
		DMESG(DMCassert, "Pseudo() failed in PutRecord()");
	ih = brl.Get(nrused).HashIndex();	// add to index hash set
	if (ih < 0L)
		bihset.ClearBitMap(true);
	else
		bihset.Set(ih % DB_HASHSETLEN);
	++nrused;
	nravail = rb->NRecordsFree();
	if (state < DBPSsync)
		state = DBPSsync;		// mark block as dirty
	ReleaseObject(objCache);		// release memory
	return true;
failure:					// go here if not enough space
	rb->nrawbytes = orig_nrawbytes;
	if (--rb->nrecords <= 0)
		DMESG(DMCassert, "Record bigger than database block");
	nravail = 0;
	ReleaseObject(objCache);		// release memory
	return false;
}

// Add new records to database
int
DBAccess::AddRecords(const DBRecord *rlist, int n)
{
	if ((rlist == NULL) | (n <= 0))
		return 0;
	if (ReadOnly()) {
		DMESG(DMCdata, "AddRecords called on read-only DB");
		return 0;
	}
	int	nadd = 0;			// add what we can to this block
	while (PutRecord(rlist[nadd]))
		if (++nadd >= n)
			return nadd;
	if (next == NULL) {			// create new block if needed
		if (rlist[nadd].GetNAlloc() > header->GetNFields()) {
			DMESG(DMCparameter, "Too many fields in AddRecords");
			return nadd;
		}
		next = new DBAccess(this);
	}
						// add to next block
	return nadd + next->AddRecords(rlist + nadd, n - nadd);
}

// Delete all matching records from database
int
DBAccess::DeleteRecords(const DBRecord *rlist, int n)
{
	if ((n <= 0) | (rlist == NULL))
		return 0;
	if (ReadOnly()) {
		DMESG(DMCdata, "DeleteRecords called on read-only DB");
		return 0;
	}
	int	ndeleted = 0;
	if (next != NULL)			// delete from other blocks
		ndeleted = next->DeleteRecords(rlist, n);
	if (NRecords() <= 0)
		return ndeleted;
	ABitMap		toDo(n);		// check if any to do here
	uint32		i;
	for (i = n; i--; )
		if (rlist[i].GetNAssigned() && !HashReject(rlist[i].HashIndex()))
			toDo.Set(i);
	if (toDo.Find() == ABMend)		// nothing in this block
		return ndeleted;
	if (InUse())
		DMESG(DMCwarning, "Deleting records from block in use");
	GetCacheObject();			// delete from this block
	DBRecordBlock *	rb = (DBRecordBlock *)objCache;
	if (rb == NULL)
		return ndeleted;		// should be error?
	DASSERT(objSize == sizeof(DBRecordBlock));
	ABitMap		toDel(rb->nrecords);	// search for records to delete
	for (int j = rb->nrecords; j--; )
		for (i = 0; toDo.Find(&i); i++) {
			if (rlist[i] != brl[j])
				continue;	// target exact matches only
			toDel.Set(j);
			++ndeleted;
			break;			// this record is going
		}
	rb->DeleteRecords(toDel);		// delete records from block
	if (!rb->nrecords) {			// like a virgin?
		VirginBlock();
		state = DBPSsync;
		brl.Init();
	} else if (nrused != rb->nrecords) {	// else update our block
		nrused = rb->nrecords;
		nravail = rb->NRecordsFree();
		if (state < DBPScompact)
			state = DBPScompact;
		ListRecords();
	}
	ReleaseObject(objCache);		// done with block
	return ndeleted;
}
