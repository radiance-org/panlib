/*
 *  dbase.h
 *  panlib
 *
 *  Depends on <iostream.h> and <string.h>
 *
 *  Classes for flat database management
 *
 *  Created by gward on Sun May 06 2001.
 *  Modified by plonghurst on Thursday Feb 15 2007.
 *  Copyright (c) 2004 Anyhere Software. All rights reserved.
 *
 */

#ifndef _DBASE_H_
#define _DBASE_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;
#include <string.h>
#include "cache.h"
#include "pstrings.h"
#include "dbhlink.h"
#include "abitmap.h"

/*
 *  This database manager uses fixed block sizes (DB_FBLEN)
 *  to simplify cache management while still allowing string
 *  and array fields.  Additions to the field list for a
 *  given database are accommodated through the DBFieldID
 *  class and its derivatives.  The names and integer IDs
 *  corresponding to previously defined fields are never
 *  changed from the time the database is created, but new
 *  ones may be added as desired.  Typically, the program
 *  opens an existing database, getting field definitions
 *  from its header.  It then checks these with its own
 *  notion of these fields and gets the corresponding IDs.
 *
 *  Separate databases with different field definitions
 *  require separate files.  Other than the unique header
 *  at the beginning of the file, there is no notion of
 *  sections in a database file -- it's just a collection
 *  of fixed-size record blocks.  Fields may be any length
 *  and have any number of members up to what will fit in
 *  a single record block (64 Kbytes).  Supported field types
 *  are 32-bit integer, 32-bit IEEE float, and nul-terminated
 *  (8-bit) character strings.  This format is NOT
 *  well-suited to binary fields, unless they are small
 *  enough to be encoded as a short hex string.  Larger data
 *  fields should go in separate files referred to by name.
 */

// Change DB_FBLEN if you want to break everything...
#define DB_FBLEN	8192			// fields allocated per block

#define DB_MAXFIELD	128			// max. fields per record
#define DB_MNAMELEN	16			// mean field name length

// The following class maps database field names to IDs
class DBFieldID {
protected:
	int		nfields;		// number of fields
	const char *	mname[DB_MAXFIELD];	// field name pointers
	short		mlexo[DB_MAXFIELD];	// sorted name order
	int		NewField(const char *mnam);
public:
			DBFieldID() {
				nfields = 0;
			}
			DBFieldID(const DBFieldID &orig) {
				nfields = 0;
				*this = orig;
			}
			~DBFieldID() {
				Clear();
			}
			// Clear our field list
	void		Clear() {
				while (nfields > 0) strunlink(mname[--nfields]);
			}
			// Return the number of fields
	int		GetNFields() const {
				return nfields;
			}
			// Get field ID from name
	int		GetID(const char *mnam) const;
			// Get field name from ID
	const char *	GetName(int mid) const {
				if ((mid < 0) | (mid >= nfields)) return NULL;
				return mname[mid];
			}
			// Get or generate field ID for name
	int		Field(const char *mnam) {
				int	mid = GetID(mnam);
				if (mid >= 0) return mid;
				return NewField(mnam);
			}
			// Rename the given field
	bool		Rename(int mid, const char *mnam);
			// Map IDs from another field list
			// Return true if map is not 1-1
	bool		MapIDs(int idmap[DB_MAXFIELD], const DBFieldID *srcp) const;
			// Quick translation for a single field ID
	int		XlateID(const DBFieldID *srcp, int mid) const {
				if ((this == NULL) | (srcp == NULL)) return -1;
				if ((mid < 0) | (mid >= srcp->nfields)) return -1;
				if (srcp == this) return mid;
				if (mid < nfields && mname[mid] == srcp->mname[mid])
					return mid;
				return GetID(srcp->mname[mid]);
			}
			// Write tab-separated field names to stream
	void		Write(ostream *os, const char tabch = '\t') const {
				if (os == NULL) return;
				for (int i = 0; i < nfields; ) {
					*os << mname[i++];
					if (i < nfields) *os << tabch;
				} *os << '\n';
			}
			// Read delimited field list from stream
	int		Read(istream *ins, char *tabcp = NULL);
			// Memberwise copy operator
	DBFieldID &	operator=(const DBFieldID &src);
};

#define	DB_MAXSORTF	15			// maximum # sort fields

// Class to track sorting order
class DBFieldSort {
	friend class	DBFieldInfo;
private:
	int		order[DB_MAXSORTF+1];	// sort priority
public:
			DBFieldSort() {
				Clear();
			}
			// Clear sort
	void		Clear() {
				memset(order, 0, sizeof(order));
			}
			// Add field to sort priority list
	bool		AddSort(int mid, bool reverse = false);
			// Push field to top of sort priority list
	bool		PushSort(int mid, bool reverse = false);
			// Get ID of field in sort list or -1 at end
	int		GetSort(int i, bool *reverse = NULL) const;
};

// Get ID of field in sort list or -1 at end
inline int
DBFieldSort::GetSort(int i, bool *reverse) const
{
	if (i >= DB_MAXSORTF || !order[i])
		return -1;
	if (order[i] < 0) {
		if (reverse != NULL) *reverse = true;
		return -order[i]-1;
	}
	if (reverse != NULL) *reverse = false;
	return order[i]-1;
}

// Field flags and masks
// Add new flags to this list, but DO NOT alter existing bit values.
enum	DBFieldFlag	{
			DBFFarray=0x1,		// array allowed/set
			DBFFownmem=0x2,		// we own aux. memory
			DBFFstrseq=0x4,		// strings in sequence
			DBFFoffset=0x8,		// array uses offset
			DBFFshow=0x100,		// display this field
			DBFFhide=0x200,		// hide field normally
			DBFFuser=0x400,		// user-alterable field
			DBFFrequired=0x800,	// required field
			DBFFindex=0x1000,	// index field
			DBFFcustom=0x2000,	// user-defined field
			DBFFstoreM=0xf,		// storage flag mask
			DBFFinfoM=0x3f01	// info flag mask
};

#define DB_MAXDESCR	64			// max. description length

class DBField;					// defined later in header

// Standard field formatting functions
char *		DBformatArray(char *, int, const DBField &);
char *		DBformatUnknown(char *, int, const DBField &);
char *		DBformatInteger(char *, int, const DBField &);
char *		DBformatFloat(char *, int, const DBField &);
char *		DBformatString(char *, int, const DBField &);

// Struct to contain field details
struct DBFDetails {
	DBDataType	dtype;			// data type
	int		flags;			// info flags
	char		descrip[DB_MAXDESCR];	// description
						// display formatting function
	char *		(*DisplayFormat)(char *buf, int len, const DBField &df);
};

// Subclass to include associated field details
class DBFieldInfo : public DBFieldID {
protected:
	DBFDetails	details[DB_MAXFIELD];	// field details
public:
	DBFieldSort	sord;			// sort priority
			DBFieldInfo() {
				memset(details, 0, sizeof(details));
			}
			DBFieldInfo(const DBFieldInfo &orig) {
				*this = orig;
			}
	virtual		~DBFieldInfo() {}
			// Clear fields
	void		Clear() {
				sord.Clear();
				DBFieldID::Clear();
				memset(details, 0, sizeof(details));
			}
			// Get details for the given field ID
	const DBFDetails *
			GetDetails(int mid) const {
				if (this == NULL) return NULL;
				if ((mid < 0) | (mid >= nfields)) return NULL;
				return &details[mid];
			}
			// Get details for the given field name
	const DBFDetails *
			GetDetails(const char *mnam) const {
				if (this == NULL) return NULL;
				return GetDetails(GetID(mnam));
			}
			// Set flag(s) on or off for given field
	void		SetFlags(int mid, int fl, bool setting = true) {
				if ((mid < 0) | (mid >= nfields)) return;
				if (!(fl &= DBFFinfoM)) return;
				if (setting) details[mid].flags |= fl;
				else details[mid].flags &= ~fl;
			}
			// Clear flag(s)
	void		ClearFlags(int mid, int fl) {
				SetFlags(mid, fl, false);
			}
			// Define or redefine a field and its details
	int		Field(const char *mnam, const DBFDetails *fd) {
				if (fd == NULL) return -1;
				int	mid = DBFieldID::Field(mnam);
				if (mid < 0) return -1;
				details[mid] = *fd;
				return mid;
			}
	int		Field(const char *mnam, DBDataType dtype,
				const char *descrip = NULL, int flags = 0,
				char *(*df)(char *, int, const DBField &) = NULL);
			// Copy corresponding formatting functions
	void		CopyFormat(const DBFieldInfo *fi);
			// Read tab-separated field list from stream
	int		Read(istream *ins, char *tabcp = NULL) {
				Clear();
				return DBFieldID::Read(ins, tabcp);
			}
			// Write field information to a stream
	virtual int	WriteInfo(ostream *ostr, int tab = 0) const;
			// Get upper bound on field information size
	virtual int	InfoSize() const {
				return 100 +
					DB_MAXFIELD*(DB_MNAMELEN+DB_MAXDESCR+64);
			}
			// Read field information from a stream
	virtual int	ReadInfo(istream *istr);
			// Copy operator
	DBFieldInfo &	operator=(const DBFieldInfo &src);
};

// Subclass for managing database header information
class DBHeader : public DBFieldInfo {
protected:
	int		hdrsiz;			// header size (bytes)
	int		bigendfile;		// DBF is bigendian?
	int		fieldsize;		// DBF's sizeof(DBFieldS)
	int		fieldsperblock;		// DBF's DB_FBLEN
	char		soft[128];		// DBF's software
	mutable int32	lastsync;		// hash at last sync
	bool		readonly;		// file is read only?
	iostream *	strm;			// database file pointer
	DBHeaderLink *	hvars;			// variable handler
public:
	static const bool
			bigendian;		// native bigendian?
			DBHeader();
	virtual		~DBHeader() {
				Sync();
				delete hvars;
			}
			// Attach header with the specified i/o stream
	int		Attach(iostream *dbstream, bool ro=false);
			// Return attached stream pointer
	iostream *	GetStream() const {
				return strm;
			}
			// Is header stream read-only?
	bool		ReadOnly() const {
				return readonly;
			}
			// Get on-disk header size
	int		GetSize() const {
				return hdrsiz;
			}
			// Is disk file byte order bigendian (Intel)?
	bool		BigEndianFile() const {
				return (bool)bigendfile;
			}
			// Is disk byte order different from memory?
	bool		NeedByteSwap() const {
				return (bigendian != (bool)bigendfile);
			}
			// Update disk to match in-memory information
	bool		Sync() const;
			// Changed flag (SetChanged(false) called after sync)
	void		SetChangedFlag(bool changed = true) {
				lastsync = (changed ? int32(-1) : Hash());
			}
			// Has the header changed since the last sync?
	bool		HasChanged() const {
				return (lastsync < 0 || lastsync != Hash());
			}
			// Compute hash value for current settings
	virtual int32	Hash() const {
				return (int32)memhash(this,
					(char *)&lastsync - (char *)this, 31);
			}
			// Write out "magic string" at file beginning
	bool		PutMagic() const;
			// Write header information to stream
	virtual int	WriteInfo(ostream *ostr, int tab = 0) const {
				return hvars->WriteVars(ostr, tab);
			}
			// Return size of information header
	virtual int	InfoSize() const {
				return hvars->GetMaxLength();
			}
			// Check "magic string" at file beginning
	bool		CheckMagic() const;
			// Read header information from stream
	virtual int	ReadInfo(istream *istr) {
				hvars->Clear();
				int	c = hvars->ScanVars(istr);
				if (c < 0) return c;
				if (!CheckValues()) return -1;
				return c;
			}
			// Check that our DB header values are legal
	virtual bool	CheckValues() const;
			// Set all values -- be sure to clear first!
	virtual void	SetDefaults();
};

// Struct to contain on-disk database field values (8 bytes + colocated memory)
struct DBFieldS {
	uint8		dt;			// data type
	uint8		fl;			// storage flags
	uint16		nv;			// number of values
	union {
		int32		mint;			// integer value
		float		mfloat;			// real value
		int32		offset;			// offset to array data
	}		v;			// field value
			DBFieldS() {
				Clear();
			};
			// Clear any assignment
	void		Clear() {
				dt = DBDTunknown; fl = 0; nv = 0; v.offset = 0;
			}
			// Assign an offset array
	bool		SetOffsetArray(DBDataType typ, int nva, const void *mp);
			// Return true if field has array or string data
	bool		HasPointer() const {
				if (!nv || (nv == 1) & (dt != DBDTstring))
					return false;
				return (v.offset != 0);
			}
			// Get generic pointer to array data
	const void *	GetPointer() const {
				if (!HasPointer())
					return NULL;
				return (const void *)
						((const char *)this + v.offset);
			}
			// Get generic pointer to ith data item
	const void *	GetPointer(int i) const;
};

// Class to contain database field values (8/12 bytes + array memory)
class DBField {
protected:
	uint8		dt;			// data type
	uint8		fl;			// storage flags
	uint16		nv;			// number of values
	union {
		int32		mint;			// integer value
		float		mfloat;			// real value
		const char *	str;			// nul-terminated string
		int32 *		aint;			// array of ints
		float *		afloat;			// array of floats
		const char **	astr;			// array of strings
		int32		offset;			// offset to array data
		const void *	ptr;			// generic pointer
	}		v;			// field value
			// Reference constructor (does not allocate arrays)
			DBField(const DBField *origp) {
				nv = 0;
				*this = origp;	// array reference
			}
			// Temp reference operator (does not allocate arrays)
	DBField &	operator=(const DBField *ref);
			// Return true if field has array or string data
	bool		HasPointer() const {
				if (!nv || (nv == 1) & (dt != DBDTstring))
					return false;
				if (fl & DBFFoffset)
					return (v.offset != 0);
				return (v.ptr != NULL);
			}
			// Get generic pointer to array data
	const void *	GetPointer() const {
				if (!HasPointer())
					return NULL;
				if (fl & DBFFoffset)
					return (const void *)
						((const char *)this + v.offset);
				return v.ptr;
			}
public:
			DBField() {
				dt = DBDTunknown; fl = 0; nv = 0; v.ptr = NULL;
			};
			DBField(const DBField &orig) {
				nv = 0;
				*this = orig;	// copies memory
			}
			DBField(const DBFieldS &orig) {
				nv = 0;
				*this = orig;	// copies memory
			}
			DBField(int32 iv) {
				nv = 0;
				*this = iv;
			}
			DBField(float fv) {
				nv = 0;
				*this = fv;
			}
			DBField(double dv) {
				nv = 0;
				*this = dv;
			}
			DBField(const char *sv) {
				nv = 0;
				*this = sv;	// string reference
			}
			DBField(const DBFieldS *origp) {
				nv = 0;
				*this = origp;	// array reference
			}
			~DBField() {
				Clear();
			}
			// Get number of values
	int		GetNV() const {
				return nv;
			}
			// Get data type
	DBDataType	GetDT() const {
				return (DBDataType)dt;
			}
			// Clear any assignment (freeing memory if ours)
	void		Clear();
			// Get generic pointer to ith data item
	const void *	GetPointer(int i) const;
			// Compute hash on field value
	unsigned long	Hash(int nbits = 32) const;
			// Type-checking value calls, return # values
	int		Get(int32 ia[], int n = 1) const;
	int		Get(float fa[], int n = 1) const;
			// Beware:  Get(const char **) shares memory!
	int		Get(const char *sa[], int n = 1) const;
			// Copy ith string
	bool		Get(char str[], int len, int i = 0) const;
			// Copy ith value into a singular field
	bool		Get(DBField *fv, int i = 0) const;
			// Get ith value calls -- type cast w/ zero fill
	int32		GetInt(int i = 0) const;
	float		GetFloat(int i = 0) const;
	const char *	GetString(int i = 0) const;
			// Typed assignment calls
	void		Set(int32 iv) {
				Clear();
				fl = 0; dt = DBDTint; nv = 1;
				v.mint = iv;
			}
	DBField &	operator=(int32 iv) {
				Set(iv); return *this;
			}
	void		Set(float fv);
	DBField &	operator=(float fv) {
				Set(fv); return *this;
			}
	void		Set(double dv) {
				Set((float)dv);
			}
	DBField &	operator=(double dv) {
				Set((float)dv); return *this;
			}
	void		Set(const char *sv) {
				Clear();
				if (sv == NULL) return;
				fl = DBFFownmem; dt = DBDTstring; nv = 1;
				v.str = strlink(sv);
			}
			// Beware:  assigning const string uses reference!
	DBField &	operator=(const char *sv) {
				Clear();
				if (sv == NULL) return *this;
				fl = 0; dt = DBDTstring; nv = 1;
				v.str = sv;
				return *this;
			}
	void		Set(const int32 ia[], int n);
	void		Set(const float fa[], int n);
	void		Set(const char * const sa[], int n);
			// Put text into string array (nl separates elements)
	int		SetText(const char *lines, const char nl = '\n');
			// Sort array values
	bool		SortArray();
			// Copy operator (duplicates memory for arrays)
	DBField &	operator=(const DBField &src);
	DBField &	operator=(const DBFieldS &src);
			// Reference operator for "short" field data
	DBField &	operator=(const DBFieldS *ref);
			// Compare to another field, >0 if this greater
	int		Compare(const DBField &that) const;
			// Compare fields for full equality (arrays must match)
	bool		operator==(const DBField &that) const;
			// Typed (range) match calls, true for any array match
	bool		MatchRange(int32 ilo, int32 ihi) const;
	bool		Match(int32 i) const {
				return MatchRange(i, i);
			}
	bool		MatchRange(float flo, float fhi) const;
	bool		Match(float f) const {
				return MatchRange(f, f);
			}
	bool		MatchRange(const char *slo, const char *shi) const;
	bool		Match(const char *s) const {
				return MatchRange(s, s);
			}
	bool		MatchRange(const DBField &lo, const DBField &hi) const;
	bool		Match(const DBField &dbf) const {
				return MatchRange(dbf, dbf);
			}
};

// Get one or more integer values
inline int
DBField::Get(int32 ia[], int n) const
{
	if (!nv | (n <= 0) || dt != DBDTint)
		return 0;
	if (n > nv)
		n = nv;
	if (nv == 1)
		ia[0] = v.mint;
	else
		memcpy(ia, GetPointer(), n*sizeof(int32));
	return n;
}

// Get one or more float values
inline int
DBField::Get(float fa[], int n) const
{
	if (!nv | (n <= 0) || dt != DBDTfloat)
		return 0;
	if (n > nv)
		n = nv;
	if (nv == 1)
		fa[0] = v.mfloat;
	else
		memcpy(fa, GetPointer(), n*sizeof(float));
	return n;
}

// Convenient operators for field comparisons
inline bool
operator<(const DBField &lf, const DBField &rt)
{
	return (lf.Compare(rt) < 0);
}
inline bool
operator>(const DBField &lf, const DBField &rt)
{
	return (lf.Compare(rt) > 0);
}
inline bool
operator<=(const DBField &lf, const DBField &rt)
{
	return (lf.Compare(rt) <= 0);
}
inline bool
operator>=(const DBField &lf, const DBField &rt)
{
	return (lf.Compare(rt) >= 0);
}
inline bool
operator!=(const DBField &lf, const DBField &rt)
{
	return !(lf == rt);
}

// Append array values from another field
DBField &	operator+=(DBField &dst, const DBField &src);
// Append different array values from another field
DBField &	operator|=(DBField &dst, const DBField &src);
// Array field intersecton operator
DBField &	operator&=(DBField &dst, const DBField &src);
// Array field difference operator
DBField &	operator-=(DBField &dst, const DBField &src);

// I/O stream operators for fields
istream &	operator>>(istream &ci, DBField &f);
ostream &	operator<<(ostream &co, const DBField &f);

#define DB_NHASHBITS	31			// # bits in index hash value

enum {DBHVunknown=-2, DBHVnone=-1};		// hash states

// Structure to hold linked record information
struct DBRLink {
	DBField *	field;			// field array pointer
	short		nfields;		// number of fields in array
	short		readonly;		// read-only fields?
	int		nlinks;			// link count (>=1)
	long		ihash;			// index hash value
			DBRLink(int siz = 0) {
				nlinks = 1; readonly = -1;
				Init(siz);
			}
			~DBRLink() {
				Init();
			}
	void		Init(int n = 0) {
				if (readonly >= 0) delete [] field;
				if (n <= 0) { field = NULL; nfields = 0; }
				else field = new DBField [nfields=n];
				readonly = 0; ihash = DBHVunknown;
			}
	void		GrowTo(int n);
};

struct DBQuery;

// Database record holder, designed to work in array for qsort()
class DBRecord {
	friend class	DBRecordList;
protected:
	DBRLink *	l;			// linked information
	const DBFieldInfo *
			finfo;			// field info pointer
	void		Empty() {
				finfo = NULL; fieldlock = false;
				l = NULL; sortord = NULL;
			}
	void		Free() {
				if (l == NULL || --l->nlinks) return;
				delete l;
			}
	long		ComputeHashIndex() const;
	bool		PrepField(int fn, DBDataType dt, int nv = 1);
public:
	bool		fieldlock;		// info pointer locked?
	const DBFieldSort *
			sortord;		// override sort order
			DBRecord() {
				Empty();
			}
			DBRecord(const DBFieldInfo *fi) {
				Empty(); Init(fi);
			}
			DBRecord(const DBRecord &orig) {
				Empty();
				*this = orig;	// copies memory
			}
			~DBRecord() {
				Free();
			}
			// (re)Initialize record with field info pointer
			// clear previous field values and make the new
			// declarations stick if the usekey boolean is set
	const DBFieldInfo *
			Init(const DBFieldInfo *fi, bool usekey = false) {
				Free();
				if (fi == NULL) { Empty(); return NULL; }
				l = new DBRLink;
				if (usekey | !fieldlock || finfo == NULL) {
					finfo = fi; fieldlock = usekey;
					sortord = NULL;
				}
				if (finfo->GetNFields() > 0)
					sortord = &finfo->sord;
				return finfo;
			}
			// No arguments uninitializes record
	void		Init() {
				Free(); Empty();
			}
			// Is record currently read-only?
	bool		ReadOnly() const {
				return (l == NULL || l->readonly);
			}
			// Set/clear read-only property (propogates across links)
	bool		SetReadOnly(bool set2 = true) {
				if (l == NULL || l->readonly < 0) return set2;
				l->readonly = set2;
				return true;
			}
			// Are we linked (to the given record)?
	int		IsLinked(const DBRecord *rp = NULL) const {
				if (l == NULL) return 0;
				if (rp != NULL && rp->l != l) return 0;
				if (rp == this) return 1;
				return l->nlinks - 1;
			}
			// Establish link with another record (fieldlock unchanged)
	bool		Link(DBRecord *rp) {
				if (rp == NULL) return false;
				if (rp == this) return true;
				Free(); finfo = rp->finfo;
				if (finfo != NULL && finfo->GetNFields() > 0)
					sortord = &finfo->sord;
				else	sortord = NULL;
				return ((l = rp->l) != NULL && l->nlinks++ > 0);
			}
			// Take ownership of another record's contents (exact copy)
	bool		Take(DBRecord *rp) {
				if (rp == NULL) return false;
				if (rp == this) return true;
				Free();
				finfo = rp->finfo; fieldlock = rp->fieldlock;
				l = rp->l; sortord = rp->sortord;
				rp->Empty();
				return (l != NULL && l->nfields > 0);
			}
			// Copy another record's contents without calling Init()
	bool		Copy(const DBRecord &rec) {
				if (ReadOnly()) return false;
				if (finfo != rec.finfo) return false;
				if ((rec.l == NULL) | (rec.l == l)) return false;
				int i = rec.l->nfields;
				l->Init(i);
				while (i-- > 0) l->field[i] = rec.l->field[i];
				l->ihash = rec.HashIndex();
				return (l->nfields > 0);
			}
			// Fill unassigned fields from another record
	int		FillFrom(const DBRecord &rec);
			// Make psuedo-record from a field array we don't own
	bool		Pseudo(const DBFieldInfo *fi, const DBField *fa, int nf) {
				fieldlock = false; Init(fi);
				if ((fi == NULL) | (fa == NULL) | (nf <= 0))
					return false;
				l->field = const_cast<DBField *>(fa);
				l->nfields = nf; l->readonly = -1;
				sortord = &fi->sord;
				return true;
			}
			// Make pseudo-record referencing short field array data
	bool		Pseudo(const DBFieldInfo *fi, const DBFieldS *fa, int nf) {
				if (sizeof(DBFieldS) == sizeof(DBField))
					return Pseudo(fi, (const DBField *)fa, nf);
				fieldlock = false; Init(fi);
				if ((fi == NULL) | (fa == NULL) | (nf <= 0))
					return false;
				l->Init(nf);
				while (nf--) l->field[nf] = &fa[nf];
				l->readonly = 1;
				sortord = &fi->sord;
				return true;
			}
			// Get the field information pointer
	const DBFieldInfo *
			GetFieldInfo() const {
				return finfo;
			}
			// Return maximum fields that can be assigned
	int		GetNFields() const {
				if (finfo == NULL) return 0;
				return finfo->GetNFields();
			}
			// Return upper limit on fields currently in use
	int		GetNAlloc() const {
				if (l == NULL) return 0;
				int	nf = GetNFields();	// paranoia
				if (l->nfields < nf) nf = l->nfields;
				return nf;
			}
			// Return number of fields actually set
	int		GetNAssigned() const {
				int	nset = 0;
				for (int i = GetNAlloc(); i--; )
					nset += (l->field[i].GetNV() > 0);
				return nset;
			}
			// Get field ID from name
	int		FieldID(const char *mnam) const {
				if (finfo == NULL) return -1;
				return finfo->GetID(mnam);
			}
			// Get field name from ID
	const char *	GetFieldName(int fn) const {
				if (finfo == NULL) return NULL;
				return finfo->GetName(fn);
			}
			// Get details for the given field ID
	const DBFDetails *
			GetFieldDetails(int fn) const {
				if (finfo == NULL) return NULL;
				return finfo->GetDetails(fn);
			}
			// Get hash on index fields (DBHVnone if incomplete)
	long		HashIndex() const {
				if (l == NULL) return DBHVnone;
				if (l->ihash >= DBHVnone) return l->ihash;
				return ComputeHashIndex();
			}
			// Get a pointer to an assigned field
	const DBField *	GetField(int fn) const {
				if (l == NULL || (fn < 0) | (fn >= l->nfields) ||
						l->field[fn].GetNV() <= 0)
					return NULL;
				return &l->field[fn];
			}
			// Array operator gets (const) field reference
	const DBField &	operator[](int fn) const {
				static const DBField	Efield;
				if (l == NULL || (fn < 0) | (fn >= l->nfields))
					return Efield;
				return l->field[fn];
			}
			// Format a field value for display
	char *		FormatField(int fn, char *buf, int len) const;
			// Clear a field value
	void		ClearField(int fn) {
				if (ReadOnly()) return;
				if ((fn < 0) | (fn >= l->nfields)) return;
				l->field[fn].Clear();
				if (finfo->GetDetails(fn)->flags & DBFFindex)
					l->ihash = DBHVnone;
			}
			// Clear all fields
	void		Clear() {
				if (!ReadOnly()) l->Init();
			}
			// Typed field setting calls, return true on success
	bool		SetField(int fn, const DBField &src) {
				if (!PrepField(fn, src.GetDT(), src.GetNV()))
					return false;
				l->field[fn] = src; return true;
			}
			// Set corresponding field from source record
	bool		SetField(int fn, const DBRecord &rec) {
				if (rec.finfo == NULL) return false;
				return SetField( fn,
					rec[rec.finfo->XlateID(finfo,fn)] );
			}
	bool		SetField(int fn, int32 lv) {
				if (!PrepField(fn, DBDTint)) return false;
				l->field[fn].Set(lv); return true;
			}
	bool		SetField(int fn, float fv) {
				if (!PrepField(fn, DBDTfloat)) return false;
				l->field[fn].Set(fv); return true;
			}
	bool		SetField(int fn, const char *sv) {
				if (sv == NULL) {
					ClearField(fn); return false;
				}
				if (!PrepField(fn, DBDTstring)) return false;
				l->field[fn].Set(sv); return true;
			}
	bool		SetField(int fn, const int32 ia[], int n);
	bool		SetField(int fn, const float fa[], int n);
	bool		SetField(int fn, const char * const sa[], int n);
	int		SetText(int fn, const char *txt, const char nl = '\n');
			// Write unformatted record to stream
	void		Write(ostream *os, const char tabch = '\t') const {
				for (int i = 0; i < l->nfields; ) {
					*os << l->field[i++];
					if (i < l->nfields) *os << tabch;
				} *os << '\n';
			}
			// Read unformatted record from stream
	int		Read(istream *ins, const char tabch = '\t');
			// Compare against another record using sort info
	int		Compare(const DBRecord &that) const;
			// See if we match a database query
	bool		MatchQuery(const DBQuery *dq) const;
			// Check for exact record match
	bool		operator==(const DBRecord &that) const;
			// Record copy assignment, translate if fieldlock
	DBRecord &	operator=(const DBRecord &src);
			// Union operator (gathers array values)
	DBRecord &	operator|=(const DBRecord &src);
			// Intersection operator (keeps only common values)
	DBRecord &	operator&=(const DBRecord &src);
};

// Compare two records for inequality
inline bool
operator!=(const DBRecord &lr, const DBRecord &rr)
{
	return !(lr == rr);
}

// Record comparison function for qsort()
int		DBrecordCmp(const void *rp1, const void *rp2);

// A convenient class for managing an array of database records
// If initialized with field info, member records are locked to it
class DBRecordList {
	const DBFieldInfo *
			finfo;			// common field information
	DBRecord *	rl;			// pointer to record array
	int		nr;			// number of records
	int		maxr;			// actual array length
	const mutable DBQuery *
			mq;			// query we're matching
	mutable int	spos;			// current search position
	static DBRecord	dummy;			// bogus return record
	void		Empty() {
				rl = NULL; nr = maxr = 0;
				finfo = NULL; spos = 0;
				sord.Clear();
			}
public:
	DBFieldSort	sord;			// overriding sort order
			DBRecordList(const DBFieldInfo *fi = NULL, int n = 0) {
				Empty();
				Init(fi, n);
			}
			DBRecordList(const DBRecordList &orig) {
				Empty();
				*this = orig;		// copies records
			}
			~DBRecordList() {
				Resize(0);
			}
			// (re)initialize
	int		Init(const DBFieldInfo *fi = NULL, int n = 0) {
				Resize(0);
				if ((finfo = fi) != NULL) sord = fi->sord;
				else sord.Clear();
				return Resize(n);
			}
			// Get field info pointer
	const DBFieldInfo *
			GetFieldInfo() const {
				return finfo;
			}
			// Get array size
	int		GetSize() const {
				return nr;
			}
			// Array pointer access
	DBRecord *	Array() {
				return rl;
			}
			// Const array pointer
	const DBRecord *
			GetArray() const {
				return rl;
			}
			// Array accessor
	DBRecord &	operator[](int i) {
				if ((i < 0) | (i >= nr)) return dummy;
				return rl[i];
			}
			// Const accessor
	const DBRecord &
			Get(int i) const {
				if ((i < 0) | (i >= nr)) return dummy;
				return rl[i];
			}
			// Find the given record (exact match)
	int		FindRecord(const DBRecord &dr) const {
				if (nr <= 0) return -1;
				int	i;
				i = (const char *)&dr - (const char *)rl;
				if (i >= 0 && i < nr*(int)sizeof(DBRecord))
					return i/sizeof(DBRecord);
				for (i = nr; i--; )
					if (dr == rl[i]) return i;
				return -1;
			}
			// Find first record matching query
	int		FindFirst(const DBQuery *dq) const {
				if ((mq = dq) == NULL) return -1;
				spos = 0; return FindNext();
			}
			// Find next query match (or -1 if no more)
	int		FindNext() const {
				for ( ; spos < nr; spos++)
					if (rl[spos].MatchQuery(mq))
						return spos++;
				mq = NULL; return -1;
			}
			// Resize array, returning new size
	int		Resize(int n);
			// Allocate new record at end of list
	DBRecord &	NewRecord() {
				if (!Resize(nr+1)) return dummy;
				return rl[nr-1];
			}
			// Eliminate uninitialized or empty records
	void		Compact(bool delEmpty = false);
			// Clear (portion of) array
	void		Clear(int i = 0, int end = 0) {
				if (i < 0) i = 0;
				if (i >= nr) return;
				if ((end <= 0) | (end > nr)) end = nr;
				while (i < end)
					rl[i++].Init(finfo, true);
			}
			// Sort records in array
	void		Sort();
			// Write records to stream
	void		Write(ostream *os, const char tabch = '\t') const {
				if (os == NULL) return;
				if (nr <= 0) return;
				if (finfo != NULL) finfo->Write(os, tabch);
				else rl[0].GetFieldInfo()->Write(os, tabch);
				for (int i = 0; i < nr; i++)
					rl[i].Write(os, tabch);
			}
			// Read records from stream, adding to our list
	int		Read(istream *ins, char *tabcp = NULL,
					DBFieldInfo *fi = NULL);
			// Link all records from another list
	int		LinkAll(DBRecordList *srcp);
			// Compute union with all records
	void		Union(DBRecord *drp, bool add2 = false) const {
				if (!add2) drp->Init(finfo);
				int i = nr;
				while (i-- > 0) *drp |= rl[i];
			}
			// Compute intersection with all records
	void		Intersection(DBRecord *drp, bool add2 = false) const {
				if (!add2) drp->Init(finfo);
				int i = nr;
				while (i-- > 0) *drp &= rl[i];
			}
			// Record array copy assignment
	DBRecordList &	operator=(const DBRecordList &src);
};

// Return prime number reasonably close to size for hashing
inline uint32
primeAbove(uint32 siz)
{
	const uint32  primetab[] = {
		31, 61, 127, 251, 509, 1021, 2039, 4093, 8191, 16381, 
		32749, 65521, 131071, 262139, 524287, 1048573, 2097143, 
		4194301, 8388593, 0
	};
	for (int i = 0; primetab[i] > 0; i++)
		if (primetab[i] >= siz)
			return primetab[i];
	return siz;
}

extern const uint32	DB_HASHSETLEN;		// index hash set size (prime)

// Database Page States:  clean, needs sync, needs compacting, error lock
enum	DBPageState	{ DBPSclean, DBPSsync, DBPScompact, DBPSerror };

// Class for accessing database record blocks
// If opened read-only, changes cannot be written to file
class DBAccess : CacheObject {
private:
	off_t		strpos;			// page offset in stream
	int		rwidth;			// record width (fields)
	int		nrused;			// number of records used
	int		nravail;		// number available (approx.)
	DBAccess *	next;			// next page in database
	ABitMap		bihset;			// block index hash set
	DBPageState	state;			// resident page state
	DBRecordList	brl;			// block pseudo record list
	virtual bool	RestoreMemory();
	virtual void	FreeMemory();
	void		VirginBlock() {
				rwidth = 1;
				nrused = 0;
				nravail = DB_FBLEN/rwidth;
				bihset.ClearBitMap(false);
			}
	void		UnknownBlock() {
				rwidth = -1;
				nrused = -1;
				nravail = -1;
				bihset.ClearBitMap(true);
			}
protected:
	const DBHeader *
			header;			// DB field information
	void		FreeLoad() {
				ReleaseObject(GetCacheObject());
			}
	int		RecordWidth() {
				if (rwidth < 0)
					FreeLoad();
				else if (!nrused && header != NULL)
					rwidth = header->GetNFields();
				return rwidth;
			}
	int		NRecords() {
				if (nrused < 0) FreeLoad();
				return nrused;
			}
	int		NRecordsAvail() {
				if (nravail < 0) FreeLoad();
				return nravail;
			}
	bool		HashReject(long ih) const {
				if (ih < 0) return false;
				return !bihset.Check(ih % DB_HASHSETLEN);
			}
	void		ListRecords();
	void		IndexRecords();
	bool		PutRecord(const DBRecord &rec);
	bool		SyncPage(bool keep = true);
			DBAccess(const DBAccess *prev, off_t dblen = 0);
public:
			DBAccess(DBHeader *hdr = NULL, iostream *dbstrm = NULL,
					bool ro = false) {
				header = NULL; next = NULL;
				Init(hdr, dbstrm, ro);
			}
	virtual		~DBAccess() {
				delete next;
				while (InUse()) ReleaseObject(objCache);
				FreeMemory();
			}
			// Get name of subclass
	virtual const char *
			ClassName() const {
				return "DBAccess";
			}
			// Initialize database with given header/stream
	bool		Init(DBHeader *hdr, iostream *dbstrm = NULL,
					bool ro = false);
			// Is database read-only?
	bool		ReadOnly() const {
				if (header == NULL) return true;
				return header->ReadOnly();
			}
			// Get database header pointer
	const DBHeader *
			GetHeader() const {
				return header;
			}
			// Is the database ready for reading?
	bool		Ready() const {
				if (state == DBPSerror) return false;
				if (next != NULL) return next->Ready();
				if (header == NULL) return false;
				if (header->GetStream() == NULL) return false;
				return (header->GetSize() > 0);
			}
			// Execute callback for each database block
	virtual int     ForEachBlock(int (*f)(const DBRecordList &, void *),
					void *client_data = NULL,
					const DBFieldSort *ord = NULL);
			// Execute callback for each database record
	int		ForEachRecord(int (*f)(const DBRecord &, void *),
					void *client_data = NULL,
					const DBFieldSort *ord = NULL);
			// Get n records starting at rn, return #found
	virtual int	GetRecords(DBRecord *rlist, int n, int rn = 0);
			// Find n records matching query starting from *rnp
	virtual int	FindRecords(const DBQuery *dq, DBRecord *rlist, int n,
					int *rnp = NULL);
			// Add n records to database
	virtual int	AddRecords(const DBRecord *rlist, int n);
			// Add a list of records to database
	int		AddRecordList(const DBRecordList &rls) {
				return AddRecords(rls.GetArray(),rls.GetSize());
			}
			// Add a single record to database
	bool		AddRecord(const DBRecord &r) {
				return (bool)AddRecords(&r, 1);
			}
			// Delete all matching records from database
	virtual int	DeleteRecords(const DBRecord *rlist, int n);
			// Delete a list of records from database
	int		DeleteRecordList(const DBRecordList &rls) {
				return DeleteRecords(rls.GetArray(),rls.GetSize());
			}
			// Delete records matching this one from database
	bool		DeleteRecord(const DBRecord &r) {
				return (bool)DeleteRecords(&r, 1);
			}
			// Return the total records in database
	virtual int	TotalRecords() {
				int	mr = 0, nr = NRecords();
				if ((nr < 0) | (state == DBPSerror)) return -1;
				if (next != NULL) {
					mr = next->TotalRecords();
					if (mr < 0) return -1;
				}
				return nr + mr;
			}
			// Write all database records to a stream
	virtual void	Write(ostream *os, const char tabch = '\t');
			// Write all records (sorted)
	void		WriteSorted(ostream *os, const char tabch = '\t');
			// Read records from stream, adding to our DB
	virtual int	Read(istream *ins, char *tabcp = NULL,
					DBFieldInfo *fi = NULL);
			// Synchronize in-memory data to disk
	virtual bool	Sync() {
				bool	aok = SyncPage();
				if (next != NULL) aok &= next->Sync();
				else if (header != NULL) aok &= header->Sync();
				return aok;
			}
			// Write out database statistics
	void		PrintStats(ostream *os);
			// Close database object, return true if OK
	virtual bool	Close();
};

/*
 * Database query selector list
 * Fields set in the same selector are logically AND'ed together.
 * Array field values are also AND'ed together, and a match
 * in any element of a queried field is considered a match,
 * e.g., to match a record with multiple keywords.
 * If matchEmpty is true, then unassigned _queried_ fields will match.
 * (Unassigned fields in the query itself match regardless.)
 * Subsequent queries in selector list are logically OR'ed with this one.
 */
struct DBQuery {
	DBRecord	rlo, rhi;		// lower and upper limits
	bool		matchEmpty;		// match unassigned fields
	DBQuery *	next;			// next selector
	const ABitMap *	hashSet;		// hash on hash index set
			DBQuery(const DBFieldInfo *fi, bool me = false) {
				matchEmpty = me; hashSet = NULL; next = NULL;
				rlo.Init(fi, true);
				rhi.Init(fi, true);
			}
			DBQuery(DBQuery *nq = NULL) {
				matchEmpty = false; hashSet = NULL;
				if ((next = nq) == NULL) return;
				matchEmpty = nq->matchEmpty;
				rlo.Init(nq->rlo.GetFieldInfo(), true);
				rhi.Init(nq->rhi.GetFieldInfo(), true);
			}
			~DBQuery() {
				delete next;
			}
			// Get related hash if rlo & rhi match index fields
	long		HashIndex() const {
				const long	hlo = rlo.HashIndex();
				if (hlo < 0) return hlo;
				if (rhi.HashIndex() == hlo) return hlo;
				return DBHVnone;
			}
			// Assign hash index set for the entire query list
	bool		HashSet(ABitMap *hs);
			// Is the given index ruled out by our hash set?
	bool		HashReject(long ih) const {
				if (hashSet == NULL) return false;
				const uint32	hlen = hashSet->Length();
				if (!hlen) return false;
				if (ih < 0) return !matchEmpty;
				return !hashSet->Check(ih % hlen);
			}
			// Return first possible query to match the given hash
	const DBQuery *	FirstPossible(long ih) const;
			// Return next possible match
	const DBQuery *	NextPossible(long ih) const {
				if ((ih < 0) & matchEmpty) return next;
				const DBQuery *	dq = next;
				while (dq != NULL) {
					const long	th = dq->HashIndex();
					if ((th < 0) | (th == ih)) break;
					dq = dq->next;
				}
				return dq;
			}
};

// Class for performing database searches
// To specify a search, set the relevant record fields and call AddSelector
class DBSearch {
	DBAccess *	dbacc;			// database access object
	DBQuery *	qlist;			// query list
	ABitMap		ihSet;			// index hash set (optimization)
	int		crpos;			// current record position
public:
			DBSearch(DBAccess *dba = NULL) {
				qlist = NULL;
				SetDB(dba);
			}
			~DBSearch() {
				Clear();
			}
			// Set the database to search
	bool		SetDB(DBAccess *dba) {
				Clear();
				return ((dbacc = dba) != NULL &&
						dbacc->GetHeader() != NULL);
			}
			// Clear the current selection
	void		Clear() {
				delete qlist; qlist = NULL; crpos = 0;
			}
			// Get the current query list
	const DBQuery *	GetQuery() const {
				return qlist;
			}
			// Add a selection range, OR'ing with previous
	bool		AddSelector(const DBRecord &r1, const DBRecord &r2,
					bool me = false);
			// Add a selection value, OR'ing with previous
	bool		AddSelector(const DBRecord &req, bool me = false) {
				return AddSelector(req, req, me);
			}
			// Reset search to beginning
	void		Reset() {
				if (crpos < 0) qlist->HashSet(&ihSet);
				crpos = 0;
			}
			// Get next n records matching search criteria
	int		FindRecords(DBRecord *rlist, int n = 1) {
				if (n <= 0) return 0;
				if (dbacc == NULL) return 0;
				if (crpos < 0) Reset();
				return dbacc->FindRecords(qlist, rlist, n,
						&crpos);
			}
			// Find all records matching search criteria
	int		FindAll(DBRecordList *rl);
			// Delete all records matching search criteria
	int		DeleteMatching();
};

#endif	// ! _DBASE_H_
