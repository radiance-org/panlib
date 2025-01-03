/*
 *  dheader.cpp
 *  panlib
 *
 *  Database header class implementations
 *
 *  Created by gward on June 12, 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#include <ctype.h>
#include "system.h"
#include "rtio.h"
#include "dbase.h"
#include "dmessage.h"

// Check machine byte order
static inline bool
CheckBigendian()
{
	union { int32 i; char c[4]; } u;
	u.i = 1;
	return (u.c[0] == 0);
}

const bool	DBHeader::bigendian = CheckBigendian();

// Create new field ID, returning -1 on error
int
DBFieldID::NewField(const char *mnam)
{
	if (mnam == NULL || !*mnam)
		return -1;
	if (nfields >= DB_MAXFIELD)		// out of fields?
		return -1;
	mname[nfields] = strlink(mnam);
	if ((mname[nfields] == NULL) |		// out of memory?
			(mname[nfields] == strNoMemory))
		return -1;
	int		c = 1;
	int		n;			// insert in sorted list
	for (n = nfields; n-- && (c = istrcmp(mnam, mname[mlexo[n]])) < 0; )
		mlexo[n+1] = mlexo[n];
	if (!c)
		return -1;			// duplicate field
	return (mlexo[n+1] = nfields++);	// return new ID
}

// Rename the given field
bool
DBFieldID::Rename(int mid, const char *mnam)
{
	if ((mid < 0) | (mid >= nfields))
		return false;
	if (mnam == NULL || !*mnam)
		return false;
	if (!istrcmp(mnam, mname[mid]))
		return true;			// no change
	if (GetID(mnam) >= 0)
		return false;			// name is taken
	strunlink(mname[mid]);
	mname[mid] = strlink(mnam);
	if ((mname[mid] == NULL) |		// out of memory?
			(mname[mid] == strNoMemory))
		return false;
	int		n;			// fix sort order
	for (n = 0; mlexo[n] != mid; ++n)
		;
	while (n < nfields-1 && istrcmp(mnam, mname[mlexo[n+1]]) > 0) {
		mlexo[n] = mlexo[n+1];
		++n;
	}
	while (n > 0 && istrcmp(mnam, mname[mlexo[n-1]]) < 0) {
		mlexo[n] = mlexo[n-1];
		--n;
	}
	mlexo[n] = mid;
	return true;
}

// Find field ID by name
int
DBFieldID::GetID(const char *mnam) const
{
	if (mnam == NULL)
		return -1;
	int	i;
	if (nfields < 12) {			// linear search
		for (i = 0; i < nfields; i++)
			if (!istrcmp(mnam, mname[i]))
				return i;
		return -1;
	}
	int		ilower = 0, iupper = nfields;
	int		c = iupper;		// binary search
	while ((i = (iupper + ilower) >> 1) != c) {
		c = istrcmp(mnam, mname[mlexo[i]]);
		if (c > 0)
			ilower = i;
		else if (c < 0)
			iupper = i;
		else
			return mlexo[i];
		c = i;
	}
	return -1;
}

// Map IDs for foreign field source, returning false if 1-1
bool
DBFieldID::MapIDs(int idmap[DB_MAXFIELD], const DBFieldID *srcp) const
{
	if (srcp == NULL)			// no fields to map
		return false;
	int		i = srcp->nfields;
	if (this == NULL) {			// nothing to map to!
		while (i--)
			idmap[i] = -1;
		return true;
	}
	if (this == srcp) {			// definite no-op
		while (i--)
			idmap[i] = i;
		return false;
	}
	bool		needmap = (srcp->nfields > nfields);
	int		j;
	for (i = j = 0; i < srcp->nfields; i++) {
		int	cm = 0;
		while (j < nfields &&
				(cm = istrcmp(srcp->mname[srcp->mlexo[i]],
					mname[mlexo[j]])) > 0) {
			needmap = true;
			j++;
		}
		if ((j >= nfields) | (cm != 0)) {
			needmap = true;
			idmap[srcp->mlexo[i]] = -1;
		} else
			idmap[srcp->mlexo[i]] = mlexo[j++];
	}
	return needmap;
}

static inline bool
isFieldName(int c)
{
	return (isalnum(c) || c == '_');
}

// Read delimiter-separated field list from stream
int
DBFieldID::Read(istream *ins, char *tabcp)
{
	if (ins == NULL) return -1;
	char		linebuf[DB_MNAMELEN*DB_MAXFIELD];
	char *		cp = linebuf;
	Clear();				// clear field IDs
	while (ins->good()) {			// read header line
		int	c = ins->get();
		if ((c == '\n') | (c == '\r') | (c == EOF))
			break;
		*cp++ = c;
		if (cp >= linebuf+sizeof(linebuf))
			return -1;
	}
	*cp = '\0';
	char		tabch = '\0';		// default delimiter (any)
	const char *	nm = cp = linebuf;
	if (tabcp == NULL)
		tabcp = &tabch;
	while (*cp) {				// add each field
		if (!*tabcp && !isFieldName(*cp)) {
			if ((*cp == ',') | (*cp == '"') | (*cp == '\\') |
					(*cp == '{') | (*cp == '}') |
					(*cp == '(') | (*cp == ')'))
				return -1;
			*tabcp = *cp;
		}
		if (*cp == *tabcp) {
			*cp++ = '\0';
			if (*nm && NewField(nm) < 0)
				return -1;
			nm = cp;
		} else
			cp++;
	}
						// final field
	if (*nm && NewField(nm) < 0)
		return -1;
	return nfields;
}

// Memberwise copy operator
DBFieldID &
DBFieldID::operator=(const DBFieldID &src)
{
	if (this == &src)
		return *this;
	Clear();
	while (nfields < src.nfields) {
		mname[nfields] = strlink(src.mname[nfields]);
		mlexo[nfields] = src.mlexo[nfields];
		++nfields;
	}
	return *this;
}

// Format integer field
char *
DBformatInteger(char *buf, int len, const DBField &f)
{
	if ((buf == NULL) | (len < 2))
		return NULL;
	int32	ival;
	if (!f.Get(&ival)) {
		buf[0] = '\0';
		return NULL;
	}
	char	mybuf[16];
	sprintf(mybuf, "%d", ival);
	int	rlen = strlen(mybuf);
	if (rlen <= len)
		strcpy(buf, mybuf);
	else if (ival > 0)
		strcpy(buf, "+");
	else
		strcpy(buf, "-");
	return buf;
}

// Format real field
char *
DBformatFloat(char *buf, int len, const DBField &f)
{
	if ((buf == NULL) | (len < 2))
		return NULL;
	float	fval;
	if (!f.Get(&fval)) {
		buf[0] = '\0';
		return NULL;
	}
	char	mybuf[32];
	int	rlen;
	int	prec = 8;
	do {
		char	fmt[8];
		sprintf(fmt, "%%.%dg", prec);
		sprintf(mybuf, fmt, fval);
		rlen = strlen(mybuf);
	} while (rlen > len && prec-- > 0);
	if (rlen <= len)
		strcpy(buf, mybuf);
	else if (fval > 0)
		strcpy(buf, "+");
	else
		strcpy(buf, "-");
	return buf;
}

// Format string field
char *
DBformatString(char *buf, int len, const DBField &f)
{
	if ((buf == NULL) | (len < 2))
		return NULL;
	const char *	sval;
	if (!f.Get(&sval)) {
		buf[0] = '\0';
		return NULL;
	}
	strlcpy(buf, sval, len);
	return buf;
}

// Format unknown (any) field
char *
DBformatUnknown(char *buf, int len, const DBField &f)
{
	if ((buf == NULL) | (len < 2) | !f.GetNV())
		return NULL;
	switch (f.GetDT()) {
	case DBDTint:
		return DBformatInteger(buf, len, f);
	case DBDTfloat:
		return DBformatFloat(buf, len, f);
	case DBDTstring:
		return DBformatString(buf, len, f);
	default:
		strlcpy(buf, "BADTYPE", len);
		break;
	}
	return NULL;
}

// Format (any) array field
char *
DBformatArray(char *buf, int len, const DBField &f)
{
	if (f.GetNV() < 2 && f.GetDT() != DBDTstring)
		return DBformatUnknown(buf, len, f);
	char *	cp = buf;
#define	rem (len - (cp-buf))
	*cp = '\0';
	for (int i = 0; i < f.GetNV(); i++) {
		const bool	more = (i < f.GetNV()-1);
		DBField		f1;
		if (!f.Get(&f1, i))
			return NULL;
		if (f1.GetDT() == DBDTstring) {
			const char	quot = strchr(f1.GetString(),'"')==NULL
						? '"' : '\'';
			*cp++ = quot;
			if (DBformatString(cp, rem, f1) == NULL)
				return NULL;
			while (*cp) ++cp;
			if (rem < 2 + 5*more)
				break;
			*cp++ = quot; *cp = '\0';
		} else {
			if (DBformatUnknown(cp, rem, f1) == NULL)
				return NULL;
			while (*cp) ++cp;
			if (rem < 1 + 3*more)
				break;
		}
		if (more) {
			*cp++ = ','; *cp++ = ' '; *cp = '\0';
		}
	}
	return buf;
#undef rem
}

// Add field to end of sort priority list
bool
DBFieldSort::AddSort(int mid, bool reverse)
{
	if (mid < 0)
		return false;
	++mid;			// list stores ID+1
	int	i;
	for (i = 0; i < DB_MAXSORTF && order[i]; i++)
		if ((order[i] == mid) | (order[i] == -mid))
			return false;
	if (i >= DB_MAXSORTF)
		return false;
	order[i++] = reverse ? -mid : mid;
	order[i] = 0;		// terminator
	return true;
}

// Push field to top of sort priority list
bool
DBFieldSort::PushSort(int mid, bool reverse)
{
	if (mid < 0)
		return false;
	++mid;			// list stores ID+1
	int	i;
	for (i = 0; i < DB_MAXSORTF; i++) {
		if (!order[i]) {
			++i; break;
		}
		if ((order[i] == mid) | (order[i] == -mid))
			break;
	}
	while (i--)		// bump previous sort
		order[i+1] = order[i];
	order[0] = reverse ? -mid : mid;
	return true;
}

// Define new field information and return ID
int
DBFieldInfo::Field(const char *mnam, DBDataType dtype,
			const char *descrip, int flags,
			char *(*df)(char *, int, const DBField &))
{
	int	mid = DBFieldID::Field(mnam);
	if (mid < 0)
		return -1;
						// check changes
	if (descrip == NULL)
		descrip = mnam;
	if (df == NULL && dtype == details[mid].dtype &&
			(flags & DBFFarray) == (details[mid].flags & DBFFarray))
		df = details[mid].DisplayFormat;
	if (df == NULL) {
		switch (dtype) {
		case DBDTunknown:	df = &DBformatUnknown;	break;
		case DBDTint:		df = &DBformatInteger;	break;
		case DBDTfloat:		df = &DBformatFloat;	break;
		case DBDTstring:	df = &DBformatString;	break;
		}
		if (flags & DBFFarray) df = &DBformatArray;
	}
						// copy details
	details[mid].dtype = dtype;
	strlcpy(details[mid].descrip, descrip, DB_MAXDESCR);
	details[mid].flags = flags & DBFFinfoM;
	details[mid].DisplayFormat = df;
	return mid;
}

// Copy corresponding formatting functions
void
DBFieldInfo::CopyFormat(const DBFieldInfo *fi)
{
	if (fi == NULL || fi->nfields <= 0)
		return;
	int	idmap[DB_MAXFIELD];
	int	i;
	MapIDs(idmap, fi);
	for (i = fi->nfields; i--; ) {
		if (idmap[i] < 0)
			continue;
		if (fi->details[i].dtype != details[idmap[i]].dtype)
			continue;
		details[idmap[i]].DisplayFormat = fi->details[i].DisplayFormat;
	}
}

// Memberwise copy operator for field information
DBFieldInfo &
DBFieldInfo::operator=(const DBFieldInfo &src)
{
	if (this == &src)
		return *this;
	DBFieldID::operator=(src);
	memcpy(details, &src.details, sizeof(details));
	sord = src.sord;
	return *this;
}

// Alter the following at your peril...
static const char	DBMagic[] = "ASDB = {\n";
static const char	DBSoftwareName[] = "Software";
static const char	DBFieldSizeName[] = "FieldSize";
static const char	DBFieldsPerBlockName[] = "FieldsPerBlock";
static const char	DBBigEndianName[] = "BigEndian";
static const char	DBHeaderSizeName[] = "HeaderSize";
static const char	DBFieldInfoName[] = "FieldInfo";
static const char	DBSortOrderName[] = "SortOrder";
static const char	DBdtypeName[4][8] = {"unknown", "int",
						"float", "string"};

// Write field information to stream
int
DBFieldInfo::WriteInfo(ostream *ostr, int tab) const
{
	int	i;
	TabIn(ostr, tab);
	*ostr << DBSortOrderName << " =";
	for (i = 0; sord.order[i]; i++)
		*ostr << ' ' << sord.order[i];
	*ostr << " 0\n";
	for (i = 0; i < nfields; i++) {
		TabIn(ostr, tab);
		*ostr << '"' << mname[i] << "\",\t";
		*ostr << details[i].flags;
		*ostr << ",\t" << DBdtypeName[details[i].dtype];
		*ostr << ",\t\"" << details[i].descrip << "\"\n";
	}
	return ostr->bad() ? -1 : nfields;
}

// Read field information from stream up until line with '}'
int
DBFieldInfo::ReadInfo(istream *istr)
{
	Clear();				// start fresh
	const char *	err = "Unexpected EOF";
	int		c;
	DBHeaderLink	sortlink(DBSortOrderName, sord.order, DB_MAXSORTF+1);
	if (sortlink.ScanVar(istr) != &sortlink) {
		err = "Expected sort order";
		goto erreturn;
	}
	while ((c = istr->peek()) != EOF) {
		if (c == '}') {			// that's all folks...
			istr->ignore(64, '\n');
			return nfields;
		}
		if (isspace(c)) {		// skip leading space
			istr->ignore();
			continue;
		}
		if (c == '/' || c == '#') {	// skip comment line
			istr->ignore(256, '\n');
			continue;
		}
		if (c != '"') {
			err = "Expected quote for field name";
			break;
		}
		istr->ignore();			// skip starting quote
		char	mnam[64];
		if (istr->get(mnam, sizeof(mnam), '"').fail()) {
			err = "Missing field name";
			break;
		}
		istr->ignore(2);		// skip end quote and comma
		if (GetID(mnam) >= 0) {
			err = "Duplicate field";
			break;
		}
		int	fl;
		*istr >> fl;			// get flags
		if (istr->fail()) {
			err = "Missing flags";
			break;
		}
		istr->ignore();			// skip comma
		while (isspace(istr->peek()))
			istr->ignore();		// skip leading space
		char		buf[DB_MAXDESCR];
		DBDataType	dt;		// get data type
		if (istr->get(buf, sizeof(buf), ',').fail()) {
			err = "Missing data type";
			break;
		}
		if (!strcmp(buf, DBdtypeName[DBDTint]))
			dt = DBDTint;
		else if (!strcmp(buf, DBdtypeName[DBDTfloat]))
			dt = DBDTfloat;
		else if (!strcmp(buf, DBdtypeName[DBDTstring]))
			dt = DBDTstring;
		else if (!strcmp(buf, DBdtypeName[DBDTunknown]))
			dt = DBDTunknown;
		else {
			err = "Unrecognized data type";
			break;
		}
		if (!istr->ignore(7, '"').good()) {
			err = "Expected quote for field description";
			break;
		}
		if (istr->get(buf, sizeof(buf), '"').fail()) {
			err = "Cannot read field description";
			break;
		}
		if (!istr->ignore(7, '\n').good()) {
			err = "Syntax error after description";
			break;
		}
		if (Field(mnam, dt, buf, fl) < 0) {
			err = "Field definition failed";
			break;
		}
	}
erreturn:
	DMESGF(DMCdata, "%s in DBFieldInfo::ReadInfo", err);
	return -1;
}

static int
DBHWriteFieldInfo(void *ho, ostream *ostr, int tab)
{
	return (*(DBHeader *)ho).DBFieldInfo::WriteInfo(ostr, tab);
}

static int
DBHReadFieldInfo(void *ho, istream *istr)
{
	return (*(DBHeader *)ho).DBFieldInfo::ReadInfo(istr);
}

// Initialize database header object
DBHeader::DBHeader()
{
	strm = NULL; SetChangedFlag(true); readonly = true;
	hvars = NULL;				// set up variable links
	soft[0] = '\0';
	hvars = new DBHeaderLink(DBSoftwareName, soft, sizeof(soft), hvars);
	fieldsize = sizeof(DBFieldS);
	hvars = new DBHeaderLink(DBFieldSizeName, &fieldsize, 1, hvars);
	fieldsperblock = DB_FBLEN;
	hvars = new DBHeaderLink(DBFieldsPerBlockName, &fieldsperblock, 1, hvars);
	bigendfile = bigendian;
	hvars = new DBHeaderLink(DBBigEndianName, &bigendfile, 1, hvars);
	hdrsiz = InfoSize();
	hvars = new DBHeaderLink(DBHeaderSizeName, &hdrsiz, 1, hvars);
	hvars = new DBHeaderLink(DBFieldInfoName, this,
			&DBHWriteFieldInfo, &DBHReadFieldInfo,
			DBFieldInfo::InfoSize(), hvars);
}

// Assign default header values
void
DBHeader::SetDefaults()
{
	hvars->Clear();				// clear linked values
	hvars->SetVal(soft, "(unknown)");
	hvars->SetVal(&fieldsize, sizeof(DBFieldS));
	hvars->SetVal(&fieldsperblock, DB_FBLEN);
	hvars->SetVal(&bigendfile, (int)bigendian);
	hvars->SetVal(&hdrsiz, InfoSize());
	Clear();				// clear fields & lie
	hvars->FindVar(DBFieldInfoName)->nset = 1;
}

// Check current assignments
bool
DBHeader::CheckValues() const
{
	if (fieldsize != sizeof(DBFieldS))
		return false;			// can't handle this
	if (fieldsperblock != DB_FBLEN)
		return false;			// blocksize must also match
	return true;				// otherwise, OK
}

// Position at beginning of file and write magic number
bool
DBHeader::PutMagic() const
{
	strm->clear();
	strm->seekp(0);
	*strm << DBMagic;
	return !strm->bad();
}

// Check magic number
bool
DBHeader::CheckMagic() const
{
	strm->clear();
	strm->seekg(0);
	return MatchInput(strm, DBMagic);
}

// Attach database header to open stream
int
DBHeader::Attach(iostream *dbstream, bool ro)
{
	Sync();					// sync old stream
	hvars->Clear();				// clear linked values
	Clear();				// clear fields
	if ((strm = dbstream) == NULL)
		return hdrsiz = 0;
	if (!CheckMagic()) {
		strm->clear();
		if (strm->tellg() > 0)		// bad magic?
			return hdrsiz = 0;
		if (ro)				// empty file
			return hdrsiz = 0;
		SetDefaults();			// start fresh
		SetChangedFlag(true);
		readonly = false;
		return hdrsiz;
	}
	if (ReadInfo(strm) < 0)			// header OK?
		return hdrsiz = 0;
	SetChangedFlag(false);			// got it
	readonly = ro;
	return hdrsiz;
}

// Synchronize header with stream
bool
DBHeader::Sync() const
{
	if ((strm == NULL) | (hdrsiz <= 0))
		return false;
	int32	newhash = Hash();		// check for changes
	if (lastsync == newhash)		// !HasChanged()
		return true;			// same hash ==> assume OK
	if (readonly)
		return false;			// altered read-only header
	int	spos;
	if (!PutMagic() || WriteInfo(strm) < 0)	// else start writing
		goto writerr;
	*strm << "}\t// ASDB\n";
	spos = (int)strm->tellp();
	if (spos > hdrsiz) {
		DMESG(DMCparameter, "Database header overflowed");
		return false;			// OOPS!
	}
	while (spos++ < hdrsiz)
		strm->put('\n');		// top off with newlines
	if (!strm->flush().bad()) {
		lastsync = newhash;		// SetChangedFlag(false)
		return true;
	}
writerr:
	DMESG(DMCsystem, "Database header write error");
	return false;
}
