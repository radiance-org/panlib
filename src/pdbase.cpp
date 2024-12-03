/*
 *  pdbase.cpp
 *  panlib
 *
 *  Pancine database routines.
 *
 *  Created by gward on Wed Jun 06 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#include <math.h>
#include "pancine.h"

#ifndef IO_CHECKLIM
#define IO_CHECKLIM	10000		// maximum records to double-check
#endif

// Software name as it appears in image database header
const char		PancineSoftwareName[] = "Pancine Photosphere";

// Do not alter the following header variable names
const char		PDBUniqBoolName[] = "UniqueFiles";
const char		PDBLastSubjName[] = "RecentSubject";
const char		PDBLastAlbmName[] = "RecentAlbum";
const char		PDBLastOwnrName[] = "RecentOwner";
const char		PDBLastKeywName[] = "RecentKeyword";

const PDBFieldInfo	PDBFInfo;		// Pancine field information

const PDBFieldInfo	PDBSFInfo;		// Pancine search field info.

// Initialize Pancine field information (order must match enum PDBFieldID)
PDBFieldInfo::PDBFieldInfo()
{
	static const char	eemsg[] = "enumeration error in PDBFieldInfo";

	if (Field("File", DBDTstring, "File Name",
			DBFFshow|DBFFrequired|DBFFindex) != PDBFfile)
		DMESG(DMCassert, eemsg);
	if (Field("Location", DBDTstring, "File Location",
			DBFFarray|DBFFrequired, &PDBformatFolder) != PDBFlocation)
		DMESG(DMCassert, eemsg);
	if (Field("NBytes", DBDTint, "File Length (bytes)",
			DBFFrequired|DBFFindex) != PDBFnbytes)
		DMESG(DMCassert, eemsg);
	if (Field("HashVal", DBDTint, "File Hash Value",
			DBFFhide|DBFFrequired|DBFFindex) != PDBFhashval)
		DMESG(DMCassert, eemsg);
	if (Field("Format", DBDTstring, "File Format",
			0) != PDBFformat)
		DMESG(DMCassert, eemsg);
	if (Field("ModTime", DBDTint, "File Modified Time",
			DBFFhide) != PDBFmodtime)
		DMESG(DMCassert, eemsg);
	if (Field("CapDate", DBDTstring, "Capture Date and Time",
			DBFFshow, &PDBformatDate) != PDBFcapdate)
		DMESG(DMCassert, eemsg);
	if (Field("GMTdate", DBDTstring, "Precise GMT at Capture",
			0, &PDBformatDate) != PDBFgmtdate)
		DMESG(DMCassert, eemsg);
	if (Field("Latitude", DBDTfloat, "Degrees North Latitude",
			DBFFuser) != PDBFlatitude)
		DMESG(DMCassert, eemsg);
	if (Field("Longitude", DBDTfloat, "Degrees East Longitude",
			DBFFuser) != PDBFlongitude)
		DMESG(DMCassert, eemsg);
	if (Field("Altitude", DBDTfloat, "Meters above Sea Level",
			DBFFuser) != PDBFaltitude)
		DMESG(DMCassert, eemsg);
	if (Field("Bearing", DBDTfloat, "Sight Bearing degrees East of North",
			DBFFuser) != PDBFbearing)
		DMESG(DMCassert, eemsg);
	if (Field("Subject", DBDTstring, "General Subject",
			DBFFuser|DBFFshow|DBFFrequired) != PDBFsubject)
		DMESG(DMCassert, eemsg);
	if (Field("Title", DBDTstring, "Title",
			DBFFuser|DBFFshow) != PDBFtitle)
		DMESG(DMCassert, eemsg);
	if (Field("Owner", DBDTstring, "Content Owner",
			DBFFuser) != PDBFowner)
		DMESG(DMCassert, eemsg);
	if (Field("History", DBDTstring, "File History",
			DBFFarray) != PDBFhistory)
		DMESG(DMCassert, eemsg);
	if (Field("Album", DBDTstring, "Album Containing File",
			DBFFuser|DBFFarray|DBFFshow) != PDBFalbum)
		DMESG(DMCassert, eemsg);
	if (Field("Position", DBDTint, "Position in Album",
			DBFFuser|DBFFarray|DBFFhide) != PDBFposition)
		DMESG(DMCassert, eemsg);
	if (Field("Caption", DBDTstring, "Caption in Album",
			DBFFuser|DBFFarray) != PDBFcaption)
		DMESG(DMCassert, eemsg);
	if (Field("Keyword", DBDTstring, "Keyword for Search",
			DBFFuser|DBFFarray|DBFFshow) != PDBFkeyword)
		DMESG(DMCassert, eemsg);
	if (Field("Comment", DBDTstring, "File Comment",
			DBFFuser|DBFFarray|DBFFshow) != PDBFcomment)
		DMESG(DMCassert, eemsg);
	if (Field("Source", DBDTint, "Image Source Class",
			0, &PDBformatSource) != PDBFsource)
		DMESG(DMCassert, eemsg);
	if (Field("Make", DBDTstring, "Source Manufacturer",
			0) != PDBFmake)
		DMESG(DMCassert, eemsg);
	if (Field("Model", DBDTstring, "Source Model",
			0) != PDBFmodel)
		DMESG(DMCassert, eemsg);
	if (Field("Version", DBDTstring, "Source Software Version",
			0) != PDBFversion)
		DMESG(DMCassert, eemsg);
	if (Field("Quality", DBDTint, "Compression Quality (0-100)",
			0) != PDBFquality)
		DMESG(DMCassert, eemsg);
	if (Field("Crop", DBDTint, "Cropping Rectangle",
			DBFFuser|DBFFarray|DBFFhide) != PDBFcrop)
		DMESG(DMCassert, eemsg);
	if (Field("Orient", DBDTint, "Rotation for Display (degrees CW)",
			DBFFuser|DBFFhide) != PDBForient)
		DMESG(DMCassert, eemsg);
	if (Field("Flip", DBDTint, "Flip Horiz. (1) or Vert. (2)",
			DBFFuser|DBFFhide) != PDBFflip)
		DMESG(DMCassert, eemsg);
	if (Field("RedEye", DBDTint, "Red-eye Removal Regions",
			DBFFuser|DBFFarray|DBFFhide) != PDBFredeye)
		DMESG(DMCassert, eemsg);
	if (Field("SpotExp", DBDTint, "Spot-metered Exposure Regions",
			DBFFuser|DBFFarray|DBFFhide) != PDBFspotexp)
		DMESG(DMCassert, eemsg);
	if (Field("XSize", DBDTint, "Image Width",
			DBFFrequired) != PDBFxsize)
		DMESG(DMCassert, eemsg);
	if (Field("YSize", DBDTint, "Image Height",
			DBFFrequired) != PDBFysize)
		DMESG(DMCassert, eemsg);
	if (Field("XDensity", DBDTfloat, "Horizontal Density (ppi)",
			0) != PDBFxdensity)
		DMESG(DMCassert, eemsg);
	if (Field("YDensity", DBDTfloat, "Vertical Density (ppi)",
			0) != PDBFydensity)
		DMESG(DMCassert, eemsg);
	if (Field("XAngle", DBDTfloat, "Horizontal View Angle (degrees)",
			0) != PDBFxangle)
		DMESG(DMCassert, eemsg);
	if (Field("YAngle", DBDTfloat, "Vertical View Angle (degrees)",
			0) != PDBFyangle)
		DMESG(DMCassert, eemsg);
	if (Field("View", DBDTstring, "View Specification (Radiance)",
			0) != PDBFview)
		DMESG(DMCassert, eemsg);
	if (Field("StoNits", DBDTfloat, "Sample to Nits Conversion Factor",
			DBFFhide) != PDBFstonits)
		DMESG(DMCassert, eemsg);
	if (Field("ExpTime", DBDTfloat, "Exposure Duration (seconds)",
			0, &PDBformatExposure) != PDBFexptime)
		DMESG(DMCassert, eemsg);
	if (Field("Aperture", DBDTfloat, "Camera F-stop",
			0, &PDBformatAperture) != PDBFaperture)
		DMESG(DMCassert, eemsg);
	if (Field("ASA", DBDTfloat, "Sensitivity (ISO)",
			0) != PDBFasa)
		DMESG(DMCassert, eemsg);
	if (Field("Flash", DBDTint, "Flash Setting",
			0, &PDBformatFlash) != PDBFflash)
		DMESG(DMCassert, eemsg);
	if (Field("WhiteBal", DBDTint, "White Balance (illuminant)",
			0, &PDBformatWhiteBal) != PDBFwhitebal)
		DMESG(DMCassert, eemsg);
	if (Field("Focus", DBDTfloat, "Focus Distance (meters)",
			0) != PDBFfocus)
		DMESG(DMCassert, eemsg);
						// set default ordering
	sord.AddSort(PDBFfile);
	sord.AddSort(PDBFcapdate);
	sord.AddSort(PDBFsubject);
	sord.AddSort(PDBFowner);
}

// Initialize Pancine principal database search field information
PDBSearchFieldInfo::PDBSearchFieldInfo()
{
	static const char	eemsg[] = "define error in PDBSearchFieldInfo";

#define PS_ADDFIELD(id)	if (Field(PDBFInfo.GetName(id), \
					PDBFInfo.GetDetails(id)) < 0) \
				DMESG(DMCassert, eemsg)

	PS_ADDFIELD(PDBFsubject);
	PS_ADDFIELD(PDBFcapdate);
	PS_ADDFIELD(PDBFalbum);
	PS_ADDFIELD(PDBFkeyword);

#undef PS_ADDFIELD
}

// Set up query to match the given record
bool
PCsetQuery(DBQuery *dq, const DBRecord &dr, bool fileOnly)
{
	if (dq == NULL)
		return false;
	if (dr.GetFieldInfo() == NULL)
		return false;
	dq->rlo.Init(dr.GetFieldInfo());
	dq->rhi.Link(&dq->rlo);
						// basic file ID match
	if (!PDBsetField(&dq->rlo, PDBFfile, dr))
		return false;
	if (!PDBsetField(&dq->rlo, PDBFhashval, dr))
		return false;
	if (!PDBsetField(&dq->rlo, PDBFnbytes, dr))
		return false;
	if (fileOnly)				// just match file?
		return true;
	int	i;				// match user fields also
	for (i = dr.GetNAlloc(); i--; )
		if ( (dr.GetFieldDetails(i)->flags & (DBFFuser|DBFFarray))
					== DBFFuser )
			dq->rlo.SetField(i, dr[i]);
						// special cases
	if ((i = PDBfieldID(dr,PDBFkeyword)) >= 0)
		dq->rlo.SetField(i, dr[i]);
	if ((i = PDBfieldID(dr,PDBFcomment)) >= 0)
		dq->rlo.SetField(i, dr[i]);
	if ((i = PDBfieldID(dr,PDBFcrop)) >= 0)
		dq->rlo.SetField(i, dr[i]);
	if ((i = PDBfieldID(dr,PDBFspotexp)) >= 0)
		dq->rlo.SetField(i, dr[i]);
	if ((i = PDBfieldID(dr,PDBFredeye)) >= 0)
		dq->rlo.SetField(i, dr[i]);
	return true;
}

// Check to see if two records' orientations match
bool
PmatchOrient(const DBRecord &dr1, const DBRecord &dr2)
{
	if (&dr1 == &dr2)
		return true;

	int	flp1 = PDBgetField(dr1,PDBFflip).GetInt();
	int	flp2 = PDBgetField(dr2,PDBFflip).GetInt();
	int	rot1 = PDBgetField(dr1,PDBForient).GetInt();
	int	rot2 = PDBgetField(dr2,PDBForient).GetInt();
	if (flp1 == 3) {
		rot1 += 180;
		flp1 = 0;
	}
	if (flp2 == 3) {
		rot2 += 180;
		flp2 = 0;
	}
	if (flp1 != flp2)
		return false;
	while (rot1 < 0) rot1 += 360;
	while (rot1 >= 360) rot1 -= 360;
	while (rot2 < 0) rot2 += 360;
	while (rot2 >= 360) rot2 -= 360;
	return (rot1 == rot2);
}

// Update a record list according to a named change
bool
DBChangeList::ReflectChange(DBRecordList *rl, bool rev) const
{
	if (rl == NULL)
		return false;
	if ((rdel.GetSize() <= 0) & (radd.GetSize() <= 0))
		return false;
	const DBRecord *	rp;
	int			n;
	if (rev) {				// delete matching records
		rp = radd.GetArray();
		n = radd.GetSize();
	} else {
		rp = rdel.GetArray();
		n = rdel.GetSize();
	}
	for ( ; n-- > 0; rp++) {
		if (!rp->GetNAssigned())
			continue;
		for (int i = rl->GetSize(); i--; )
			if (rl->Get(i) == *rp) {
				(*rl)[i].Init();
				break;		// delete one at most
			}
	}
	rl->Compact();				// clear deadwood
	if (rev) {				// add new records
		rp = rdel.GetArray();
		n = rdel.GetSize();
	} else {
		rp = radd.GetArray();
		n = radd.GetSize();
	}
	int	last = rl->GetSize();
	if (rl->Resize(last+n) < last+n) {
		DMESG(DMCdata, "Resize failed in ReflectChange");
		return false;
	}
	for ( ; n-- > 0; rp++)
		if (rp->GetNAssigned())
			(*rl)[last++] = *rp;
	rl->Resize(last);
	if (prev == NULL || prev->name[0])
		return true;			// at the end of this change

	return prev->ReflectChange(rl, rev);	// else continue up list
}

// Apply pending change and create new one if name != NULL
bool
PDBAccess::Checkpoint(const char *name, PSnapshot *snap)
{
					// check special cases
	if (pending != NULL && (pending->radd.GetSize() <= 0) &
				(pending->rdel.GetSize() <= 0)) {
		if (pending->name[0])
			DMESGF(DMCwarning, "Empty change '%s' ignored",
					pending->name);
		delete pending;
		pending = NULL;
	}
	if (name != NULL) {		// delete redo list
		delete redoList;
		redoList = NULL;
	}
	if (pending == NULL) {
		if ((name != NULL) | (snap != NULL))
			pending = new DBChangeList(header, name, snap);
		return true;
	}
					// apply pending change
	int	nr;
	nr = DBAccess::DeleteRecords(pending->rdel.GetArray(), pending->rdel.GetSize());
	nunsynced += nr;
	if (nr < pending->rdel.GetSize())
		DMESGF(DMCwarning, "Unable to delete %d records in Checkpoint",
				pending->rdel.GetSize() - nr);
	nr = DBAccess::AddRecords(pending->radd.GetArray(), pending->radd.GetSize());
	nunsynced += nr;
	if (nr < pending->radd.GetSize()) {
		DMESGF(DMCdata, "Unable to add %d records in Checkpoint",
				pending->radd.GetSize() - nr);
		delete snap;
		return false;
	}
					// push change onto undo stack
	pending->prev = undoList;
	undoList = pending;
					// create new list
	if ((name != NULL) | (snap != NULL))
		pending = new DBChangeList(header, name, snap);
	else
		pending = NULL;

	return true;
}

// Undo last named change and push onto redo stack
const DBChangeList *
PDBAccess::Undo()
{
	if (pending != NULL) {			// move change to redo stack
		pending->prev = redoList;
		redoList = pending;
		pending = NULL;
		redoList->Restore(true);	// restore state
		if (redoList->name[0])
			return redoList;
	}
	DBChangeList *	cl = undoList;
	if (cl == NULL)
		return NULL;			// nothing to undo

	int	nr;				// apply change in reverse
	nr = DBAccess::DeleteRecords(cl->radd.GetArray(), cl->radd.GetSize());
	nunsynced += nr;
	if (nr != cl->radd.GetSize())
		DMESGF(DMCwarning, "Unable to unadd %d records in Undo",
				cl->radd.GetSize() - nr);
	nr = DBAccess::AddRecords(cl->rdel.GetArray(), cl->rdel.GetSize());
	nunsynced += nr;
	if (nr != cl->rdel.GetSize()) {
		DMESGF(DMCdata, "Unable to undelete %d records in Undo",
				cl->rdel.GetSize() - nr);
		return NULL;
	}	
	undoList = cl->prev;			// move from undo to redo
	cl->prev = redoList;
	redoList = cl;
	cl->Restore(true);			// restore state

	if (cl->name[0])			// was it named?
		return cl;

	return Undo();				// else try next one
}

// Redo last named change we undid
// User can undo<->redo all day and all we do is swap pointers
const DBChangeList *
PDBAccess::Redo()
{
	DBChangeList *	cl = redoList;
	if (cl == NULL)
		return NULL;			// nothing to redo

	if (!Checkpoint())			// apply pending changes
		return NULL;

	redoList = cl->prev;			// pop top off redo list
	cl->prev = NULL;
	pending = cl;				// make it our next change
	cl->Restore(false);			// restore state

	if (cl->name[0])			// was it named?
		return cl;

	return Redo();				// else try next one
}

// Pretend to add records to database by adding them to our change list
int
PDBAccess::AddRecords(const DBRecord *rlist, int n)
{
	if ((n <= 0) | (rlist == NULL))
		return 0;
	if (pending == NULL)
		pending = new DBChangeList(header);
	int	i;
	i = pending->radd.GetSize();
	if (!pending->radd.Resize(i + n)) {
		DMESG(DMCparameter, "Resize failed in AddRecords");
		return 0;
	}
	DBRecord *	rp = pending->radd.Array() + i;
	for (i = 0; i < n; i++)		// copy each record
		*rp++ = rlist[i];
	return n;
}
	
// Pretend to delete matching records from database by putting in change list
int
PDBAccess::DeleteRecords(const DBRecord *rlist, int n)
{
	if ((n <= 0) | (rlist == NULL))
		return 0;
	if (pending == NULL)
		pending = new DBChangeList(header);
	ABitMap	rMap(n);
	int	i;
	if (n * pending->radd.GetSize() <= IO_CHECKLIM) {
		int	j;			// short-circuit add->delete
		for (i = 0; i < n; i++)
			if ((j = pending->radd.FindRecord(rlist[i])) >= 0) {
				pending->radd[j].Init();
				rMap.Set(i);
			}
		pending->radd.Compact();
	}
	i = pending->rdel.GetSize();
	if (!pending->rdel.Resize(i + rMap.SumTotal(false))) {
		DMESG(DMCparameter, "Resize failed in DeleteRecords");
		return 0;
	}
	DBRecord *	rp = pending->rdel.Array() + i;
	for (i = 0; i < n; i++)		// copy each record
		if (!rMap.Check(i))
			*rp++ = rlist[i];
	return n;
}

static const int	Pjdate[12] = {0, 31, 59, 90, 120, 151,
					181, 212, 243, 273, 304, 334};

// Convert "YYYY:MM:DD HH:MM:SS" to seconds from midnight Jan 1, P_FIRSTYR
unsigned long
PsecsFromDate(const PDate dt)
{
	int		yr, mo, da, hr, mn, sc;
	if (dt == NULL)
		return 0L;
	if (strlen(dt) != 19)
		return 0L;
	if (sscanf(dt, "%d:%d:%d %d:%d:%d", &yr, &mo, &da, &hr, &mn, &sc) != 6)
		return 0L;
	if ((yr < P_FIRSTYR) | (mo < 1) | (mo > 12) | (da < 1) | (da > 31) |
			(hr < 0) | (hr > 23) | (mn < 0) | (mn > 59) | (sc < 0) | (sc > 59))
		return 0L;
	unsigned long	sum;
	sum = sc + 60L*(mn + 60L*(hr + 24L*(da-1 + Pjdate[mo-1])));
	sum += 365L*24L*3600L*(yr - P_FIRSTYR);
	sum += 24L*3600L*(1 + ((yr - P_FIRSTYR)>>2));
	if (!(yr & 3) & (mo <= 2))
		sum -= 24L*3600L;
	return sum;
}

#if 0
// Reverse conversion gets back same date as went into PSecsFromDate()
char *
PdateFromSecs(PDate dt, unsigned long sec)
{
	if ((sec <= 0) | (dt == NULL))
		return NULL;
	
	yr = P_FIRSTYR + sec/(365L*24L*3600L);
	
	int	yr, mo, da, hr, mn, sc;
	sprintf(dt, "%04d:%02d:%02d %02d:%02d:%02d", yr, mo, da, hr, mn, sc);
	return dt;
}
#endif

// Format folder name, emphasizing last part of path
char *
PDBformatFolder(char *buf, int len, const DBField &f)
{
	if ((buf == NULL) | (len < 4))
		return NULL;
	const char *	sval;
	if (!f.Get(&sval)) {
		buf[0] = '\0';
		return NULL;
	}
	int	dlen = strlen(sval);
	if (len >= dlen)
		return strcpy(buf, sval);
	strcpy(buf, "...");
	strcpy(buf+3, sval+(dlen-len+3));
	return buf;
}

// Format date and time into buffer of specified length
char *
PDBformatDate(char *buf, int len, const DBField &f)
{
	static const char	Month[12][6] = {"Jan", "Feb", "March", "April",
						"May", "June", "July", "Aug",
						"Sept", "Oct", "Nov", "Dec"};
	static const char *	format[] = {
		"Month D YYYY h:mm:ss@",
		"Month D YYYY h:mm@",
		"M/DD/YYYY h:mm@",
		"M/DD/YY HH:mm",
		"M/DD/YYYY",
		"M/DD/YY",
		"M/YY",
		NULL
	};
	if ((buf == NULL) | (len < 2))
		return NULL;
	buf[0] = '\0';
	const char *	din;
	if (!f.Get(&din))
		return NULL;
	int	yr, mo, da, hr, mn, sc;
	if (strlen(din) != 19 || sscanf(din, "%d:%d:%d %d:%d:%d",
					&yr, &mo, &da, &hr, &mn, &sc) != 6 ||
			(yr < P_FIRSTYR) | (mo < 1) | (mo > 12) | (da < 1) | (da > 31) |
			(hr < 0) | (hr > 23) | (mn < 0) | (mn > 59) | (sc < 0) | (sc > 59))
		return NULL;
	char *		dout = buf;
	const char *	fmt;
	int		i;
	for (i = 0; format[i] != NULL && (int)strlen(format[i]) > len; i++)
		;
	if ((fmt = format[i]) == NULL)
		return buf;			// nothing fits!
	while (*fmt) {
		int	cnt = 1;
		switch (*fmt) {
		case 'Y':					// year
			while (fmt[1] == 'Y') { fmt++; cnt++; }
			if (cnt == 2)
				sprintf(dout, "%02d", yr % 100);
			else if (cnt == 4)
				sprintf(dout, "%04d", yr);
			break;
		case 'M':					// month
			if (!strncmp(fmt, "Month", 5)) {
				fmt += 4;
				strcpy(dout, Month[mo-1]);
				break;
			}
			while (fmt[1] == 'M') { fmt++; cnt++; }
			if (cnt == 1)
				sprintf(dout, "%d", mo);
			else if (cnt == 2)
				sprintf(dout, "%02d", mo);
			break;
		case 'D':					// day
			while (fmt[1] == 'D') { fmt++; cnt++; }
			if (cnt == 1)
				sprintf(dout, "%d", da);
			else if (cnt == 2)
				sprintf(dout, "%02d", da);
			break;
		case 'H':					// 24-hour
			while (fmt[1] == 'H') { fmt++; cnt++; }
			if (cnt == 1)
				sprintf(dout, "%d", hr);
			else if (cnt == 2)
				sprintf(dout, "%02d", hr);
			break;
		case 'h':					// 12-hour
			while (fmt[1] == 'h') { fmt++; cnt++; }
			if (cnt == 1)
				sprintf(dout, "%d", (hr==0) ? 12 :
						(hr>12) ? hr-12 : hr);
			else if (cnt == 2)
				sprintf(dout, "%02d", (hr==0) ? 12 :
						(hr>12) ? hr-12 : hr);
			break;
		case 'm':					// minute
			while (fmt[1] == 'm') { fmt++; cnt++; }
			if (cnt == 1)
				sprintf(dout, "%d", mn);
			else if (cnt == 2)
				sprintf(dout, "%02d", mn);
			break;
		case 's':					// second
			while (fmt[1] == 's') { fmt++; cnt++; }
			if (cnt == 1)
				sprintf(dout, "%d", sc);
			else if (cnt == 2)
				sprintf(dout, "%02d", sc);
			break;
		case '@':					// am/pm
			while (fmt[1] == '@') { fmt++; cnt++; }
			if (cnt == 1)
				strcpy(dout, hr>=12 ? "p" : "a");
			else if (cnt == 2)
				strcpy(dout, hr>=12 ? "pm" : "am");
			break;
		default:					// punctuation?
			dout[0] = *fmt; dout[1] = '\0';
			break;
		}
		while (*dout) dout++;
		fmt++;
	}
	return buf;
}

// Format exposure time
char *
PDBformatExposure(char *buf, int len, const DBField &f)
{
	if ((buf == NULL) | (len < 7))
		return NULL;
	float	fval;
	if (!f.Get(&fval)) {
		buf[0] = '\0';
		return NULL;
	}
	if (fval > .5f) {
		int	maxSD = len;
		if (maxSD > 4)
			maxSD = 4;
		if (DBformatFloat(buf, maxSD, f) == NULL)
			return NULL;
		if (len - strlen(buf) > 4)
			strcat(buf, " sec");
		return buf;
	}
	fval = 1.f/fval;				// round-off
	if (fabs(fval - int(fval+.5f)) < .01f*fval)
		fval = int(fval+.5f);
	DBField	recipr = fval;
	buf[0] = '1'; buf[1] = '/';
	if (DBformatFloat(buf+2, len-2, recipr) == NULL) {
		buf[0] = '\0';
		return NULL;
	}
	if (len - strlen(buf) > 4)
		strcat(buf, " sec");
	return buf;
}

// Format lens aperture (f-stop)
char *
PDBformatAperture(char *buf, int len, const DBField &f)
{
	if (buf == NULL || len < 5)
		return NULL;
	if (len > 8)
		len = 8;				// forces round-off
	float	fval;
	if (!f.Get(&fval)) {
		buf[0] = '\0';
		return NULL;
	}
	buf[0] = 'f'; buf[1] = '/';
	if (DBformatFloat(buf+2, len-2, fval) == NULL) {
		buf[0] = '\0';
		return NULL;
	}
	return buf;
}

// Translate flash setting
char *
PDBformatFlash(char *buf, int len, const DBField &f)
{
	if ((buf == NULL) | (len < 4))
		return NULL;
	int32	flash;
	if (!f.Get(&flash)) {
		buf[0] = '\0';
		return NULL;
	}
	const char *	fmode = "UNKNOWN";
	switch (flash & IFMmask) {
	case IFMoff:
		fmode = "off";
		break;
	case IFMon:
	case IFMstrobe:
		fmode = "on";
		break;
	case IFMnostrobe:
		fmode = "unreturned";
		break;
	}
	strlcpy(buf, fmode, len);
	return buf;
}

// Translate white balance setting
char *
PDBformatWhiteBal(char *buf, int len, const DBField &f)
{
	if ((buf == NULL) | (len < 2))
		return NULL;
	int32	illum;
	if (!f.Get(&illum)) {
		buf[0] = '\0';
		return NULL;
	}
	const char *	illumtype = "UNKNOWN";
	switch (illum) {
	case IWBauto:
		illumtype = "auto";
		break;
	case IWBdaylight:
	case IWBdaylight1:
		illumtype = "daylight";
		break;
	case IWBfluor:
	case IWBfluor1:
		illumtype = "fluorescent";
		break;
	case IWBincand:
	case IWBincand1:
		illumtype = "incandescent";
		break;
	case IWBovercast:
		illumtype = "overcast";
		break;
	case IWBflash:
		illumtype = "flash";
		break;
	case IWBshadow:
		illumtype = "shadow";
		break;
	case IWBstilla:
		illumtype = "Illum A";
		break;
	case IWBstillb:
		illumtype = "Illum B";
		break;
	case IWBstillc:
		illumtype = "Illum C";
		break;
	case IWBd50:
		illumtype = "5000K daylight";
		break;
	case IWBd55:
		illumtype = "5500K daylight";
		break;
	case IWBd65:
		illumtype = "6500K daylight";
		break;
	case IWBd75:
		illumtype = "7500K daylight";
		break;
	case IWBother:
		illumtype = "other";
		break;
	}
	strlcpy(buf, illumtype, len);
	return buf;
}

// Translate source ID
char *
PDBformatSource(char *buf, int len, const DBField &f)
{
	if ((buf == NULL) | (len < 2))
		return NULL;
	int32	sid;
	if (!f.Get(&sid)) {
		buf[0] = '\0';
		return NULL;
	}
	bool		shortName = (len < 18);
	const char *	stype = "UNKNOWN";
	switch (sid) {
	case ISTdigicam:
		stype = shortName ? "digicam" : "digital camera";
		break;
	case ISTfilm:
		stype = "film";
		break;
	case ISTflatbed:
		stype = shortName ? "scanner" : "flatbed scanner";
		break;
	case ISTrangefinder:
		stype = shortName ? "laser" : "laser range finder";
		break;
	case ISTeditor:
		stype = shortName ? "editor" : "image editor";
		break;
	case ISTrender:
		stype = shortName ? "renderer" : "rendering software";
		break;
	}
	strlcpy(buf, stype, len);
	return buf;
}
