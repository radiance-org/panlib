/*
 *  pimage.cpp
 *  panlib
 *
 *  Pancine image processing routines.
 *
 *  Created by gward on Wed Jun 13 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#include "pancine.h"
#include <ctype.h>
#include <math.h>

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

#define POI_MAXD	16			// max. directories to check

PFileMigration *	PFMigrations = NULL;	// our file migration list

// Compute change in directories between two paths
bool
PFileMigration::GetChange(char *d1, char *d2, const char *p1, const char *p2)
{
	if ((d1 == NULL) | (d2 == NULL) | (p1 == NULL) | (p2 == NULL))
		return false;
	char	*ocp, *ncp;
						// copy paths
	for (ocp = d1; (*ocp = *p1++); ++ocp)
		;
	for (ncp = d2; (*ncp = *p2++); ++ncp)
		;
						// eliminate common tails
	while ((ocp > d1) & (ncp > d2) && *--ocp == *--ncp)
		;
	if ((ocp <= d1) & (ncp <= d2))		// check exact match
		return false;
						// check special cases
	if ((ocp <= d1) | (ncp <= d2) || IS_DIRSEP(*ocp) ^ IS_DIRSEP(*ncp)) {
		++ocp; ++ncp;
	}
	while (*ocp && !IS_DIRSEP(*ocp))
		++ocp;
	*ocp = '\0';
	while (*ncp && !IS_DIRSEP(*ncp))
		++ncp;
	*ncp = '\0';
	return true;
}

// Initialize new migration entry -- op and np must be different!
PFileMigration::PFileMigration(const char *op, const char *np,
				PFileMigration *nfm)
{
	char	olddir[1024], newdir[1024];
						// save different prefixes
	next = nfm;
	if (!GetChange(olddir, newdir, op, np))
		DMESG(DMCassert, "False file migration!");
	dold = strlink(olddir);
	dnew = strlink(newdir);
	doldlen = strlen(dold);
	dnewlen = strlen(dnew);
	sprintf(dmessage_buf, "Added migration path: '%s' -> '%s'", dold, dnew);
	DMESG(DMCinfo, dmessage_buf);

}

// Add new migration path to end of list
PFileMigration *
PFileMigration::AddMigration(const char *op, const char *np)
{
	char	olddir[1024], newdir[1024];
						// see if it's in our list
	if (!GetChange(olddir, newdir, op, np))
		DMESG(DMCassert, "False file migration!");
	PFileMigration *	fmp = this;
	for ( ; ; ) {
		if (!strcmp(fmp->dold, olddir) && !strcmp(fmp->dnew, newdir))
			return fmp;
		if (fmp->next == NULL)
			break;
		fmp = fmp->next;
	}
	return (fmp->next = new PFileMigration(olddir, newdir));
}

// Promote the given member to this position
PFileMigration *
PFileMigration::Promote(const PFileMigration *fmp)
{
	if (fmp == this)
		return this;
	PFileMigration *	flast = this;
	PFileMigration *	fsrch;
	while ((fsrch = flast->next) != NULL) {
		if (fsrch == fmp) {
			flast->next = fsrch->next;
			fsrch->next = this;
			return fsrch;
		}
		flast = fsrch;
	}
	return NULL;		// not in list!
}

// Find next migration entry that matches this path prefix
const PFileMigration *
PFileMigration::Match(const char *op) const
{
	const PFileMigration *	fmp;
	for (fmp = this; fmp != NULL; fmp = fmp->next) {
		const char *	mcp = fmp->dold;
		const char *	tcp = op;
		if (mcp == tcp)
			break;
		while (*mcp && *mcp == *tcp++)
			++mcp;
		if (!*mcp && (IS_DIRSEP(*tcp) || !*tcp))
			break;
	}
	return fmp;
}

// Find next migration entry that matches this path pair
const PFileMigration *
PFileMigration::Match(const char *op, const char *np) const
{
	char	olddir[1024], newdir[1024];
	if (!GetChange(olddir, newdir, op, np))
		return NULL;
	const PFileMigration *	fmp;
	for (fmp = this; fmp != NULL; fmp = fmp->next)
		if (!strcmp(fmp->dold, olddir) && !strcmp(fmp->dnew, newdir))
			break;
	return fmp;
}

// Generate a new path name using matching prefix substitution
bool
PFileMigration::Migrate(char *np, int len, const char *op) const
{
	if ((np == NULL) | (op == NULL) | (np == op))
		return false;
	if (strncmp(op, dold, doldlen))
		return false;
	if (!IS_DIRSEP(op[doldlen]) && op[doldlen] != '\0')
		return false;
	if (len < (int)strlen(op) - doldlen + dnewlen)
		return false;
	strcpy(np, dnew);
	strcpy(np+dnewlen, op+doldlen);
	return true;
}

// Record new file name and location (& modified time)
static void
recordFilePath(DBRecord *dr, const char *fpath)
{
	time_t		modtime = getfiletime(fpath);
	const char *	fnpos = PgetFilename(fpath);
						// set file modified time
	if (modtime != 0)
		PDBsetField(dr, PDBFmodtime, (int32)modtime);
						// set file name
	PDBsetField(dr, PDBFfile, fnpos);
						// set directory location
	int	dirlen = fnpos - fpath;
	if (dirlen <= 0)
		return;
	const char *	dnl[POI_MAXD];
	int		ndn = PDBgetField(*dr,PDBFlocation).Get(dnl,POI_MAXD);
	char		dirname[1024];
	int		i;
	if (dirlen > 1)				// remove final DIRSEP
		dirlen--;
	if (dirlen >= (int)sizeof(dirname))
		dirlen = sizeof(dirname)-1;
	strncpy(dirname, fpath, dirlen);
	dirname[dirlen] = '\0';
						// insert at head of list
	for (i = 0; i < ndn; i++)
		if (!strcmp(dnl[i], dirname))
			break;
	if (i >= ndn)
		if (ndn < POI_MAXD)
			++ndn;
		else
			--i;
	while (i--)
		dnl[i+1] = dnl[i];
	dnl[0] = dirname;
	PDBarrayField(dr, PDBFlocation, dnl, ndn);
}

// Locate Pancine image file and confirm size and hash value
PLocateFileRes
PlocateFile(char *pnbuf, int blen, const DBRecord &dr)
{
	if ((pnbuf == NULL) | (blen <= 0))
		return PLFRfail;
	if (!dr.GetNAlloc())
		return PLFRfail;
	const char *	fn;
	const char *	dnl[POI_MAXD];
	int		ndn = PDBgetField(dr,PDBFlocation).Get(dnl,POI_MAXD);
	if (ndn <= 0 || !PDBgetField(dr,PDBFfile).Get(&fn))
		return PLFRfail;
	char		mypbuf[1024];
	int32		dhval = 0, dnbyt = 0;
	int32		fhval = 0, fnbyt = 0;
	int32		*fhvalp = &fhval, *fnbytp = &fnbyt;
	int32		dmtim;
	int		di;
						// get files & directories
	pnbuf[0] = '\0';
	if (!PDBgetField(dr,PDBFnbytes).Get(&dnbyt))
		fnbytp = NULL;
	if (!PDBgetField(dr,PDBFhashval).Get(&dhval))
		fhvalp = NULL;
	if (!PDBgetField(dr,PDBFmodtime).Get(&dmtim))
		dmtim = 0;
						// check all combinations
	PLocateFileRes	rval = PLFRfail;
	for (di = 0; di < ndn; di++) {
		snprintf(mypbuf, blen, "%s%c%s", dnl[di], DIRSEP, fn);
		if ((int)strlen(mypbuf) >= blen)
			continue;		// should be error?
		int32	nmodtime = (int32)getfiletime(mypbuf);
		if ((dmtim != 0) & (dmtim == nmodtime)) {
			fhval = dhval; fnbyt = dnbyt;
		} else if (!PhashFile(fhvalp, fnbytp, mypbuf))
			continue;
		if ((fhval == dhval) & (fnbyt == dnbyt)) {
			strcpy(pnbuf, mypbuf);
			rval = (!di & (dmtim==nmodtime)) ? PLFRok : PLFRputback;
			break;			// found exact match
		}
		if (rval == PLFRfail) {		// found file at least
			strcpy(pnbuf, mypbuf);
			rval = PLFRchanged;
		}
	}
	if (rval != PLFRfail)			// got something?
		return rval;
						// check for file migration
	for (const PFileMigration *fmp = PFMigrations->Match(dnl[0]);
			fmp != NULL; fmp = fmp->GetNext()->Match(dnl[0])) {
		if (!fmp->Migrate(mypbuf, sizeof(mypbuf), dnl[0]))
			continue;
		char *	cp = mypbuf;
		while (*cp) cp++;
		*cp++ = DIRSEP;
		strcpy(cp, fn);
		if ((int)strlen(mypbuf) >= blen)
			continue;
		if (!PhashFile(fhvalp, fnbytp, mypbuf))
			continue;
		if ((fhval == dhval) & (fnbyt == dnbyt)) {
			strcpy(pnbuf, mypbuf);
			rval = PLFRmigrated;
			PFMigrations = PFMigrations->Promote(fmp);
			DASSERT(PFMigrations != NULL);
			break;			// found exact match
		}
		if (rval == PLFRfail) {		// found file at least
			strcpy(pnbuf, mypbuf);
			rval = PLFRchanged;
		}
	}
	return rval;
}

// Locate and record image file location, enlisting user's help if necessary
PLocateFileRes
PlocateFile(char *pnbuf, int blen, DBRecord *dr,
			bool (*queryf)(char *, int, DBRecord *))
{
	if ((pnbuf == NULL) | (blen <= 0) | (dr == NULL))
		return PLFRfail;
	if (!dr->GetNAlloc())
		return PLFRfail;
	PLocateFileRes	rval = PlocateFile(pnbuf, blen, *dr);
	switch (rval) {
	case PLFRok:
		return rval;
	case PLFRputback:
	case PLFRrenamed:
	case PLFRmigrated:
		recordFilePath(dr, pnbuf);
		return rval;
	case PLFRchanged:
		PgrokFile(dr, pnbuf);
		return rval;
	case PLFRfail:
		break;
	}
	if (queryf == NULL)			// need user assist
		return rval;
	char	prevpath[1024];			// remember previous path
	char *	fname = PcurrentFile(prevpath, sizeof(prevpath), *dr);
	if (!(*queryf)(pnbuf, blen, dr))	// call inquiry function
		return PLFRfail;
	int32	hval, flen;
	if (!PhashFile(&hval, &flen, pnbuf))	// check file and contents
		return PLFRfail;
	int32	phval, pflen;
	if (PDBgetField(*dr,PDBFhashval).Get(&phval) &&
			PDBgetField(*dr,PDBFnbytes).Get(&pflen) &&
			(hval != phval) | (flen != pflen)) {
		if (strcmp(fname, PgetFilename(pnbuf)))
			return PLFRfail;	// new name AND new contents?!
		PgrokFile(dr, pnbuf);
		return PLFRchanged;
	}
	PDBsetField(dr, PDBFhashval, hval);
	PDBsetField(dr, PDBFnbytes, flen);
	if (!strcmp(pnbuf, prevpath))
		return PLFRok;			// original reappeared
						// new path, same contents
	recordFilePath(dr, pnbuf);
						// check for renamed file
	if (strcmp(fname, PgetFilename(pnbuf)))
		return PLFRrenamed;
						// record migration if new
	if (PFMigrations->Match(prevpath, pnbuf) == NULL)
		PFMigrations = new PFileMigration(prevpath, pnbuf, PFMigrations);
	return PLFRmigrated;
}

// Get current file path without any checks at all (return filename pointer)
char *
PcurrentFile(char *pnbuf, int blen, const DBRecord &dr)
{
	if ((pnbuf == NULL) | (blen <= 0))
		return NULL;
	if (!dr.GetNAlloc())
		return NULL;
	char *	fname = pnbuf;
	if (PDBgetField(dr,PDBFlocation).Get(pnbuf, blen-1)) {
		while (*fname) fname++;
		*fname++ = DIRSEP;
	}
	*fname = '\0';
	if (PDBgetField(dr,PDBFfile).Get(fname, pnbuf+blen-fname))
		return fname;
	return NULL;
}

// Get hash value and length for file (do NOT alter this code)
bool
PhashFile(int32 *hval, int32 *flen, const char *fname)
{
	if (fname == NULL)
		return false;
	int	fd = open(fname, O_RDONLY);
	if (fd < 0)
		return false;
	if ((hval == NULL) & (flen == NULL)) {	// checking readability
		close(fd);
		return true;
	}
	SET_FD_BINARY(fd);
						// get file length
	off_t	fileLen = lseek(fd, (off_t)0, SEEK_END);
	if (flen != NULL)
		*flen = (int32)fileLen;
						// hash file
	if (hval != NULL) {
		int		maxReads = 5;	// limit read calls
		off_t		skipAmount = 0;
		char		hashbuf[P_NHASHBITS*512];
		if (fileLen > (off_t)sizeof(hashbuf)*maxReads)
			skipAmount = 1 + ( fileLen -
					(off_t)sizeof(hashbuf)*maxReads ) /
						(maxReads - 1);
		lseek(fd, (off_t)0, SEEK_SET);
		*hval = 0;
		while (maxReads--) {
			int	nr = read(fd, hashbuf, sizeof(hashbuf));
			if (nr <= 0)
				break;
			*hval ^= memhash(hashbuf, nr, P_NHASHBITS);
			if (nr < (int)sizeof(hashbuf))
				break;
			if ((maxReads > 0) & (skipAmount > 0))
				lseek(fd, skipAmount, SEEK_CUR);
		}
	}
	close(fd);
	return true;
}

// Assign record fields based on information from file
static void
recordFileInfo(DBRecord *dr, const char *fpath)
{
						// set file path and modified time
	recordFilePath(dr, fpath);
						// set file length & hash
	int32	hval, flen;
	if (PhashFile(&hval, &flen, fpath)) {
		PDBsetField(dr, PDBFhashval, hval);
		PDBsetField(dr, PDBFnbytes, flen);
	}
	char	datime[20];
	if (!PDBgetField(*dr,PDBFcapdate).GetNV() && getfiledate(datime,fpath) != NULL)
		PDBsetField(dr, PDBFcapdate, datime);
}

// Assign record fields based on info from open image
static void
recordImageInfo(DBRecord *dr, ImgReader *ir)
{
						// set size and format
	PDBsetField(dr, PDBFxsize, (int32)ir->xres);
	PDBsetField(dr, PDBFysize, (int32)ir->yres);
	char		fmtdesc[128];
	if (ir->ri->suffixList != NULL)
		strlcpy(fmtdesc, ir->ri->suffixList, sizeof(fmtdesc));
	else
		fmtdesc[0] = '\0';
	char *	cp = fmtdesc;
	while ((*cp != '\0') & (*cp != '.'))
		cp++;
	*cp++ = ' ';
	if (ir->encoding != NULL)		// explicit encoding?
		strcpy(cp, ir->encoding);
	else					// else describe color space
		PdescribeCS(&ir->cs, cp);

	PDBsetField(dr, PDBFformat, fmtdesc);
						// set info parameters
	ImgInfo	info;
	if (IRgetInfo(ir, &info) != IREnone)
		return;
	if (info.flags & IIForientation)	// image orientation
		switch (info.orientation) {
		case IOtopleft:
			/* Ignore default setting...
			PDBsetField(dr, PDBForient, (int32)0);
			PDBsetField(dr, PDBFflip, (int32)0);
			*/
			break;
		case IOtopright:
			PDBsetField(dr, PDBForient, (int32)0);
			PDBsetField(dr, PDBFflip, (int32)1);
			break;
		case IObotright:
			PDBsetField(dr, PDBForient, (int32)180);
			PDBsetField(dr, PDBFflip, (int32)0);
			break;
		case IObotleft:
			PDBsetField(dr, PDBForient, (int32)180);
			PDBsetField(dr, PDBFflip, (int32)1);
			break;
		case IOlefttop:
			PDBsetField(dr, PDBForient, (int32)270);
			PDBsetField(dr, PDBFflip, (int32)2);
			break;
		case IOrighttop:
			PDBsetField(dr, PDBForient, (int32)90);
			PDBsetField(dr, PDBFflip, (int32)0);
			break;
		case IOrightbot:
			PDBsetField(dr, PDBForient, (int32)90);
			PDBsetField(dr, PDBFflip, (int32)2);
			break;
		case IOleftbot:
			PDBsetField(dr, PDBForient, (int32)270);
			PDBsetField(dr, PDBFflip, (int32)0);
			break;
		}
	if (info.flags & IIFquality)		// compression quality
		PDBsetField(dr, PDBFquality, (int32)info.quality);
						// crop rectangle
	if (info.flags & IIFcrop && PlegalRect(&info.crop,ir->xres,ir->yres)) {
		int32	crect[4];
		crect[0] = info.crop.xleft;
		crect[1] = info.crop.ytop;
		crect[2] = info.crop.xright;
		crect[3] = info.crop.ybottom;
		PDBarrayField(dr, PDBFcrop, crect, 4);
	}
	if (info.flags & IIFcapdate)		// capture date
		PDBsetField(dr, PDBFcapdate, info.capdate);
	if (info.flags & IIFgmtdate)		// precise GMT capture time
		PDBsetField(dr, PDBFgmtdate, info.gmtdate);
	if (info.flags & IIFhdensity) {		// density
		PDBsetField(dr, PDBFxdensity, float(info.hdensity));
		PDBsetField(dr, PDBFydensity,
					float(info.hdensity*ir->pixAspect));
	} else if ((ir->pixAspect < .98f) | (ir->pixAspect > 1.02f)) {
		PDBsetField(dr, PDBFxdensity, 72.f);
		PDBsetField(dr, PDBFydensity, float(72.f*ir->pixAspect));
	}
	if (info.flags & IIFhvangle) {		// image angles
		PDBsetField(dr, PDBFxangle, float(info.hvangle));
		if (info.hvangle < 125.)
			PDBsetField(dr, PDBFyangle,
				float(360./M_PI*atan(
					ir->yres/(ir->pixAspect*ir->xres) *
					tan(M_PI/360.*info.hvangle) )));
		else
			PDBsetField(dr, PDBFyangle,
					float(ir->yres/(ir->pixAspect*ir->xres) *
						info.hvangle));	// XXX probably wrong
	}
	if (info.flags & IIFsource) {		// image source
		PDBsetField(dr, PDBFsource, (int32)info.source.stype);
		PDBsetField(dr, PDBFmake, info.source.make);
		PDBsetField(dr, PDBFmodel, info.source.model);
		PDBsetField(dr, PDBFversion, info.source.vers);
	}
	if (info.flags & IIFstonits)		// sample to Nits conversion
		PDBsetField(dr, PDBFstonits, info.stonits);
	if (info.flags & IIFfocus)		// focus distance
		PDBsetField(dr, PDBFfocus, info.focus);
	if (info.flags & IIFlatlong) {		// latitude/longitude
		PDBsetField(dr, PDBFlatitude, info.latlong[0]);
		PDBsetField(dr, PDBFlongitude, info.latlong[1]);
	}
	if (info.flags & IIFaltitude)		// altitude
		PDBsetField(dr, PDBFaltitude, info.altitude);
	if (info.flags & IIFbearing)		// sight bearing
		PDBsetField(dr, PDBFbearing, info.bearing);
	if (info.flags & IIFflash)		// flash setting
		PDBsetField(dr, PDBFflash, (int32)info.flash);
	if (info.flags & IIFexptime)		// exposure seconds
		PDBsetField(dr, PDBFexptime, info.exptime);
	if (info.flags & IIFaperture)		// f-stop
		PDBsetField(dr, PDBFaperture, info.aperture);
	if (info.flags & IIFasa)		// ISO ASA
		PDBsetField(dr, PDBFasa, info.asa);
	if (info.flags & IIFwhitebal)		// white balance setting
		PDBsetField(dr, PDBFwhitebal, (int32)info.whitebal);
	if (info.flags & IIFview)		// Radiance view parameters
		PDBsetField(dr, PDBFview, info.view);
	if (info.flags & IIFowner)		// Content owner
		PDBsetField(dr, PDBFowner, info.owner);
						// Comments & Params -> history
	DBField	tsum = PDBgetField(*dr, PDBFhistory);
	DBField	tnew;
	if (info.flags & IIFcomments && tnew.SetText(info.comments))
		tsum |= tnew;
	if (info.flags & IIFparams && tnew.SetText(info.params))
		tsum |= tnew;
	if (info.flags & (IIFcomments|IIFparams))
		PDBsetField(dr, PDBFhistory, tsum);
}

// Set record fields based on information from file
void
PgrokFile(DBRecord *dr, const char *fpath)
{
	if (dr == NULL)
		return;
	if (dr->GetFieldInfo() == NULL)
		dr->Init(&PDBFInfo);
						// set file info
	recordFileInfo(dr, fpath);
						// set image info
	ImgReader *	ir = PopenImageF(fpath, true, NULL);
	if (ir != NULL) {
		if (ir->errCode == IREnone)	// always true
			recordImageInfo(dr, ir);
		IRclose(ir);
	}
}

// Set record fields based on info from open image
void
PgrokImage(DBRecord *dr, ImgReader *ir)
{
	if ((dr == NULL) | (ir == NULL))
		return;
	if (dr->GetFieldInfo() == NULL)
		dr->Init(&PDBFInfo);
						// set file info
	recordFileInfo(dr, ir->file);
						// set image info
	recordImageInfo(dr, ir);
}

// Open an image from the given Pancine image record
ImgReader *
PopenImageD(const DBRecord &pr, bool quiet)
{
	if (!pr.GetNAlloc())
		return NULL;
	const char *	fn;			// get file name
	if (!PDBgetField(pr,PDBFfile).Get(&fn)) {
		if (!quiet)
			DMESG(DMCdata, "Missing image file name");
		return NULL;
	}
						// get reader list
	const ImgReaderInterface *	irilist[P_MAXREADER];
	int				nri;
	char				suffix[64];
	int				i;
	if (PDBgetField(pr,PDBFformat).Get(suffix, sizeof(suffix))) {
		char *	sp = suffix;		// designated format reader(s)
		while (*sp && !isspace(*sp))
			sp++;
		*sp = '\0';
	} else {				// else use file name
		const char *	sp = PgetSuffix(fn);
		if (sp != NULL)
			strcpy(suffix, sp);
		else
			suffix[0] = '\0';
	}
						// get reader list
	nri = PtargetReader(irilist, suffix);
	if (nri <= 0) {
		if (!quiet)
			DMESGF(DMCdata, "%s: unknown image format", fn);
		return NULL;
	}
	char		fname[1024];		// get pathname
	PLocateFileRes	lfr = PlocateFile(fname, sizeof(fname), pr);
	if (lfr == PLFRfail) {
		if (!quiet)
			DMESGF(DMCresource, "%s: cannot find image", fn);
		return NULL;
	}
	for (i = 0; i < nri; i++) {		// try each interface
		ImgReader *	ir = PopenImageF(fname, quiet, irilist[i]);
		if (ir != NULL)
			return ir;		// success!
	}
	if (!quiet)
		DMESGF(DMCdata, "%s: image format error", fname);
	return NULL;				// return failure
}

// Set image information based on database record and return pixel aspect
float
PsetInfo(ImgInfo *iip, const DBRecord &pr)
{
	if (iip == NULL)
		return .0;
	*iip = defImgInfo;
	if (!pr.GetNAlloc())
		return .0;
	float	pixAspect = 1.f;
	int32	ival;
						// sample-to-nits factor
	if (PDBgetField(pr,PDBFstonits).Get(&iip->stonits))
		iip->flags |= IIFstonits;
						// image source
	if (PDBgetField(pr,PDBFsource).Get(&ival) &&
			PDBgetField(pr,PDBFmake).Get(iip->source.make,
				sizeof(iip->source.make)) &&
			PDBgetField(pr,PDBFmodel).Get(iip->source.model,
					sizeof(iip->source.model))) {
		iip->source.stype = ival;
		if (!PDBgetField(pr,PDBFversion).Get(iip->source.vers,
				sizeof(iip->source.vers)))
			iip->source.vers[0] = '\0';
		iip->flags |= IIFsource;
	}
						// capture date
	if (PDBgetField(pr,PDBFcapdate).Get(iip->capdate, sizeof(iip->capdate)))
		iip->flags |= IIFcapdate;
						// GMT date
	if (PDBgetField(pr,PDBFgmtdate).Get(iip->gmtdate, sizeof(iip->gmtdate)))
		iip->flags |= IIFgmtdate;
						// latitude/longitude
	if (PDBgetField(pr,PDBFlatitude).Get(&iip->latlong[0]) &&
			PDBgetField(pr,PDBFlongitude).Get(&iip->latlong[1]))
		iip->flags |= IIFlatlong;
						// altitude
	if (PDBgetField(pr,PDBFaltitude).Get(&iip->altitude))
		iip->flags |= IIFaltitude;
						// sight bearing
	if (PDBgetField(pr,PDBFbearing).Get(&iip->bearing))
		iip->flags |= IIFbearing;
						// view angle
	if (PDBgetField(pr,PDBFxangle).Get(&iip->hvangle))
		iip->flags |= IIFhvangle;
						// view specification
	if (PDBgetField(pr,PDBFview).Get(iip->view, sizeof(iip->view)))
		iip->flags |= IIFview;
						// compression quality
	if (PDBgetField(pr,PDBFquality).Get(&ival)) {
		iip->quality = ival;
		iip->flags |= IIFquality;
	}
	int32	rect[4];			// crop rectangle
	if (PDBgetField(pr,PDBFcrop).Get(rect,4) == 4) {
		iip->crop.xleft = rect[0];
		iip->crop.ytop = rect[1];
		iip->crop.xright = rect[2];
		iip->crop.ybottom = rect[3];
		iip->flags |= IIFcrop;
	}
	POrient	por(&pr);			// image orientation
	iip->orientation = IOtopleft;
	switch (por.xyswap << 2 | por.vflip << 1 | por.hflip) {
	case 01: iip->orientation = IOtopright; break;
	case 02: iip->orientation = IObotleft; break;
	case 03: iip->orientation = IObotright; break;
	case 04: iip->orientation = IOlefttop; break;
	case 05: iip->orientation = IOleftbot; break;
	case 06: iip->orientation = IOrighttop; break;
	case 07: iip->orientation = IOrightbot; break;
	}
	if (iip->orientation != IOtopleft)
		iip->flags |= IIForientation;
						// density
	if (PDBgetField(pr,PDBFxdensity).Get(&iip->hdensity)) {
		float	vdens;
		if (PDBgetField(pr,PDBFydensity).Get(&vdens))
			pixAspect = vdens/iip->hdensity;
		iip->flags |= IIFhdensity;
	}
						// focus distance
	if (PDBgetField(pr,PDBFfocus).Get(&iip->focus))
		iip->flags |= IIFfocus;
						// flash setting
	if (PDBgetField(pr,PDBFflash).Get(&ival)) {
		iip->flash = (ImgFlashMode)ival;
		iip->flags |= IIFflash;
	}
						// exposure seconds
	if (PDBgetField(pr,PDBFexptime).Get(&iip->exptime))
		iip->flags |= IIFexptime;
						// f-stop
	if (PDBgetField(pr,PDBFaperture).Get(&iip->aperture))
		iip->flags |= IIFaperture;
						// ISO ASA
	if (PDBgetField(pr,PDBFasa).Get(&iip->asa))
		iip->flags |= IIFasa;
						// white balance setting
	if (PDBgetField(pr,PDBFwhitebal).Get(&ival)) {
		iip->whitebal = (ImgWhiteBal)ival;
		iip->flags |= IIFwhitebal;
	}
						// content owner
	if (PDBgetField(pr,PDBFowner).Get(iip->owner, sizeof(iip->owner)))
		iip->flags |= IIFowner;
	char *		cp = iip->comments;	// history -> comments & params
	char *		pp = iip->params;
	const char *	hp;
	for (int i = 0; *(hp = PDBgetField(pr,PDBFhistory).GetString(i)); i++) {
		const int	hlen = strlen(hp);
		if (isalpha(hp[0]) && strchr(hp, '=') != NULL) {
						// parameter
			if (hlen >= iip->params+(sizeof(iip->params)-1) - pp)
				continue;
			while (*hp)
				*pp++ = *hp++;
			*pp++ = '\n';
			continue;
		}
						// else comment
		if (hlen >= iip->comments+(sizeof(iip->comments)-1) - cp)
			continue;
		while (*hp)
			*cp++ = *hp++;
		*cp++ = '\n';
	}
	if (pp > iip->params) {
		*pp = '\0';
		iip->flags |= IIFparams;
	}
	if (cp > iip->comments) {
		*cp = '\0';
		iip->flags |= IIFcomments;
	}
	return pixAspect;			// return pixel aspect ratio
}

// Set orientation transform according to record values
void
POrient::SetRecord(const DBRecord *ir)
{
	xres = yres = 1;
	hflip = vflip = false;
	xyswap = false;
	if (ir == NULL || !ir->GetNAlloc())
		return;
	PDBgetField(*ir,PDBFxsize).Get(&xres);
	PDBgetField(*ir,PDBFysize).Get(&yres);
	int32	iv;
	if (PDBgetField(*ir,PDBForient).Get(&iv))
		RotateCW(iv);
	if (PDBgetField(*ir,PDBFflip).Get(&iv)) {
		if (iv & 1) HFlip();
		if (iv & 2) VFlip();
	}
}

// Reorient image according to transform (color space and size must match)
bool
POrient::OrientImage(ImgStruct *dst, const ImgStruct *src) const
{
	if (dst == NULL || dst->csp == NULL)
		return false;
	if (src == NULL)
		src = dst;
	else if (src->csp == NULL)
		return false;
	if (src->img == NULL)
		return false;
	if ((src->xres != xres) | (src->yres != yres)) {
		DMESG(DMCparameter, "Mismatched source resolution in OrientImage");
		return false;
	}
	const int	psiz = ImgPixelSize(src->csp);
	const bool	vflipOnly = !(xyswap | hflip);
	int		xs, ys, xd, yd;
	const uby8 *	sptr;
					// check for different source & dest.
	if (dst != src) {
		if (!PmatchColorSpace(dst->csp, src->csp, PICMall)) {
			DMESG(DMCparameter, "Color space mismatch in OrientImage");
			PfreeImage(dst);
			return false;
		}
		if (dst->img != NULL) {
			if ((dst->xres != GetWidth()) | (dst->yres != GetHeight())) {
				DMESG(DMCparameter,
					"Size mismatch in OrientImage");
				PfreeImage(dst);
				return false;
			}
			if (PimagesOverlap(dst, src)) {
				DMESG(DMCparameter,
					"Cannot orient linked images in place");
				return false;
			}
		} else {
			dst->xres = GetWidth(); dst->yres = GetHeight();
		}
		if (IsNoop())
			return PmapImage(dst, src, 1.f);
		if (!PnewImage(dst, .0))
			return false;
		for (ys = 0; ys < yres; ys++) {
			sptr = ProwPtr(src,ys);
			if (vflipOnly) {
				memcpy(ProwPtr(dst,yres-1-ys), sptr, xres*psiz);
				continue;
			}
			for (xs = 0; xs < xres; xs++) {
				MapToCorrect(&xd, &yd, xs, ys);
				memcpy(PpixP(dst,xd,yd,psiz), sptr, psiz);
				sptr += psiz;
			}
		}
		return true;
	}
					// In-place copy (src == dst)
	if (IsNoop())
		return true;
	if (vflipOnly && dst->rowsize == -xres*psiz) {
		dst->img += (ssize_t)(yres-1)*dst->rowsize;
		dst->rowsize = -dst->rowsize;
		return true;
	}
	ImgStruct	idst;
	if (src->rowsize != xres*psiz) {
		idst.img = NULL;	// Need to allocate copy
		idst.csp = src->csp;
		if (!OrientImage(&idst, src))
			return false;
		dst->xres = idst.xres;
		dst->yres = idst.yres;
		dst->rowsize = dst->xres*psiz;
		PmapImage(dst, &idst, 1.f);
		PfreeImage(&idst);
		return true;
	}
	/*
	 * The following code copies pixels in natural loops,
	 * but requires that the source and destination images
	 * occupy exactly the same memory footprint.
	 */
	ABitMap2	pmark(GetWidth(), GetHeight());
	uby8		ptmp[6*sizeof(double)];
	DASSERT((int)sizeof(ptmp) >= psiz);
	idst = *dst;
	idst.xres = GetWidth(); idst.yres = GetHeight();
	idst.rowsize = idst.xres*psiz;
	for (xd = yd = 0; pmark.Find(&xd, &yd, false); xd++) {
					// start copy loop
		memcpy(ptmp, PpixP(&idst,xd,yd,psiz), psiz);
		pmark.Set(xd, yd);
		int	xd1 = xd, yd1 = yd;
		int	xd0, yd0;
		for ( ; ; ) {		// run full circle
			xd0 = xd1; yd0 = yd1;
			MapToOrig(&xs, &ys, xd0, yd0);
			sptr = PpixP(src,xs,ys,psiz);
			if (!PpixPos(&xd1, &yd1, &idst, sptr))
				DMESG(DMCassert, "PpixPos failed!");
			if (!pmark.TestAndSet(xd1, yd1))
				break;
			memcpy(PpixP(&idst,xd0,yd0,psiz), sptr, psiz);
		}
					// complete loop
		memcpy(PpixP(&idst,xd0,yd0,psiz), ptmp, psiz);
	}
	*dst = idst;
	return true;
}

// Close cached image
void
PCacheImage::CloseImage()
{
	isrc.csp = NULL;
	if (irdr == NULL)
		return;			// already closed
	DASSERT(!InUse());
	FreeMemory();			// free in-core image
	objSize = 0;
	if (rdrMine)
		IRclose(irdr);		// finish with reader
	irdr = NULL;
	ircd.Init();			// erase data record
}

// Set up cached image handler from image file
bool
PCacheImage::SetImage(const char *fpath, bool quiet)
{
	CloseImage();			// close old image
	if (fpath == NULL)
		return false;
	irdr = PopenImageF(fpath, quiet, NULL);
	if (irdr == NULL)
		return false;
	rdrMine = true;
	PgrokImage(&ircd, irdr);
	return true;
}

// Set up cached image from DB record
bool
PCacheImage::SetImage(const DBRecord &pr, bool quiet)
{
	CloseImage();			// close old image
	irdr = PopenImageD(pr, quiet);
	if (irdr == NULL)
		return false;
	rdrMine = true;
	ircd = pr;			// copy data record
	if (!PDBgetField(ircd,PDBFxsize).GetNV())
		PgrokImage(&ircd, irdr);
	return true;
}

// Set up cached image from image reader
bool
PCacheImage::SetImage(ImgReader *ir)
{
	CloseImage();			// close old image
	if (ir == NULL)
		return false;
	irdr = ir;
	rdrMine = false;
	PgrokImage(&ircd, irdr);
	return true;
}

// Get subimage region from cache image
bool
PCacheImage::GetSubimage(ImgStruct *ims, const ImgRect *r)
{
	if (ims == NULL || ims->csp == NULL || !Ready())
		return false;
	if (objCache == NULL) {			// use color space of new alloc
		if (ims->csp->cstatic)
			isrc.csp = ims->csp;
		else {
			PcopyCS(&myCS, ims->csp);
			isrc.csp = &myCS;
		}
	}
	if (GetCacheObject() == NULL)		// check/load image
		return false;
	if ((ims->mbase != NULL) & (ims->mbase == isrc.mbase))
		PfreeImage(ims);		// don't overwrite ourselves
	ImgRect		rdflt;
	if (r == NULL) {
		rdflt.xleft = 0; rdflt.xright = isrc.xres;
		rdflt.ytop = 0; rdflt.ybottom = isrc.yres;
		r = &rdflt;
	}
	if ((ims->xres <= 0) | (ims->yres <= 0)) {
		ims->xres = r->xright - r->xleft;
		ims->yres = r->ybottom - r->ytop;
	}
	bool	ok;				// possible link operation?
	if (ims->img == NULL && (ims->xres == r->xright - r->xleft) &
				(ims->yres == r->ybottom - r->ytop) &&
			PmatchColorSpace(ims->csp, isrc.csp, PICMall)) {
		ok = PlinkSubimage(ims, &isrc, r);
	} else {
		ImgStruct	isub;		// else render from subimage
		ok = PlinkSubimage(&isub, &isrc, r) &&
				PrenderImageI(ims, &isub);
		PfreeImage(&isub);
	}
	ReleaseObject(objCache);		// done with cache
	return ok;
}

// Load image into memory cache
bool
PCacheImage::RestoreMemory()
{
	if (IsResident())
		return true;
	DASSERT(irdr != NULL);
	DASSERT(isrc.csp != NULL);
	isrc.xres = irdr->xres; isrc.yres = irdr->yres;
	isrc.rowsize = isrc.xres*ImgPixelSize(isrc.csp);
	objSize = (size_t)isrc.yres*isrc.rowsize;
	objCache = Cmalloc(objSize);
	if (objCache == NULL) {
		objSize = 0;
		return false;
	}
	isrc.img = (uby8 *)objCache;
	isrc.mbase = NULL;
						// allocate and load image
	DMESGF(DMCtrace, "Loading image '%s'", irdr->file);
	if (!PrenderImageR(&isrc, irdr, false)) {
		Cfree(objCache); objCache = NULL;
		objSize = 0;
		return false;
	}
	isrc.mbase = this;			// only transfers retention
	return true;				// got it!
}

// Free image cache memory
void
PCacheImage::FreeMemory()
{
	if (!IsPurgeable())
		return;
	DASSERT(objCache == (void *)isrc.img);
	DMESGF(DMCtrace, "Unloading image '%s'", irdr->file);
	isrc.img = NULL;
	isrc.mbase = NULL;
	Cfree(objCache); objCache = NULL;
}
