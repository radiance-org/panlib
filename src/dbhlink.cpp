/*
 *  dbhlink.cpp
 *  panlib
 *
 *  Routines for reading and writing linked database header variables
 *
 *  Created by Greg Ward on Wed Apr 14 2004.
 *  Copyright (c) 2004 Anyhere Software. All rights reserved.
 *
 */

#include <stdlib.h>
#include <iostream>
using namespace std;
#include <sstream>
#include "rtio.h"
#include "pstrings.h"
#include "dbhlink.h"
#include "dmessage.h"

// Check to see if input matches the given string
bool
MatchInput(istream *istr, const char *txt)
{
	int	c;
	while (*txt && (c = istr->get()) != EOF)
		if (c != *txt++)
			return false;
	return (!*txt);
}

// Compute maximum header length for a list of variables
int
DBHeaderLink::GetMaxLength() const
{
	int	len = 0;
						// name + spaces + '='
	len += strlen(vname) + 12;
						// value
	switch (vtype) {
	case DBDTunknown:
	case DBDTstring:
		len += maxlen + 12;		// assumes one assignment
		break;
	case DBDTint:
		len = maxlen*(len + 14);	// account for multiple =
		break;
	case DBDTfloat:
		len = maxlen*(len + 16);	// account for multiple =
		break;
	}
						// other variables
	if (next != NULL)
		len += next->GetMaxLength();
						// return total
	return len;
}

// Set a string value using sequential storage (e.g., "this\0that\0")
int
DBHeaderLink::SetVal(char *sval, const char *sset)
{
	if ((vtype != DBDTstring) | (sval != link.sval)) {
		if (next == NULL)
			return 0;
		return next->SetVal(sval, sset);
	}
	const int	len1 = strlen(sset) + 1;
	if (len1 >= maxlen) {			// room for just this
		strlcpy(sval, sset, maxlen);
		return nset = 1;
	}
	const char *	svend = sval + maxlen;
	char *		lastsv = sval;
	char *		cp = sval;
	int		n;			// find value to set
	for (n = 0; n < nset; n++) {
		while (*cp++)			// next in sequence
			;
		if (svend-cp < len1)		// enough space left?
			break;
		lastsv = cp;
	}
	strcpy(lastsv, sset);			// copy value to position
	return nset = n + 1;			// new # set
}

// Read value from stream for this variable starting on non-white
int
DBHeaderLink::ReadVal(istream *istr)
{
	int	c;
						// other types (incl. opaque)
	char	linebuf[2048];			// read value to '\n'
	istr->getline(linebuf, sizeof(linebuf));
	if (istr->fail() || !linebuf[0]) {
		DMESGF(DMCdata, "Could not read value for %s",
				vname);
		return -1;
	}
	istringstream	iss(string(linebuf, strlen(linebuf)));
	switch (vtype) {
	case DBDTunknown:			// opaque type
		if (linebuf[0] != '{') {
			DMESGF(DMCdata, "Expected '{' for %s in settings",
					vname);
			return nset = -1;
		}
		if ((link.op.o == NULL) | (link.op.rdr == NULL)) {
			DMESGF(DMCparameter, "Missing opaque reader for %s",
					vname);
			return nset = -1;
		}
		if (nset < 0) nset = 0;
		c = (*link.op.rdr)(link.op.o, istr);	// call opaque reader
		if (c < 0)
			return nset = c;
		nset += c;			// should we do this?
		return c;
	case DBDTint:				// integer value(s)
		if (nset < 0) nset = 0;
		c = 0;
		while (!iss.eof()) {
			int	iv;
			iss >> iv;
			if (iss.fail())
				break;
			if (nset >= maxlen) nset--;
			link.ival[nset++] = iv;
			c++;
		}
		if (!c || !iss.eof()) {
			DMESGF(DMCdata, "Expected integer value for %s",
					vname);
			return nset = -1;
		}
		return c;
	case DBDTfloat:				// real value(s)
		if (nset < 0) nset = 0;
		c = 0;
		while (!iss.eof()) {
			float	fv;
			iss >> fv;
			if (iss.fail())
				break;
			if (nset >= maxlen) nset--;
			link.fval[nset++] = fv;
			c++;
		}
		if (!c || !iss.eof()) {
			DMESGF(DMCdata, "Expected float value for %s",
					vname);
			return nset = -1;
		}
		return c;
	case DBDTstring: {			// string value(s)
		int		n = 0;
		istream *	myistr = &iss;	// start with line buffer
		char		sbuf[4096];
		if (nset < 0) nset = 0;
		while (!myistr->eof()) {	// read to end of line
			char *	sp;
			c = myistr->get();	// check next character
			switch (c) {
			case '\n':
			case EOF:		// end of line or file
				break;
			case ' ':
			case '\t':
			case '\r':		// ignored white space
				continue;
			case '\'':
			case '"':		// quoted string
				sp = sbuf; *sp = '\0';
			quote_continue:
				myistr->get(sp, sizeof(sbuf)-(sp-sbuf), char(c));
				myistr->clear();
				if (myistr->get() != c) {
					if (myistr != istr) {
						// switch streams mid-boat
						myistr = istr;
						while (*sp) sp++;
						goto quote_continue;
					}
					DMESGF(DMCdata,
					"Missing quote for %s setting",
							vname);
					return nset = -1;
				}
				SetVal(link.sval, sbuf);
				n++;
				continue;
			default:		// unquoted word
				myistr->putback(c);
				*myistr >> sbuf;
				SetVal(link.sval, sbuf);
				n++;
				continue;
			}
			break;			// finished this setting
		}
		return n;			// return number read
		}
	}
	DMESG(DMCparameter, "Missing type in DBHeaderLink::CheckVar");
	return -1;
}

// Skip the next value in the given input stream
void
DBHeaderLink::SkipVal(istream *istr)
{
	int	c = istr->get();
	switch (c) {				// decide on first non-white
	case EOF:				// nothing
	case '\n':
		return;
	case '"':				// quoted string
	case '\'':
		istr->ignore(4096, c);
		break;
	case '{':				// opaque value
		while ((c = istr->get()) != EOF && c != '}')
			if ((c == '"') | (c == '\''))
				istr->ignore(4096, c);	// skip quoted string
			else if (c == '{') {		// nested opaque value?
				istr->putback(c);
				SkipVal(istr);
			}
		if (c != '}')				// found end brace
			DMESG(DMCwarning, "Expected '}' in settings");
		break;
	}
	istr->ignore(1024, '\n');	// eat rest of line
}

// Scan variable name from stream, digesting '=' and white space
char *
DBHeaderLink::ScanVarName(char *vname, int vnl, istream *istr)
{
	int	c, i;
	
	if ((vname == NULL) | (vnl < 2) | (istr == NULL))
		return NULL;
readloop:
	while ((c = istr->get()) != EOF)
		if (!isspace(c))		// skip leading white space
			break;
	switch (c) {
	case EOF:				// unexpected end-of-file
		DMESG(DMCwarning, "Unexpected EOF in settings");
		return NULL;
	case '}':				// end of header section
		istr->ignore(256, '\n');	// eat rest of line
		return NULL;
	case '=':
		DMESG(DMCdata, "Missing variable name in settings");
		return NULL;
	case '/':
	case '#':				// assume it's a comment
		istr->ignore(1024, '\n');	// so we ignore it
		goto readloop;
	default:
		if (!isalpha(c) && c != '_') {
			DMESGF(DMCdata,
				"Unexpected character ('%c') in settings", c);
			return NULL;
		}
		break;				// else part of variable name
	}
	vname[0] = c;				// read variable name
	for (i = 1; i < vnl; i++) {
		c = istr->get();
		if (c == EOF || isspace(c))
			break;
		if (c == '=') {
			istr->putback(c);
			break;
		}
		vname[i] = c;
	}
	if (i >= vnl) {
		DMESG(DMCparameter, "Variable name too long in settings");
		return NULL;
	}
	vname[i] = '\0';
						// scan in '='
	while ((c = istr->get()) != EOF)
		if ((c == '=') | (c == '\n'))
			break;
	if (c != '=') {
		DMESGF(DMCdata, "Missing '=' in assignment for %s",
				vname);
		return NULL;
	}
	while ((c = istr->get()) != EOF && c != '\n')
		if (!isspace(c))		// skip white space
			break;
	if (c == EOF) {
		DMESG(DMCdata, "Unexpected EOF in settings");
		return NULL;
	}
	istr->putback(c);			// replace first value character
	return vname;				// return loaded variable name
}

// Read header variable=value from stream, returning matched link
DBHeaderLink *
DBHeaderLink::ScanVar(istream *istr)
{
	DBHeaderLink *	hl = NULL;
	char		varname[128];
	int		c;
						// read until we get assignment
	for ( ; ; ) {
		if (ScanVarName(varname, sizeof(varname), istr) == NULL)
			return NULL;		// EOF or LHS syntax
		hl = FindVar(varname);		// find matching link
		if (hl != NULL) {
			c = hl->ReadVal(istr);	// read RHS
			if (c < 0)
				hl = NULL;	// RHS syntax
			break;
		}
		SkipVal(istr);			// skip RHS for stranger
		DMESGF(DMCwarning, "Unrecognized variable: %s", varname);
	}
	return hl;				// return assignment
}

// Find problem character for writing string, return EOF if none
static inline int
ProblemChar(const char *s)
{
	if (!*s)
		return 0;			// starting nul is a problem!
	int	pc = EOF;
	while (*++s)
		switch (*s) {
		case '\'':
		case '"':
			pc = *s;		// contains quote
		case ' ':
		case '\t':
		case '\r':
		case '\n':
			if (pc == EOF)
				pc = ' ';	// use space for any white
		}
	return pc;				// return worst found
}

// Write out database header variable=value for this link
int
DBHeaderLink::WriteVar(ostream *ostr, int tab) const
{
	int	i;
	if ((vname == NULL) | (nset <= 0))
		return 0;
	switch (vtype) {
	case DBDTunknown:			// write opaque data
		if ((link.op.o == NULL) | (link.op.wrt == NULL)) {
			DMESGF(DMCparameter, "Missing opaque writer for %s",
					vname);
			return 0;
		}
		TabIn(ostr, tab);
		*ostr << vname << " = {\n";
		i = (*link.op.wrt)(link.op.o, ostr, tab+1);
		TabIn(ostr, tab);
		*ostr << "} // " << vname << '\n';
		return i;
	case DBDTint:				// write integer(s)
		TabIn(ostr, tab);
		*ostr << vname << " =";
		for (i = 0; i < nset; i++)
			*ostr << ' ' << link.ival[i];
		*ostr << '\n';
		break;
	case DBDTfloat:				// write float(s)
		TabIn(ostr, tab);
		*ostr << vname << " =";
		for (i = 0; i < nset; i++)
			*ostr << ' ' << link.fval[i];
		*ostr << '\n';
		break;
	case DBDTstring: {			// write string(s)
		if (link.sval == NULL)
			return 0;
		const char *	sv = link.sval;
		for (i = 0; i < nset; i++) {
			TabIn(ostr, tab);
			*ostr << vname << " = ";
			switch (ProblemChar(sv)) {
			case EOF:		// no problem characters
				*ostr << sv << '\n';
				break;
			case '"':		// contains double-quote
				*ostr << '\'' << sv << "'\n";
				break;
			default:		// whitespace or single-quote
				*ostr << '"' << sv << "\"\n";
				break;
			}
			while (*sv++)		// next in sequence
				;
		}
		} break;
	default:				// botched type!
		DMESG(DMCparameter, "Missing type in DBHeaderLink::WriteVar");
		return -1;
	}
	return ostr->bad() ? -1 : nset;
}

// Insert string value at the specified point in list
bool
DBHeaderLink::InsertSVal(const char *val, int n)
{
	if (val == NULL)
		return false;
	if (vtype != DBDTstring)
		return false;
	if (n > nset)
		return false;
	if (nset <= 0) {
		if ((int)strlen(val) >= maxlen)
			return false;
		strcpy(link.sval, val);
		nset = 1;
		return true;
	}
	int		i = nset;
	char *		cp = link.sval;
	char *		cnew = NULL;
	while (i--) {
		if (!n--)
			cnew = cp;
		while (*cp++)
			;
	}
	if (cnew == NULL) {
		if (cp - link.sval + (int)strlen(val) >= maxlen)
			return false;
		strcpy(cp, val);
		++nset;
		return true;
	}
	const char *	cp1 = cp;
	cp += strlen(val) + 1;
	if (cp - link.sval > maxlen)
		return false;
	while (cp1 > cnew)
		*--cp = *--cp1;
	strcpy(cnew, val);
	++nset;
	return true;
}

// Delete the indicated string value (-1 == last)
void
DBHeaderLink::DeleteSVal(int n)
{
	if ((vtype != DBDTstring) | (nset <= 0))
		return;
	if (n >= nset)
		return;
	--nset;
	if (n < 0)
		return;
	int	nafter = nset - n;
	if (!nafter)
		return;
	char *	cp0 = link.sval;
	while (n--)
		while (*cp0++)
			;
	const char *	cp1 = cp0;
	while (*cp1++)
		;
	while (nafter--)
		while ((*cp0++ = *cp1++))
			;
}
