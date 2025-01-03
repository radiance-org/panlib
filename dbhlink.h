/*
 *  dbhlink.h
 *  panlib
 *
 *  Include after <iostream.h> and "pstrings.h"
 *
 *  Class declarations for linked database header variables
 *
 *  Created by Greg Ward on Wed Apr 14 2004.
 *  Copyright (c) 2004 Anyhere Software. All rights reserved.
 *
 */

#ifndef _DBHLINK_H_
#define _DBHLINK_H_

// Data types:  unknown, 32-bit integer, 32-bit float, nul-terminated string
enum	DBDataType	{ DBDTunknown, DBDTint, DBDTfloat, DBDTstring };

// Check if input matches the given string
bool    MatchInput(istream *istr, const char *txt);

// Tab in the specified number of tabs
inline void
TabIn(ostream *ostr, int tab = 1)
{
	while (tab-- > 0)
		ostr->put('\t');
}

// Header key-value handler
// Use this class to read in header variable lines of the form:
//	Var1 = val0 val1 val2 \n
//	Var2 = val0 \n
//	Var1 = val3 \n
//	...
// where all values for a particular variable are of one type.
// Multiple string values are stored sequentially in memory
// (e.g., "this\0that\0"), and may be quoted to contain spaces
// or newlines.  Otherwise,  one string will be assigned per word.
// More complex assignments can be handled with the syntax:
//	Var = {
//		(anything...)
//	} (ignore to newline)
// Everything up to and including the line with '}' is handled by a
// a specified function pointer, and the given function may,
// if desired, use another DBHeaderLink list to read the contents.
// Variables that don't match anyone are skipped (including {}'s).
// Lines whose first non-white character is '#' or '/' are ignored.
// Reading stops after finding a line with an unmatched '}', which
// presumably matches the opening '{' in the header's beginning
// "magic string" (not read by us)
struct DBHeaderLink {
	const char *	vname;			// variable name in header
	DBDataType	vtype;			// variable type
	int		maxlen;			// array/string length
	int		nset;			// # values assigned
	union {
		int *		ival;			// integer value(s)
		float *		fval;			// float value(s)
		char *		sval;			// string value(s)
		struct {
			void *	o;			// header object
					// write opaque value(s) ending with '}'
			int		(*wrt)(void *, ostream *, int);
					// read opaque value(s) ending with '}'
			int		(*rdr)(void *, istream *);
		}		op;			// opaque value handler
	}		link;			// value handle(r)
	DBHeaderLink *	next;			// next in list
			// variable names aren't copied, so best if static
			DBHeaderLink(const char *vn, int *iv, int len,
					DBHeaderLink *np = 0) {
				vname = vn; vtype = DBDTint; maxlen = len;
				link.ival = iv; nset = 0; next = np;
			}
			DBHeaderLink(const char *vn, float *fv, int len,
					DBHeaderLink *np = 0) {
				vname = vn; vtype = DBDTfloat; maxlen = len;
				link.fval = fv; nset = 0; next = np;
			}
			DBHeaderLink(const char *vn, char *sv, int len,
					DBHeaderLink *np = 0) {
				vname = vn; vtype = DBDTstring; maxlen = len;
				link.sval = sv; nset = 0; next = np;
			}
			DBHeaderLink(const char *vn, void *ho,
					int (*wo)(void *, ostream *, int),
					int (*ro)(void *, istream *), int len,
					DBHeaderLink *np = 0) {
				vname = vn; vtype = DBDTunknown; maxlen = len;
				link.op.o = ho; link.op.wrt = wo;
				link.op.rdr = ro;
				nset = 0; next = np;
			}
			~DBHeaderLink() {
				delete next;
			}
			// Clear all variable assignment counts
	void		Clear() {
				if (next != 0) next->Clear();
				nset = 0;
			}
			// Typed setting routines
	int		SetVal(int *ival, int iset) {
				if ((vtype == DBDTint & ival == link.ival)) {
					if (nset >= maxlen) nset--;
					ival[nset++] = iset;
					return nset;
				}
				if (!next) return 0;
				return next->SetVal(ival, iset);
			}
	int		SetVal(float *fval, float fset) {
				if ((vtype == DBDTfloat & fval == link.fval)) {
					if (nset >= maxlen) nset--;
					fval[nset++] = fset;
					return nset;
				}
				if (!next) return 0;
				return next->SetVal(fval, fset);
			}
	int		SetVal(char *sval, const char *sset);
			// Insert string value at the specified point in list
	bool		InsertSVal(const char *val, int n = 0);
			// Delete the indicated string value (-1 == last)
	void		DeleteSVal(int n = -1);
			// Get nth string value in a sequence
	static const char *
			GetSVal(const char *sp, int n) {
				if ((n < 0) | !sp) return 0;
				while (n--) while (*sp++) ;
				return sp;
			}
			// Get nth string value in this sequence
	const char *	GetSVal(int n = 0) const {
				if (!this || vtype != DBDTstring) return 0;
				if (n >= nset) return 0;
				return GetSVal(link.sval, n);
			}
			// Get single (final) string setting
	const char *	GetFinalSVal(const char *sval = 0) const {
				if (!this) return 0;
				if (!sval) {
					if (vtype != DBDTstring) return 0;
					sval = link.sval;
				}
				if ((vtype == DBDTstring) & (sval == link.sval)) {
					if (nset <= 0) return 0;
					return GetSVal(sval, nset-1);
				}
				if (!next) return 0;
				return next->GetFinalSVal(sval);
			}
			// Find index of matching string (-1 if none)
	int		FindSVal(const char *val) const {
				if (!val) return -1;
				if (!this || vtype != DBDTstring) return -1;
				const char *	sp = link.sval;
				for (int i = 0; i < nset; i++) {
					if (!strcmp(sp, val)) return i;
					while (*sp++) ;
				}
				return -1;
			}
			// Compute maximum header length in bytes
	int		GetMaxLength() const;
			// Find a variable link by its name
	DBHeaderLink *	FindVar(const char *vn) {
				if (!this) return 0;
				if (!istrcmp(vn, vname)) return this;
				return next->FindVar(vn);
			}
			// Const version of FindVar
	const DBHeaderLink *
			GetVar(const char *vn) const {
				if (!this) return 0;
				if (!istrcmp(vn, vname)) return this;
				return next->GetVar(vn);
			}
			// Scan next variable assignment name (eats '=')
	static char *	ScanVarName(char *vname, int vnl, istream *istr);
			// Scan for next header link & return match
	DBHeaderLink *	ScanVar(istream *istr);
			// Read value from stream for this variable
	int		ReadVal(istream *istr);
			// Skip the next value in the input stream
	static void	SkipVal(istream *istr);
			// Get total values currently set or -1 if error
	int		GetNSet() const {
				if (nset < 0) return -1;
				if (!next) return nset;
				int	ns = next->GetNSet();
				if (ns < 0) return -1;
				return ns + nset;
			}
			// Keep reading until an unmatched '}'
	int		ScanVars(istream *istr) {
				Clear();
				while (ScanVar(istr) != 0)
					;
				return GetNSet();
			}
			// Write this link if assigned, negative on error
	int		WriteVar(ostream *ostr, int tab = 0) const;
			// Output assigned header variable(s)
	int		WriteVars(ostream *ostr, int tab = 0) const {
				int	nw = 0;
				if (next != 0)
					nw = next->WriteVars(ostr, tab);
				if (nw < 0) return nw;
				int	mynw = WriteVar(ostr, tab);
				if (mynw < 0) return mynw;
				return nw + mynw;
			}
};

#endif  // ! _DBHLINK_H_
