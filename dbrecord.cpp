/*
 *  dbrecord.cpp
 *  panlib
 *
 *  Database record class implementations.
 *
 *  Created by gward on June 12, 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */
 
#include <math.h>
#include <ctype.h>
#include "system.h"
#include "rtio.h"
#include "dbase.h"
#include "dmessage.h"

#ifndef DB_MAXARR
#define DB_MAXARR	512			// biggest array we handle
#endif

DBRecord	DBRecordList::dummy;

// Test two reals for (near) equality
template <class real>
static inline bool
RealClose(real a, real b, real eps=real(1e-6))
{
	if (isnan(a))
		return isnan(b);
	if (isinf(a))
		return isinf(b);
	if (b != real(0))
		a = a/b - real(1);
	return ((-eps <= a) & (a <= eps));
}

// Compute hash on field value
// WARNING: Result depends on machine byte ordering!
unsigned long
DBField::Hash(int nbits) const
{
	if (!nv)
		return 0L;
	if (nv == 1)			// single value
		switch (dt) {
		case DBDTint:
			return memhash(&v.mint, sizeof(int32), nbits);
		case DBDTfloat:
			return memhash(&v.mfloat, sizeof(float), nbits);
		}
					// array or string value
	const void *	ptr = GetPointer();
	switch (dt) {
	case DBDTint:
		return memhash(ptr, nv*sizeof(int32), nbits);
	case DBDTfloat:
		return memhash(ptr, nv*sizeof(float), nbits);
	case DBDTstring:
		if (nv == 1)
			return strhash_nc((const char *)ptr, nbits);
		if (fl & DBFFstrseq) {
			unsigned long	hval = 0;
			const char *	cp = (const char *)ptr;
			int		i = nv;
			while (i--) {
				hval ^= strhash_nc(cp, nbits);
				while (*cp++)
					;
			}
			return hval;
		} else {
			unsigned long	hval = 0;
			const char **	sa = (const char **)ptr + nv;
			while (sa-- > (const char **)ptr)
				hval ^= strhash_nc(*sa, nbits);
			return hval;
		}
	default:
		DMESG(DMCparameter, "Botched type in DBField::GetHash");
		break;
	}
	return 0L;			// error return
}

// Clear a field object, freeing memory if allocated
void
DBField::Clear()
{
						// check for no allocation
	if (!(fl & DBFFownmem) || !HasPointer())
		goto clean;
						// check for string (or sequence)
	if (dt == DBDTstring && (nv == 1 || fl & DBFFstrseq)) {
		if (fl & DBFFoffset)		// XXX should never happen
			free(((char *)this + v.offset));
		else if (nv == 1)
			strunlink(v.str);
		else
			free(const_cast<void *>(v.ptr));
		goto clean;
	}
	if (fl & DBFFoffset)			// XXX should never happen
		v.ptr = (void *)((char *)this + v.offset);

	switch (dt) {				// free appropriate type
	case DBDTint:
		delete [] v.aint;
		break;
	case DBDTfloat:
		delete [] v.afloat;
		break;
	case DBDTstring: {
		int	i = nv;
		while (i--) strunlink(v.astr[i]);
		delete [] v.astr;
		} break;
	default:
		DMESG(DMCparameter, "Botched type in DBField::Clear");
		break;
	}
clean:						// clear the rest
	fl = 0; nv = 0; dt = DBDTunknown; v.ptr = NULL;
}

// Get generic pointer to ith data item
const void *
DBField::GetPointer(int i) const
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
		return (const void *)((const char **)cp)[i];
	default:
		DMESG(DMCparameter, "Botched type in DBField::Get");
		break;
	}
	return NULL;
}

// Copy ith string value
bool
DBField::Get(char str[], int len, int i) const
{
	if ((str == NULL) | (len < 2) | (dt != DBDTstring))
		return false;
	const char *	sp = (const char *)GetPointer(i);
	if (sp == NULL)
		return false;
	strlcpy(str, sp, len);
	return true;
}

// Get one or more string values, sharing memory
int
DBField::Get(const char *sa[], int n) const
{
	if (!nv | (n <= 0) || dt != DBDTstring)
		return 0;
	if (n > nv)
		n = nv;
	const char *	cp = (const char *)GetPointer();
	if (nv == 1)
		sa[0] = cp;
	else if (fl & DBFFstrseq) {
		int			i = n;
		while (i--) {
			*sa++ = cp;
			while (*cp++)
				;
		}
	} else
		memcpy(sa, cp, n*sizeof(const char *));
	return n;
}

// Copy ith value to a singular field
bool
DBField::Get(DBField *fv, int i) const
{
	if ((fv == NULL) | (i < 0) | (i >= nv))
		return false;
	switch (dt) {
	case DBDTint:
		fv->Set(GetInt(i));
		return true;
	case DBDTfloat:
		fv->Set(GetFloat(i));
		return true;
	case DBDTstring:
		fv->Set(GetString(i));
		return true;
	default:
		DMESG(DMCparameter, "Botched type in DBField::Get");
		break;
	}
	return false;
}

// Get ith integer, casting type and filling zeroes as necessary
int32
DBField::GetInt(int i) const
{
	const void *	p = GetPointer(i);
	if (p == NULL)
		return 0;
	switch (dt) {
	case DBDTint:
		return *(const int32 *)p;
	case DBDTfloat:
		return int(*(const float *)p + .5f) -
				(*(const float *)p <= -.5f);
	case DBDTstring:
		return (int32)atol((const char *)p);
	default:
		DMESG(DMCparameter, "Botched type in DBField::GetInt");
		break;
	}
	return 0;
}

// Get ith float, casting type and filling zeroes as necessary
float
DBField::GetFloat(int i) const
{
	const void *	p = GetPointer(i);
	if (p == NULL)
		return .0f;
	switch (dt) {
	case DBDTfloat:
		return *(const float *)p;
	case DBDTint:
		return float(*(const int32 *)p);
	case DBDTstring:
		return (float)atof((const char *)p);
	default:
		DMESG(DMCparameter, "Botched type in DBField::GetFloat");
		break;
	}
	return .0f;
}

// Get ith string, casting type and filling empties as necessary
const char *
DBField::GetString(int i) const
{
	static char	cvbuf[32][20];		// XXX thread tolerant
	static int	freeb = 0;
	int		myb;
	const void *	p = GetPointer(i);
	if (p == NULL)
		return "";
	switch (dt) {
	case DBDTstring:
		return (const char *)p;
	case DBDTint:
		freeb *= (freeb < 32);
		myb = freeb++;
		sprintf(cvbuf[myb], "%ld", long(*(const int32 *)p));
		return cvbuf[myb];
	case DBDTfloat:
		freeb *= (freeb < 32);
		myb = freeb++;
		sprintf(cvbuf[myb], "%.7e", *(const float *)p);
		return cvbuf[myb];
	default:
		DMESG(DMCparameter, "Botched type in DBField::GetString");
		break;
	}
	return "";
}

// Set an array of integer values
void
DBField::Set(const int32 ia[], int n)
{
	if (ia == NULL)
		return;
	if (n == 1) {
		Set(ia[0]);
		return;
	}
	Clear();
	if (n <= 0)
		return;
	fl = DBFFarray|DBFFownmem;
	dt = DBDTint;
	nv = n;
	v.aint = new int32 [nv];
	memcpy(v.aint, ia, n*sizeof(int32));
}

// Assign a single floating-point value
void
DBField::Set(float fv)
{
	Clear();
	fl = 0; dt = DBDTfloat; nv = 1;
	if (isnan(fv) | isinf(fv)) {
		v.mfloat = 0;
		DMESG(DMCwarning, "Bad float zeroed in DBField::Set()");
	} else
		v.mfloat = fv;
}

// Set an array of float values
void
DBField::Set(const float fa[], int n)
{
	if (fa == NULL)
		return;
	if (n == 1) {
		Set(fa[0]);
		return;
	}
	Clear();
	if (n <= 0)
		return;
	fl = DBFFarray|DBFFownmem;
	dt = DBDTfloat;
	nv = n;
	v.afloat = new float [nv];
	memcpy(v.afloat, fa, n*sizeof(float));
	int	badfloats = 0;
	for (int i = nv; i--; )
		if (isnan(v.afloat[i]) | isinf(v.afloat[i])) {
			v.afloat[i] = 0;
			++badfloats;
		}
	if (badfloats)
		DMESGF(DMCwarning, "%d bad float(s) zeroed in DBField::Set()", badfloats);
}

// Set an array of strings, copying memory
void
DBField::Set(const char * const sa[], int n)
{
	if (sa == NULL)
		return;
	if (n == 1) {
		Set(sa[0]);
		return;
	}
	Clear();
	if (n <= 0)
		return;
	fl = DBFFarray|DBFFownmem;
	dt = DBDTstring;
	nv = n;
	v.astr = new const char * [nv];
	for (int i = nv; i--; )
		v.astr[i] = strlink(sa[i]);
}

// Put line(s) into string array field (newline separates elements)
int
DBField::SetText(const char *lines, const char nl)
{
	if (lines == NULL || !*lines) {		// special case
		Clear();
		return 0;
	}
	char *		buf = new char [strlen(lines)+1];
	char *		sarr[DB_MAXARR];
	int		n;
	char *		cp;
	strcpy(buf, lines);			// copy lines
	for (sarr[0] = cp = buf, n = 1; *cp && n < DB_MAXARR; cp++)
		if (*cp == nl) {
			*cp = '\0';
			if (cp[1])
				sarr[n++] = cp + 1;
		}
	Set(sarr, n);				// set array
	delete [] buf;
	return n;
}

static int
cmp_int32(const void *ip1, const void *ip2)
{
	return *(const int32 *)ip1 - *(const int32 *)ip2;
}

static int
cmp_float(const void *fp1, const void *fp2)
{
	float	diff = *(const float *)fp1 - *(const float *)fp2;
	return (diff > 0 ? 1 : diff < 0 ? -1 : 0);
}

static int
cmp_string(const void *sp1, const void *sp2)
{
	return istrcmp(*(const char * const *)sp1, *(const char * const *)sp2);
}

// Sort array values if we can
bool
DBField::SortArray()
{
						// check for no array
	if (nv <= 1)
		return (bool)nv;
						// make sure array is sortable
	if ((fl & (DBFFownmem|DBFFstrseq)) != DBFFownmem) {
		DBField		fo;
		fo.dt = dt; fo.fl = fl;
		fo.nv = nv; fo.v = v;
		nv = 0; v.ptr = NULL;
		*this = fo;
	} else if (fl & DBFFoffset) {
		v.ptr = (void *)((char *)this + v.offset);
		fl &= ~DBFFoffset;
	}
						// call qsort()
	switch (dt) {
	case DBDTint:
		qsort(v.aint, nv, sizeof(int32), cmp_int32);
		break;
	case DBDTfloat:
		qsort(v.afloat, nv, sizeof(float), cmp_float);
		break;
	case DBDTstring:
		qsort(v.astr, nv, sizeof(char *), cmp_string);
		break;
	default:
		DMESG(DMCparameter, "Botched type in DBField::SortArray");
		return false;
	}
	return true;
}

// Field copy operator (duplicates memory for arrays)
DBField &
DBField::operator=(const DBField &src)
{
	if (this == &src)
		return *this;
	Clear();				// clear previous value
	if (!src.nv)
		return *this;
	fl = src.fl & ~DBFFoffset;
	dt = src.dt;
	nv = src.nv;
	v = src.v;
						// check for memory use
	const void *	ptr = src.GetPointer();
	if (ptr == NULL)
		return *this;
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
		DMESG(DMCparameter, "Botched type in DBField copy");
		v.ptr = NULL;
		return *this;
	}
	fl |= DBFFownmem;
	return *this;
}

// Field reference operator (does not allocate arrays)
DBField &
DBField::operator=(const DBField *ref)
{
	if (this == ref)
		return *this;
	Clear();				// clear previous value
	if (!ref->nv)
		return *this;
	fl = ref->fl & ~(DBFFoffset|DBFFownmem);
	dt = ref->dt;
	nv = ref->nv;
	if ((v.ptr = ref->GetPointer()) == NULL)
		switch (dt) {
		case DBDTint:
			v.mint = ref->v.mint;
			break;
		case DBDTfloat:
			v.mfloat = ref->v.mfloat;
			break;
		default:
			DMESG(DMCparameter, "Botched type in DBField reference");
		}
	return *this;				// all done
}

// Field comparison function (arrays compare until mismatch)
int
DBField::Compare(const DBField &that) const
{
	if (this == &that)
		return 0;
	if (!that.nv)
		return nv;			// empty is "smaller"
	if (!nv)
		return -1;
	if (dt != that.dt)
		return dt - that.dt;		// sort by type(?!)
	const int	n = (nv < that.nv) ? nv : that.nv;
	int		i;
	if (dt == DBDTstring) {			// string compare
		const char	*sa1[DB_MAXARR], *sa2[DB_MAXARR];
		int		cmp;
		if (n > DB_MAXARR)
			DMESG(DMCparameter, "Too many strings in field compare");
		Get(sa1, n);
		that.Get(sa2, n);
		for (i = 0; i < n; i++)
			if ((cmp = strcmp_nc(sa1[i], sa2[i])))
				return cmp;
		return nv - that.nv;
	}
	const void *	ptr1 = GetPointer(0);
	const void *	ptr2 = that.GetPointer(0);
	if (ptr1 == ptr2)
		return 0;
	if (ptr2 == NULL)
		return 1;
	if (ptr1 == NULL)
		return -1;
	switch (dt) {
	case DBDTint: {
		const int32 *	ia1 = (const int32 *)ptr1;
		const int32 *	ia2 = (const int32 *)ptr2;
		for (i = 0; i < nv; i++)
			if (ia1[i] != ia2[i])
				return ia1[i] - ia2[i];
		} return nv - that.nv;
	case DBDTfloat: {
		const float *	fa1 = (const float *)ptr1;
		const float *	fa2 = (const float *)ptr2;
		for (i = 0; i < nv; i++)
			if (!RealClose(fa1[i], fa2[i]))
				return (fa1[i] > fa2[i]) ? 1 : -1;
		} return nv - that.nv;
	}
	DMESG(DMCparameter, "Unknown data type in field compare");
	return 0;
}

// check fields for equality (all array values must match exactly)
bool
DBField::operator==(const DBField &that) const
{
	if (this == &that)			// parameter matching
		return true;
	if (nv != that.nv)
		return false;
	if (!nv)
		return true;
	if (dt != that.dt)
		return false;
	int	i;
	if (dt == DBDTstring) {			// string compare
		const char	*sa1[DB_MAXARR], *sa2[DB_MAXARR];
		i = Get(sa1, DB_MAXARR);
		if (i < nv)
			DMESG(DMCparameter, "Too many strings in field compare");
		that.Get(sa2, DB_MAXARR);
		while (i--)
			if (istrcmp(sa1[i], sa2[i]))
				return false;
		return true;
	}
	if (nv == 1)				// single value
		switch (dt) {
		case DBDTint:
			return v.mint == that.v.mint;
		case DBDTfloat:
			return RealClose(v.mfloat, that.v.mfloat);
		}
						// array value
	switch (dt) {
	case DBDTint: {
		const int32 *	ia1 = (const int32 *)GetPointer();
		const int32 *	ia2 = (const int32 *)that.GetPointer();
		if (ia1 == ia2) return true;
		for (i = nv; i--; )
			if (*ia1++ != *ia2++)
				return false;
		} return true;
	case DBDTfloat: {
		const float *	fa1 = (const float *)GetPointer();
		const float *	fa2 = (const float *)that.GetPointer();
		if (fa1 == fa2) return true;
		for (i = nv; i--; )
			if (!RealClose(*fa1++, *fa2++))
				return false;
		} return true;
	}
	DMESG(DMCparameter, "Unknown data type in field compare");
	return false;
}

// Check for range match in integer field
bool
DBField::MatchRange(int32 ilo, int32 ihi) const
{
	if (!nv || dt != DBDTint)
		return false;
	if (nv == 1)
		return ((ilo <= v.mint) & (v.mint <= ihi));

	const int32 *	ip = (const int32 *)GetPointer();
	DASSERT(ip != NULL);
	for (int i = nv; i--; ip++)
		if (ilo <= *ip && *ip <= ihi)
			return true;
	return false;
}

// Check for range match in real field
bool
DBField::MatchRange(float flo, float fhi) const
{
	if (!nv || dt != DBDTfloat)
		return false;
	if (RealClose(flo, fhi)) {
		float	eps = 1e-6*fabs(flo+fhi);
		if (eps == 0) eps = 1e-16f;
		flo -= eps;
		fhi += eps;
	}
	if (nv == 1)
		return ((flo <= v.mfloat) & (v.mfloat <= fhi));

	const float *	fp = (const float *)GetPointer();
	DASSERT(fp != NULL);
	for (int i = nv; i--; fp++)
		if (flo <= *fp && *fp <= fhi)
			return true;
	return false;
}

// Check if the given string lies between the others
static inline bool
betweenStrings(const char *slo, const char *s, const char *shi)
{
	int	lores = strcmp_nc(slo, s);
	if (lores > 0)
		return false;
	if (shi == slo)
		return (lores == 0);
	return (strcmp_nc(s, shi) <= 0);
	
}

// Check for range match in string field (caseless string comparison)
bool
DBField::MatchRange(const char *slo, const char *shi) const
{
	if (!nv || dt != DBDTstring || (slo == NULL) | (shi == NULL))
		return false;
	int	i;
	if (nv == 1 || (fl & DBFFstrseq)) {
		const char *	cp = (const char *)GetPointer();
		DASSERT(cp != NULL);
		for (i = nv; i--; ) {
			if (betweenStrings(slo, cp, shi))
				return true;
			while (*cp++)
				;
		}
		return false;
	}
	const char * const *	cpp = (const char * const *)GetPointer();
	for (i = nv; i--; cpp++)
		if (*cpp != NULL && betweenStrings(slo, *cpp, shi))
			return true;
	return false;
}

// Check for range match in arbitrary field
// Array field must match all lo-hi pairs
bool
DBField::MatchRange(const DBField &lo, const DBField &hi) const
{
	if (!nv | !lo.nv | (lo.nv != hi.nv))
		return false;
	if ((dt != lo.dt) | (dt != hi.dt))
		return false;
	switch (dt) {				// match according to type
	case DBDTint:
		if (lo.nv == 1)
			return MatchRange(lo.v.mint, hi.v.mint);
		{
			const int32 *	ia1 = (const int32 *)lo.GetPointer();
			const int32 *	ia2 = (const int32 *)hi.GetPointer();
			int		i = lo.nv;
			while (i--)
				if (!MatchRange(ia1[i], ia2[i]))
					return false;
			return true;
		}
	case DBDTfloat:
		if (lo.nv == 1)
			return MatchRange(lo.v.mfloat, hi.v.mfloat);
		{
			const float *	fa1 = (const float *)lo.GetPointer();
			const float *	fa2 = (const float *)hi.GetPointer();
			int		i = lo.nv;
			while (i--)
				if (!MatchRange(fa1[i], fa2[i]))
					return false;
			return true;
		}
	case DBDTstring:
		{
			const char	*sa1[DB_MAXARR], *sa2[DB_MAXARR];
			int		i = lo.Get(sa1, DB_MAXARR);
			hi.Get(sa2, DB_MAXARR);
			while (i--)
				if (!MatchRange(sa1[i], sa2[i]))
					return false;
			return true;
		}
	}
	DMESG(DMCparameter, "Unknown data type in range match");
	return false;
}

// Find common members of integer array
static int
arrayMatches(int32 *ia0, int n0, ABitMap *bm0, int32 *ia1, int n1)
{
	if ((n0 <= 0) | (n1 <= 0))
		return 0;
	int	i, j;
	bm0->NewBitMap(n0);
	for (i = n0; i--; )
		for (j = n1; j--; )
			if (ia0[i] == ia1[j])
				bm0->Set(i);
	return bm0->SumTotal();
}

// Find common members of float array
static int
arrayMatches(float *fa0, int n0, ABitMap *bm0, float *fa1, int n1)
{
	if ((n0 <= 0) | (n1 <= 0))
		return 0;
	int	i, j;
	bm0->NewBitMap(n0);
	for (i = n0; i--; )
		for (j = n1; j--; )
			if (RealClose(fa0[i], fa1[j]))
				bm0->Set(i);
	return bm0->SumTotal();
}

// Find common members of string array
static int
arrayMatches(const char **sa0, int n0, ABitMap *bm0, const char **sa1, int n1)
{
	if ((n0 <= 0) | (n1 <= 0))
		return 0;
	int	i, j;
	bm0->NewBitMap(n0);
	for (i = n0; i--; )
		for (j = n1; j--; )
			if (!istrcmp(sa0[i], sa1[j]))
				bm0->Set(i);
	return bm0->SumTotal();
}

// Append array values from another field
DBField &
operator+=(DBField &dst, const DBField &src)
{
	if (src.GetNV() <= 0)
		return dst;
	if (dst.GetNV() <= 0)
		return dst = src;
	if (src.GetDT() != dst.GetDT()) {
		DMESG(DMCparameter, "Mismatched types in DBField::operator+=()");
		return dst;
	}
	int		nv;
	switch (dst.GetDT()) {
	case DBDTint:
		{
			int32		iv[DB_MAXARR];
			nv = dst.Get(iv, DB_MAXARR);
			nv += src.Get(iv+nv, DB_MAXARR-nv);
			dst.Set(iv, nv);
		}
		break;
	case DBDTfloat:
		{
			float		fv[DB_MAXARR];
			nv = dst.Get(fv, DB_MAXARR);
			nv += src.Get(fv+nv, DB_MAXARR-nv);
			dst.Set(fv, nv);
		}
		break;
	case DBDTstring:
		{
			DBField		dstdup(dst);
			const char 	*sv[DB_MAXARR];
			nv = dstdup.Get(sv, DB_MAXARR);
			nv += src.Get(sv+nv, DB_MAXARR-nv);
			dst.Set(sv, nv);
		}
		break;
	default:
		break;
	}
	return dst;
}

// Append different array values from another field
DBField &
operator|=(DBField &dst, const DBField &src)
{
	if (src.GetNV() <= 0)
		return dst;
	if (dst.GetNV() <= 0)
		return dst = src;
	if (src.GetDT() != dst.GetDT()) {
		DMESG(DMCparameter, "Mismatched types in DBField::operator|=()");
		return dst;
	}
	int		nv0, nv1;
	ABitMap		com1;
	switch (dst.GetDT()) {
	case DBDTint:
		{
			int32		iv0[DB_MAXARR], iv1[DB_MAXARR];
			nv0 = dst.Get(iv0, DB_MAXARR);
			nv1 = src.Get(iv1, DB_MAXARR);
			if (arrayMatches(iv1, nv1, &com1, iv0, nv0)) {
				for (uint32 i = 0; com1.Find(&i,false); i++) {
					if (nv0 >= DB_MAXARR) break;
					iv0[nv0++] = iv1[i];
				}
			} else {
				if (nv0 + nv1 > DB_MAXARR)
					nv1 = DB_MAXARR - nv0;
				memcpy(iv0+nv0, iv1, nv1*sizeof(int32));
				nv0 += nv1;
			}
			dst.Set(iv0, nv0);
		}
		break;
	case DBDTfloat:
		{
			float		fv0[DB_MAXARR], fv1[DB_MAXARR];
			nv0 = dst.Get(fv0, DB_MAXARR);
			nv1 = src.Get(fv1, DB_MAXARR);
			if (arrayMatches(fv1, nv1, &com1, fv0, nv0)) {
				for (uint32 i = 0; com1.Find(&i,false); i++) {
					if (nv0 >= DB_MAXARR) break;
					fv0[nv0++] = fv1[i];
				}
			} else {
				if (nv0 + nv1 > DB_MAXARR)
					nv1 = DB_MAXARR - nv0;
				memcpy(fv0+nv0, fv1, nv1*sizeof(float));
				nv0 += nv1;
			}
			dst.Set(fv0, nv0);
		}
		break;
	case DBDTstring:
		{
			DBField		dstdup(dst);
			const char 	*sv0[DB_MAXARR], *sv1[DB_MAXARR];
			nv0 = dstdup.Get(sv0, DB_MAXARR);
			nv1 = src.Get(sv1, DB_MAXARR);
			if (arrayMatches(sv1, nv1, &com1, sv0, nv0)) {
				for (uint32 i = 0; com1.Find(&i,false); i++) {
					if (nv0 >= DB_MAXARR) break;
					sv0[nv0++] = sv1[i];
				}
			} else {
				if (nv0 + nv1 > DB_MAXARR)
					nv1 = DB_MAXARR - nv0;
				memcpy(sv0+nv0, sv1, nv1*sizeof(const char *));
				nv0 += nv1;
			}
			dst.Set(sv0, nv0);
		}
		break;
	default:
		break;
	}
	return dst;
}

// Array field intersecton operator
DBField &
operator&=(DBField &dst, const DBField &src)
{
	if (src.GetNV() <= 0) {
		dst.Clear();
		return dst;
	}
	if (dst.GetNV() <= 0)
		return dst;
	if (src.GetDT() != dst.GetDT()) {
		DMESG(DMCparameter, "Mismatched types in DBField::operator&=()");
		return dst;
	}
	int		nv0, nv1;
	ABitMap		com0;
	switch (dst.GetDT()) {
	case DBDTint:
		{
			int32		iv0[DB_MAXARR], iv1[DB_MAXARR];
			nv0 = dst.Get(iv0, DB_MAXARR);
			nv1 = src.Get(iv1, DB_MAXARR);
			if (arrayMatches(iv0, nv0, &com0, iv1, nv1)) {
				nv1 = 0;
				for (uint32 i = 0; com0.Find(&i,true); i++)
					iv1[nv1++] = iv0[i];
				dst.Set(iv1, nv1);
			} else
				dst.Clear();
		}
		break;
	case DBDTfloat:
		{
			float		fv0[DB_MAXARR], fv1[DB_MAXARR];
			nv0 = dst.Get(fv0, DB_MAXARR);
			nv1 = src.Get(fv1, DB_MAXARR);
			if (arrayMatches(fv0, nv0, &com0, fv1, nv1)) {
				nv1 = 0;
				for (uint32 i = 0; com0.Find(&i,true); i++)
					fv1[nv1++] = fv0[i];
				dst.Set(fv1, nv1);
			} else
				dst.Clear();
		}
		break;
	case DBDTstring:
		{
			DBField		dstdup(dst);
			const char 	*sv0[DB_MAXARR], *sv1[DB_MAXARR];
			nv0 = dstdup.Get(sv0, DB_MAXARR);
			nv1 = src.Get(sv1, DB_MAXARR);
			if (arrayMatches(sv0, nv0, &com0, sv1, nv1)) {
				nv1 = 0;
				for (uint32 i = 0; com0.Find(&i,true); i++)
					sv1[nv1++] = sv0[i];
				dst.Set(sv1, nv1);
			} else
				dst.Clear();
		}
		break;
	default:
		break;
	}
	return dst;
}

// Array field difference operator
DBField &
operator-=(DBField &dst, const DBField &src)
{
	if (src.GetNV() <= 0)
		return dst;
	if (dst.GetNV() <= 0)
		return dst;
	if (src.GetDT() != dst.GetDT()) {
		DMESG(DMCparameter, "Mismatched types in DBField::operator-=()");
		return dst;
	}
	int		nv0, nv1;
	ABitMap		com0;
	switch (dst.GetDT()) {
	case DBDTint:
		{
			int32		iv0[DB_MAXARR], iv1[DB_MAXARR];
			nv0 = dst.Get(iv0, DB_MAXARR);
			nv1 = src.Get(iv1, DB_MAXARR);
			if (arrayMatches(iv0, nv0, &com0, iv1, nv1)) {
				nv1 = 0;
				for (uint32 i = 0; com0.Find(&i,false); i++)
					iv1[nv1++] = iv0[i];
				dst.Set(iv1, nv1);
			}
		}
		break;
	case DBDTfloat:
		{
			float		fv0[DB_MAXARR], fv1[DB_MAXARR];
			nv0 = dst.Get(fv0, DB_MAXARR);
			nv1 = src.Get(fv1, DB_MAXARR);
			if (arrayMatches(fv0, nv0, &com0, fv1, nv1)) {
				nv1 = 0;
				for (uint32 i = 0; com0.Find(&i,false); i++)
					fv1[nv1++] = fv0[i];
				dst.Set(fv1, nv1);
			}
		}
		break;
	case DBDTstring:
		{
			DBField		dstdup(dst);
			const char 	*sv0[DB_MAXARR], *sv1[DB_MAXARR];
			nv0 = dstdup.Get(sv0, DB_MAXARR);
			nv1 = src.Get(sv1, DB_MAXARR);
			if (arrayMatches(sv0, nv0, &com0, sv1, nv1)) {
				nv1 = 0;
				for (uint32 i = 0; com0.Find(&i,false); i++)
					sv1[nv1++] = sv0[i];
				dst.Set(sv1, nv1);
			}
		}
		break;
	default:
		break;
	}
	return dst;
}

#define LAB		'{'		// left array bracket
#define SEP		','		// array separator
#define RAB		'}'		// right array bracket
#define QUO		'"'		// string quote character
#define NQO		'\\'		// escape character

/*
 * Read a field in from a stream.
 * Stream must be started on field -- no leading whitespace
 * Type is determined by formatting:
 *	Strings must be in double quotes
 *	Floating point needs decimal point or exponent
 * Arrays values must separated by commas and no whitespace,
 * and may optionally be enclosed in curly braces {}
 */
istream &
operator>>(istream &ci, DBField &f)
{
	char		buf[30000];
	char *		argp[DB_MAXARR+1];
	int		nargs;
	DBDataType	dt;
	char *		cp;
	int		c;
	argp[nargs=0] = NULL;			// get field value(s)
	cp = buf;
	dt = DBDTunknown;
	while (ci.good()) {
		if (cp >= buf+sizeof(buf)-1) {
			argp[nargs] = NULL;	// buffer overflow!
			ci.clear(ios::failbit);
			break;
		}
		c = ci.get();			// get next character
		if (dt == DBDTunknown) {	// determine type
			if (c == QUO)
				dt = DBDTstring;
			else if (isdigit(c) || (c == '+') | (c == '-'))
				dt = DBDTint;	// could also be float
			else if (c == '.')
				dt = DBDTfloat;
			else if (c == LAB)
				continue;
		}
		switch (dt) {			// continue checking value
		case DBDTunknown:		// bogus value
			break;
		case DBDTstring:		// string value
			if (c == QUO) {		// beginning or end
				if (argp[nargs] == NULL) {
					argp[nargs] = cp;
					continue;
				}
				c = ci.get();	// unquote (end)
				break;
			}
			if (c == NQO) {		// literal next or special
				if (argp[nargs] == NULL)
					argp[nargs] = cp;
				c = ci.get();
				switch (c) {
				case EOF: break;
				case '\n': break;
				case '\r': break;
				case 'n': *cp++ = '\n'; break;
				case 'r': *cp++ = '\r'; break;
				case 't': *cp++ = '\t'; break;
				default: *cp++ = c; break;
				}
				continue;
			}
			if (c == EOF || argp[nargs] == NULL)
				break;		// missing quote
			*cp++ = c;		// continue in string
			continue;
		case DBDTint:			// integer value
			if ((c == '.') | (c == 'E') | (c == 'e'))
				dt = DBDTfloat;	// must be float after all
			else if (!isdigit(c) && (c != '+') & (c != '-'))
				break;		// end of integer
			if (argp[nargs] == NULL)
				argp[nargs] = cp;
			*cp++ = c;		// continue read
			continue;
		case DBDTfloat:			// float value
			if (!isdigit(c) && (c != '.') & (c != 'E') & (c != 'e') &
					(c != '+') & (c != '-'))
				break;		// end of float
			if (argp[nargs] == NULL)
				argp[nargs] = cp;
			*cp++ = c;		// continue in float
			continue;
		}
		if (argp[nargs] == NULL) {	// bad value?
			if (c != EOF)
				ci.putback(c);
			dt = DBDTunknown;
			break;
		}
		*cp++ = '\0';			// terminate value
		if (c == EOF)
			ci.clear(ios::eofbit);	// unset failbit
		argp[++nargs] = NULL;
		if ((nargs >= DB_MAXARR) | (c != SEP)) {
			if ((c != EOF) & (c != RAB))
				ci.putback(c);
			if (c == SEP)		// too many values
				ci.clear(ios::failbit);
			break;			// end of values
		}
	}
	switch (dt) {				// assign field
	case DBDTunknown:
		f.Clear();
		ci.clear(ios::failbit);
		break;
	case DBDTstring:
		f.Set(argp, nargs);
		break;
	case DBDTint:
		{
			int32	ia[DB_MAXARR];
			for (int i = nargs; i--; )
				ia[i] = (int32)atol(argp[i]);
			f.Set(ia, nargs);
		}
		break;
	case DBDTfloat:
		{
			float	fa[DB_MAXARR];
			for (int i = nargs; i--; )
				fa[i] = (float)atof(argp[i]);
			f.Set(fa, nargs);
		}
		break;
	}
	return ci;				// return stream
}

// Output quoted string
static void
quoteString(ostream *os, const char *str)
{
	os->put(QUO);
	while (*str) {
		switch (*str) {
		case '\n':
			*os << "\\n";
			break;
		case '\r':
			*os << "\\r";
			break;
		case '\t':
			*os << "\\t";
			break;
		case QUO:
		case NQO:
			os->put(NQO);
		    // fall through...
		default:
			os->put(*str);
			break;
		}
		str++;
	}
	os->put(QUO);
}

// Write a field out to a stream
// Strings are in quotes
// Arrays values are enclosed in curly braces {} and separated by commas
ostream &
operator<<(ostream &co, const DBField &f)
{
	if (f.GetNV() <= 0)
		return co;
	switch (f.GetDT()) {
	case DBDTint:
		{
			int32	iv[DB_MAXARR];
			int	nv = f.Get(iv, DB_MAXARR);
			if (nv > 1) {
				co << LAB << iv[0];
				for (int i = 1; i < nv; i++)
					co << SEP << iv[i];
				co << RAB;
			} else
				co << iv[0];
		}
		break;
	case DBDTfloat:
		{
			ios_base::fmtflags	oflg = co.flags();
			int			oprc = co.precision();
			float			fv[DB_MAXARR];
			int			nv = f.Get(fv, DB_MAXARR);
			co.setf(ios::scientific, ios::floatfield);
			co.precision(6);
			if (nv > 1) {
				co << LAB << fv[0];
				for (int i = 1; i < nv; i++)
					co << SEP << fv[i];
				co << RAB;
			} else
				co << fv[0];
			co.precision(oprc);
			co.flags(oflg);
		}
		break;
	case DBDTstring:
		{
			const char *	sv[DB_MAXARR];
			int		nv = f.Get(sv, DB_MAXARR);
			if (nv > 1) {
				co << LAB;
				quoteString(&co, sv[0]);
				for (int i = 1; i < nv; i++) {
					co.put(SEP); quoteString(&co, sv[i]);
				}
				co << RAB;
			} else
				quoteString(&co, sv[0]);
		}
		break;
	default:
		DMESG(DMCparameter, "Unknown data type in operator<<");
		break;
	}
	return co;
}

#undef SEP	
#undef QUO	
#undef NQO

// Grow linked field array to hold n items
void
DBRLink::GrowTo(int n)
{
	if (readonly)
		return;
	if (n <= nfields)
		return;
	if (nfields <= 0) {
		field = new DBField [nfields=n];
		return;
	}
					// the following happens rarely
	DBField	*	newfa = new DBField [n];
	for (int i = nfields; i--; )
		newfa[i] = field[i];
	delete [] field;
	field = newfa; nfields = n;
}

// Read unformatted record line from stream
int
DBRecord::Read(istream *ins, const char tabch)
{
	if (ReadOnly()) {
		DMESG(DMCparameter, "Attempt to read into read-only record");
		return 0;
	}
	if ((ins == NULL) | (finfo == NULL))
		return 0;
	int	nfset = 0;
	int	delim = tabch;
	l->Init(finfo->GetNFields());		// clear previous values
	ins->clear();				// read each field
	for (int i = 0; (i < l->nfields) & (delim == tabch); i++) {
		delim = ins->get();		// check next character
		if ((delim == EOF) | (delim == '\n') | (delim == '\r') |
				(delim == '}'))
			break;			// end of input
		if (delim == tabch)
			continue;		// empty field
		ins->putback(delim);		// must be field char
		*ins >> l->field[i];		// get field
		if (ins->fail())
			goto failure;
		delim = ins->get();		// get delimiter
		if (l->field[i].GetNV() <= 0)
			continue;		// empty field?
						// check type
		const DBFDetails *	fd = finfo->GetDetails(i);
		if (fd->dtype != DBDTunknown) {
			if (l->field[i].GetDT() != fd->dtype)
				goto failure;
			if (!(fd->flags & DBFFarray) && l->field[i].GetNV() > 1)
				goto failure;
		} else if (l->field[i].GetDT() == DBDTunknown)
			goto failure;
		++nfset;			// field seems OK
	}
	if (delim == tabch)
		do				// skip extra fields to endl
			delim = ins->get();
		while ((delim != EOF) & (delim != '\n') & (delim != '\r'));
	else if ((delim != EOF) & (delim != '\n') & (delim != '\r'))
		ins->putback(delim);		// replace mystery character
	if (!nfset && (delim != '\n') & (delim != '\r'))
		return -1;			// reading stymied/EOF
	return nfset;
failure:
	Clear();
	return -1;
}

// Compute hash on index fields (DBHVnone if incomplete)
long
DBRecord::ComputeHashIndex() const
{
	static const long	DB_HASHMASK = ((unsigned long)1<<DB_NHASHBITS)-1L;
	int	i;
	for (i = l->nfields; i--; )
		if (finfo->GetDetails(i)->flags & DBFFindex &&
				l->field[i].GetNV() <= 0)
			return l->ihash = DBHVnone;
	l->ihash = 0;			// build hash value
	for (i = l->nfields; i--; )
		if (finfo->GetDetails(i)->flags & DBFFindex) {
			const long	fhv = l->field[i].Hash(DB_NHASHBITS);
			const char *	fnm = finfo->GetName(i);
			const int	shift = (int)(strhash(fnm) % DB_NHASHBITS);
			l->ihash ^= (fhv<<shift & DB_HASHMASK) |
					fhv>>(DB_NHASHBITS-shift);
		}
	return l->ihash;
}

// Compare against another record using our sort info
// Uninitialized records and empty fields are "smaller"
int
DBRecord::Compare(const DBRecord &that) const
{
	if (that.finfo == NULL)
		return (finfo != NULL);
	if (finfo == NULL)
		return -1;
	if (sortord == NULL)
		return 0;
	int	i, cmp;
	bool	rev;
						// check each field in sort list
	for (i = 0; (cmp = sortord->GetSort(i, &rev)) >= 0; i++) {
		const DBField *	f1 = GetField(cmp);
		const DBField *	f2 = that.GetField(
					that.finfo->XlateID(finfo, cmp) );
		if ((f1 == NULL) & (f2 == NULL))
			continue;
		if (f1 == NULL)
			cmp = -1;
		else if (f2 == NULL)
			cmp = 1;
		else
			cmp = f1->Compare(*f2);
		if (cmp)			// found difference?
			return rev ? -cmp : cmp;
	}
	return 0;				// sort fields compared equal
}

// See if this record matches a database query selector list
bool
DBRecord::MatchQuery(const DBQuery *dq) const
{
	if ((dq == NULL) | (finfo == NULL))
		return false;
						// accept any selector match
	const DBFieldInfo *	mappedinf = NULL;
	int			idmap[DB_MAXFIELD];
	for (dq = dq->FirstPossible(HashIndex()); dq != NULL;
				dq = dq->NextPossible(l->ihash)) {
						// compare assigned fields
		if (dq->rlo.GetFieldInfo() != mappedinf)
			finfo->MapIDs(idmap, mappedinf = dq->rlo.GetFieldInfo());
		int	i = dq->rlo.GetNAlloc();
		if (i <= 0 || dq->rlo.GetNAssigned() <= 0)
			continue;
		while (i--) {
			const DBField *	f1 = dq->rlo.GetField(i);
			if (f1 == NULL) continue;
			const DBField *	f2 = dq->rhi.GetField(i);
			if (f2 == NULL) f2 = f1;
			int	j = idmap[i];
			if ((j < 0) | (j >= l->nfields))
				break;
			if (l->field[j].GetNV() <= 0) {
				if (dq->matchEmpty) continue;
				break;
			}
			if (!l->field[j].MatchRange(*f1, *f2))
				break;
		}
		if (i < 0)
			return true;		// got one that matched!
	}
	return false;
}

// Check for exact match of assigned record fields
bool
DBRecord::operator==(const DBRecord &that) const
{
	if (this == &that)
		return true;
	if (l == that.l)
		return true;
	if (GetNAlloc() <= 0)
		return (that.GetNAlloc() <= 0);
	if (that.GetNAlloc() <= 0)
		return false;
	if (HashIndex() != that.HashIndex())
		return false;
	int	nmatched = 0;
	int	idmap[DB_MAXFIELD];		// map field IDs
	finfo->MapIDs(idmap, that.finfo);
	int	i = l->nfields;
	if (i > that.l->nfields)
		i = that.l->nfields;
	else
		while (i < that.l->nfields)
			idmap[i++] = -1;
	while (i--)
		if (that.l->field[i].GetNV() > 0) {
			if (idmap[i] < 0 || l->field[idmap[i]] != that.l->field[i])
				return false;
			++nmatched;
		}
	return (nmatched == GetNAssigned());
}

// Record comparison function for qsort()
int
DBrecordCmp(const void *rp1, const void *rp2)
{
	return (*(const DBRecord *)rp1).Compare(*(const DBRecord *)rp2);
}

// Copy database record holder, translating fields if we're locked
// Final lock state is same as before
DBRecord &
DBRecord::operator=(const DBRecord &src)
{
	if (this == &src)
		return *this;
	Init(src.finfo, false);
	int	i = src.GetNAlloc();
	if (i <= 0)
		return *this;
	int		idmap[DB_MAXFIELD];
	finfo->MapIDs(idmap, src.finfo);	// map field IDs
	l->Init(finfo->GetNFields());
	while (i--) {
		if (src.l->field[i].GetNV() <= 0)
			continue;		// no value to transfer
		if (finfo != src.finfo) {	// check type info
			const DBFDetails *	fd = GetFieldDetails(idmap[i]);
			if (fd == NULL)
				continue;	// unknown field
			if (fd->dtype != DBDTunknown &&
					fd->dtype != src.l->field[i].GetDT())
				continue;	// type mismatch
			if (!(fd->flags & DBFFarray) &&
					src.l->field[i].GetNV() > 1) {
				src.l->field[i].Get(&l->field[idmap[i]]);
				continue;	// copied first entry
			}
		}
						// exact field copy
		l->field[idmap[i]] = src.l->field[i];
	}
	return *this;
}

// Append new values from source record (union may add to arrays)
DBRecord &
DBRecord::operator|=(const DBRecord &src)
{
	if (src.GetNAlloc() <= 0)
		return *this;
	if (GetNAlloc() <= 0)
		return *this = src;
	if (ReadOnly())
		return *this;
	int		idmap[DB_MAXFIELD];	// map our fields to theirs
	int		i;
	src.finfo->MapIDs(idmap, finfo);
	l->GrowTo(finfo->GetNFields());
	for (i = l->nfields; i--; ) {
		if (idmap[i] < 0)
			continue;
		const DBField &		sfld = src.l->field[idmap[i]];
		const DBFDetails *	fd = finfo->GetDetails(i);
		if (sfld.GetDT() != fd->dtype)
			continue;
		if (l->field[i].GetNV() <= 0) {
			SetField(i, sfld);
			continue;
		}
		if (sfld.GetNV() > 0 && fd->flags & DBFFarray) {
			l->field[i] |= sfld;
			if (fd->flags & DBFFindex)
				l->ihash = DBHVunknown;
			continue;
		}
	}
	return *this;
}

// Leave values that match source record assignments (intersection)
DBRecord &
DBRecord::operator&=(const DBRecord &src)
{
	if (src.GetNAlloc() <= 0)
		return *this;
	if (GetNAlloc() <= 0)			// first assignment exception
		return *this = src;
	if (ReadOnly())
		return *this;
	int		idmap[DB_MAXFIELD];	// map our fields to theirs
	int		i;
	src.finfo->MapIDs(idmap, finfo);
	l->GrowTo(finfo->GetNFields());
	for (i = l->nfields; i--; ) {
		if (idmap[i] < 0)
			continue;
		if (l->field[i].GetNV() <= 0)
			continue;
		const DBField &		sfld = src.l->field[idmap[i]];
		const DBFDetails *	fd = finfo->GetDetails(i);
		if (sfld.GetDT() != fd->dtype)
			continue;
		if (fd->flags & DBFFarray) {
			l->field[i] &= sfld;
			if (fd->flags & DBFFindex)
				l->ihash = DBHVunknown;
			continue;
		}
		if (sfld.GetNV() <= 0 || l->field[i] != sfld) {
			ClearField(i);
			continue;
		}
	}
	return *this;
}

// Fill unassigned fields from another record
int
DBRecord::FillFrom(const DBRecord &rec)
{
	if (ReadOnly())
		return -1;
	int		nfilled = 0;
	int		idmap[DB_MAXFIELD];	// map our fields to theirs
	int		i;
	rec.finfo->MapIDs(idmap, finfo);
	l->GrowTo(finfo->GetNFields());
	for (i = l->nfields; i--; )
		if (l->field[i].GetNV() <= 0 && idmap[i] >= 0)
			nfilled += SetField(i, rec.l->field[idmap[i]]);
	return nfilled;
}

// Prepare a record field for assignment (internal)
bool
DBRecord::PrepField(int fn, DBDataType dt, int nv)
{
	if (ReadOnly())
		return false;
	const DBFDetails *	fd = GetFieldDetails(fn);
	if (fd == NULL)
		return false;
	if (fn < l->nfields)
		l->field[fn].Clear();
	if (fd->flags & DBFFindex)
		l->ihash = DBHVunknown;
	if (dt == DBDTunknown)
		return false;
	if ((fd->dtype != DBDTunknown) & (fd->dtype != dt))
		return false;
	if (nv > 1 && !(fd->flags & DBFFarray))
		return false;
	if (fn >= l->nfields)
		l->GrowTo(finfo->GetNFields());
	return (nv > 0);
}

// Set an array of integer values in a record field
bool
DBRecord::SetField(int fn, const int32 ia[], int n)
{
	if ((n < 1) | (ia == NULL)) {
		ClearField(fn);
		return false;
	}
	if (!PrepField(fn, DBDTint, n))
		return false;
	l->field[fn].Set(ia, n);
	return true;
}

// Set an array of float values in a record field
bool
DBRecord::SetField(int fn, const float fa[], int n)
{
	if ((n < 1) | (fa == NULL)) {
		ClearField(fn);
		return false;
	}
	if (!PrepField(fn, DBDTfloat, n))
		return false;
	l->field[fn].Set(fa, n);
	return true;
}

// Set an array of string values in a record field
bool
DBRecord::SetField(int fn, const char * const sa[], int n)
{
	if ((n < 1) | (sa == NULL)) {
		ClearField(fn);
		return false;
	}
	if (!PrepField(fn, DBDTstring, n))
		return false;
	l->field[fn].Set(sa, n);
	return true;
}

// Set a string array from text with lines separated by the given character
int
DBRecord::SetText(int fn, const char *txt, const char nl)
{
	if (txt == NULL || !*txt) {
		ClearField(fn);
		return 0;
	}
	if (!PrepField(fn, DBDTstring, 2))
		return 0;
	return l->field[fn].SetText(txt, nl);
}

// Format field value for display
char *
DBRecord::FormatField(int fn, char *buf, int len) const
{
	const DBField *		f = GetField(fn);
	if (f == NULL)
		return NULL;
	const DBFDetails *	fd = GetFieldDetails(fn);
	if (fd == NULL)
		return NULL;
	if (fd->DisplayFormat == NULL) {
		if (f->GetNV() > 1)
			return DBformatArray(buf, len, *f);
		return DBformatUnknown(buf, len, *f);
	}
	return (*fd->DisplayFormat)(buf, len, *f);
}

// Resize record array
int
DBRecordList::Resize(int n)
{
	int	i;
	if (n < 0)
		n = 0;
	if (n == nr)			// check special cases
		return nr;
	if ((maxr <= 0) | (rl == NULL)) {
		if (!n)
			return 0;
		maxr = n;		// initial allocation
		rl = (DBRecord *)malloc(maxr*sizeof(DBRecord));
		if (rl == NULL) {
			DMESG(DMCmemory, "malloc() failed in Resize");
			return nr = maxr = 0;
		}
		for (i = maxr; i--; )
			rl[i].Empty();
		if (finfo != NULL)	// initialize new records
			for (i = n; i--; )
				rl[i].Init(finfo, true);
		return nr = n;		// return fresh array
	}
					// resize: free unwanted records
	while (nr > n)
		rl[--nr].Init();

	if (!nr) {			// freed all?
		free(rl);
		rl = NULL;
		return nr = maxr = 0;
	}
					// reallocate array if necessary
	if (n > maxr || maxr > 2*n + 32) {
		const int	nold = maxr;
		maxr = n + 32;		// give a little room to grow
		rl = (DBRecord *)realloc(rl, maxr*sizeof(DBRecord));
		if (rl == NULL) {
			DMESG(DMCmemory, "realloc() failed in Resize");
			return nr = maxr = 0;
		}
		for (i = maxr; i-- > nold; )
			rl[i].Empty();
	}
	if (finfo != NULL)		// initialize new records
		for (i = n; i-- > nr; )
			rl[i].Init(finfo, true);
	return nr = n;			// all set
}

// Eliminate empty records from list
void
DBRecordList::Compact(bool delEmpty)
{
	if ((nr <= 0) | (rl == NULL))
		return;
	DBRecord * const	term = rl + nr;
	DBRecord *		dp = rl;
	DBRecord *		sp;
	for (sp = rl; sp < term; sp++) {
		if (sp->GetNFields() <= 0)
			continue;
		if (delEmpty && !sp->GetNAssigned()) {
			sp->Init();
			continue;
		}
		if (dp < sp)
			dp->Take(sp);
		dp++;
	}
	Resize(dp - rl);		// fix array size
}

// Sort record array
void
DBRecordList::Sort()
{
	if ((nr <= 0) | (rl == NULL))
		return;
	for (int i = nr*(sord.GetSort(0) >= 0); i--; )
		rl[i].sortord = &sord;
	qsort(rl, nr, sizeof(DBRecord), DBrecordCmp);
}

// Link all records from another list
int
DBRecordList::LinkAll(DBRecordList *srcp)
{
	if (this == srcp)
		return nr;
	if (srcp == NULL) {
		Init();
		return 0;
	}
	if (!Init(NULL, srcp->nr))
		return 0;
	finfo = srcp->finfo;
	for (int i = nr; i--; ) {
		if (!rl[i].Link(&srcp->rl[i])) {
			Init();
			return 0;
		}
		rl[i].fieldlock = (finfo != NULL);
	}
	if (finfo != NULL)
		sord = finfo->sord;
	return nr;
}

// Record array copy assignment
// If unlocked, we take on lock of src
DBRecordList &
DBRecordList::operator=(const DBRecordList &src)
{
	if (this == &src)
		return *this;
	if (finfo == NULL)
		Init(src.finfo, src.nr);
	else
		Resize(src.nr);
	for (int i = nr; i--; )
		rl[i] = src.Get(i);
	return *this;
}

// Read records from stream, adding to our list.
// If set, fi contains the field names, else they will be read from a header line. 
// If our list has no standard field list, then fi must outlive result.
int
DBRecordList::Read(istream *ins, char *tabcp, DBFieldInfo *fi)
{
	if ((finfo == NULL) & (fi == NULL))
		return -1;			// no header storage!
	if (ins == NULL)
		return 0;
	char		tabch = '\0';		// default delimiter (any)
	DBFieldInfo	tinf;
	int		nadd;
	if (tabcp == NULL)			// need tab holder
		tabcp = &tabch;
	if (fi == NULL)				// need header holder
		fi = &tinf;
						// read header if none given
	if (!fi->GetNFields() && (nadd = fi->Read(ins, tabcp)) <= 0)
		return nadd;
	if (!*tabcp)				// need delimiter
		*tabcp = '\t';
	DBRecord	rec(fi);		// init record holder
	rec.fieldlock = true;			// not really necessary
	nadd = 0;
	while (ins->good()) {			// read & add each record
		int	nf = rec.Read(ins, *tabcp);
		if (nf < 0)
			break;			// EOF or read error
		if (nf == 0)
			continue;		// empty record
		NewRecord() = rec;		// add record to list
		++nadd;
	}
	return nadd;
}

// Assign hash index set for the entire query list
bool
DBQuery::HashSet(ABitMap *hs)
{
	DBQuery *	qp;
	if (hs == NULL || !hs->Length())
		goto reset;
	if (next == NULL) {
		hs->NewBitMap(0);		// don't bother
		hashSet = NULL;
		return false;
	}
	hs->ClearBitMap(false);
	for (qp = this; qp != NULL; qp = qp->next) {
		long	ih = qp->HashIndex();
		if (ih < 0) {
			hs->NewBitMap(0);	// missing index
			goto reset;
		}
		hs->Set(ih % hs->Length());
		qp->hashSet = hs;
	}
	return true;
reset:						// no index
	for (qp = this; qp != NULL; qp = qp->next)
		qp->hashSet = NULL;
	return false;
}

// Return first possible query to match the given index hash value
const DBQuery *
DBQuery::FirstPossible(long ih) const
{
	if (this == NULL)
		return NULL;
	if ((ih < 0) & matchEmpty)
		return this;
	if (HashReject(ih))
		return NULL;
	const long	th = HashIndex();
	if ((th < 0) | (th == ih))
		return this;
	return NextPossible(ih);
}

// Add new selector to search criteria
bool
DBSearch::AddSelector(const DBRecord &r1, const DBRecord &r2, bool me)
{
	int	n = r1.GetNAlloc();
	if (n <= 0)
		return false;
	int	idmap[DB_MAXFIELD];
	if (dbacc == NULL || dbacc->GetHeader() == NULL ||
			r1.GetFieldInfo()->MapIDs(idmap, r2.GetFieldInfo()))
		return false;
	crpos = -1;				// reset search & create new selector
	if (qlist == NULL) {
		qlist = new DBQuery(dbacc->GetHeader());
	} else if (qlist->rlo.GetNAssigned()) {
		qlist = new DBQuery(qlist);
		if (ihSet.Length() != DB_HASHSETLEN)
			ihSet.NewBitMap(DB_HASHSETLEN);
	}
						// set empty match boolean
	qlist->matchEmpty = me;
						// map field IDs
	qlist->rlo.GetFieldInfo()->MapIDs(idmap, r1.GetFieldInfo());
	bool	gotfield = false;
	while (n--) {
		if (idmap[n] < 0)
			continue;		// XXX we should issue warning?
		const DBField *	f1 = r1.GetField(n);
		const DBField *	f2 = r2.GetField(n);
		if (f1 == NULL)
			continue;
		if (f2 == NULL)
			f2 = f1;
		if (f1->GetDT() != f2->GetDT()) {
			DMESGF(DMCwarning, "Type mismatch for query field %s",
					r1.GetFieldName(n));
			continue;
		}
		if (f1->GetNV() != f2->GetNV()) {
			DMESGF(DMCwarning, "Pair mismatch in query field %s",
					r1.GetFieldName(n));
			continue;
		}
		int	i;
		switch (f1->GetDT()) {
		case DBDTint: {			// transfer integer range(s)
			int32	ia1[DB_MAXARR], ia2[DB_MAXARR];
			int	ni = f1->Get(ia1, DB_MAXARR);
			f2->Get(ia2, DB_MAXARR);
			for (i = ni; i--; )
				if (ia1[i] > ia2[i]) {
					int32	t = ia1[i];
					ia1[i] = ia2[i]; ia2[i] = t;
				}
			if (!qlist->rlo.SetField(idmap[n], ia1, ni) ||
					!qlist->rhi.SetField(idmap[n], ia2, ni)) {
				DMESGF(DMCwarning, "Cannot assign query field %s",
						r1.GetFieldName(n));
				qlist->rlo.ClearField(idmap[n]);
				qlist->rhi.ClearField(idmap[n]);
				continue;
			}
			} break;
		case DBDTfloat: {		// transfer real range(s)
			float	fa1[DB_MAXARR], fa2[DB_MAXARR];
			int	nf = f1->Get(fa1, DB_MAXARR);
			f2->Get(fa2, DB_MAXARR);
			for (i = nf; i--; ) {
				if (fa1[i] > fa2[i]) {
					float	t = fa1[i];
					fa1[i] = fa2[i]; fa2[i] = t;
				}
			}
			if (!qlist->rlo.SetField(idmap[n], fa1, nf) ||
					!qlist->rhi.SetField(idmap[n], fa2, nf)) {
				DMESGF(DMCwarning, "Cannot assign query field %s",
						r1.GetFieldName(n));
				qlist->rlo.ClearField(idmap[n]);
				qlist->rhi.ClearField(idmap[n]);
				continue;
			}
			} break;
		case DBDTstring: {		// transfer string range(s)
			const char	*sa1[DB_MAXARR], *sa2[DB_MAXARR];
			int		ns = f1->Get(sa1, DB_MAXARR);
			f2->Get(sa2, DB_MAXARR);
			for (i = ns; i--; )
				if (strcmp_nc(sa1[i], sa2[i]) > 0) {
					const char *	t = sa1[i];
					sa1[i] = sa2[i]; sa2[i] = t;
				}
			if (!qlist->rlo.SetField(idmap[n], sa1, ns) ||
					!qlist->rhi.SetField(idmap[n], sa2, ns)) {
				DMESGF(DMCwarning, "Cannot assign query field %s",
						r1.GetFieldName(n));
				qlist->rlo.ClearField(idmap[n]);
				qlist->rhi.ClearField(idmap[n]);
				continue;
			}
			} break;
		default:			// botched type!
			DMESG(DMCparameter, "Missing type in AddSelector");
			return false;
		}
		gotfield = true;
	}
	return gotfield;
}

// Find all records matching search criteria
int
DBSearch::FindAll(DBRecordList *rl)
{
	if (rl == NULL)
		return 0;
	if (rl->GetSize() <= 0 && !rl->Resize(64))
		return 0;
	Reset();
						// find all we can
	int	nfound = 0;
	int	n;
	while ((n = FindRecords(rl->Array()+nfound, rl->GetSize()-nfound)) > 0)
		if ((nfound += n) >= rl->GetSize())
			rl->Resize(rl->GetSize()+128);
						// match final list length
	return rl->Resize(nfound);
}

// Find and delete matching records from database
int
DBSearch::DeleteMatching()
{
	int		nfnd = 0;
	int		ndel = 0;
	DBRecord	rlist[512];
	int		nr;
	
	Reset();
	while ((nr = FindRecords(rlist, 512)) > 0) {
		nfnd += nr;
		nr = dbacc->DeleteRecords(rlist, nr);
		ndel += nr;
		if ((crpos -= nr) < 0)
			crpos = 0;
	}
	if (nfnd > ndel)
		DMESGF(DMCwarning, "%d records missing in DeleteMatching", nfnd-ndel);
	return ndel;
}
