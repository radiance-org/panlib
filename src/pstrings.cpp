/*
 *  pstrings.cpp
 *  panlib
 *
 *  Pancine strings, string sets, and databaase header info.
 *
 *  Created by Greg Ward on Wed Apr 14 2004.
 *  Copyright (c) 2004 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "pstrings.h"
#include "dmessage.h"

#if defined(_WIN32) || defined(_WIN64)
#define atoll	_atoi64			// substitute under Windows
#endif

#ifndef  S_NHASH
#define  S_NHASH	131071		// string hash table size (prime)
#endif

// Caseless string compare (also special treatment of numbers)
int
strcmp_nc(const char *so1, const char *so2)
{
	if (so1 == so2) return 0;
	if (so1 == NULL) return -1;
	if (so2 == NULL) return 1;
	const char *	s1 = so1;
	const char *	s2 = so2;
	int		c;
	long long	idif;
more:
	for ( ; !(c = tolower(*s1) - tolower(*s2)); s1++)
		if (!*s2++)
			return 0;
					// special number check
	if (isdigit(*s1) | isdigit(*s2)) {
		if (!isdigit(*s2))
			return 1;
		if (!isdigit(*s1))
			return -1;
		while (s1 > so1 && isdigit(s1[-1])) {
			--s1; --s2;	// we know they're the same
		}
		if (!(idif = atoll(s1) - atoll(s2))) {
			do ++s1; while (isdigit(*s1));
			do ++s2; while (isdigit(*s2));
			goto more;
		}
		c = (idif > 0) ? 1 : -1;
	}
	return c;
}

// Caseless find string in string
const char *
strstr_nc(const char *big, const char *little)
{
	if (little == NULL) return NULL;
	if (big == NULL) return NULL;
	if (little == big) return big;
	if (!*little) return big;
	const char	fc = tolower(*little);
	for ( ; *big; ++big)
		if (tolower(*big) == fc) {
			const char *	lp = little;
			const char *	bp = big;
			do {
				if (!*++lp)
					return big;
				++bp;
			} while (tolower(*bp) == tolower(*lp));
		}
	return NULL;
}

// Compute nbits-wide hash value from memory (do NOT alter this code)
unsigned long
memhash(const void *p, int len, int nbits)
{
	static const unsigned char	shuffle[256] = {
		0, 157, 58, 215, 116, 17, 174, 75, 232, 133, 34,
		191, 92, 249, 150, 51, 208, 109, 10, 167, 68, 225,
		126, 27, 184, 85, 242, 143, 44, 201, 102, 3, 160,
		61, 218, 119, 20, 177, 78, 235, 136, 37, 194, 95,
		252, 153, 54, 211, 112, 13, 170, 71, 228, 129, 30,
		187, 88, 245, 146, 47, 204, 105, 6, 163, 64, 221,
		122, 23, 180, 81, 238, 139, 40, 197, 98, 255, 156,
		57, 214, 115, 16, 173, 74, 231, 132, 33, 190, 91,
		248, 149, 50, 207, 108, 9, 166, 67, 224, 125, 26,
		183, 84, 241, 142, 43, 200, 101, 2, 159, 60, 217,
		118, 19, 176, 77, 234, 135, 36, 193, 94, 251, 152,
		53, 210, 111, 12, 169, 70, 227, 128, 29, 186, 87,
		244, 145, 46, 203, 104, 5, 162, 63, 220, 121, 22,
		179, 80, 237, 138, 39, 196, 97, 254, 155, 56, 213,
		114, 15, 172, 73, 230, 131, 32, 189, 90, 247, 148,
		49, 206, 107, 8, 165, 66, 223, 124, 25, 182, 83,
		240, 141, 42, 199, 100, 1, 158, 59, 216, 117, 18,
		175, 76, 233, 134, 35, 192, 93, 250, 151, 52, 209,
		110, 11, 168, 69, 226, 127, 28, 185, 86, 243, 144,
		45, 202, 103, 4, 161, 62, 219, 120, 21, 178, 79,
		236, 137, 38, 195, 96, 253, 154, 55, 212, 113, 14,
		171, 72, 229, 130, 31, 188, 89, 246, 147, 48, 205,
		106, 7, 164, 65, 222, 123, 24, 181, 82, 239, 140,
		41, 198, 99
	};
	static const int		MAXBITS = int(8*sizeof(unsigned long));
	static const short		shiftinc[65] = {
		1, 1, 1, 1, 1, 1, 1, 1, 1,
		1, 3, 3, 5, 5, 5, 7, 7,
		7, 7, 7, 7, 11, 13, 11, 11,
		11, 11, 11, 11, 11, 11, 11, 11,
		13, 13, 13, 13, 13, 13, 17, 17,
		17, 17, 17, 17, 17, 17, 17, 17,
		17, 17, 19, 19, 19, 19, 19, 19,
		23, 23, 23, 23, 23, 23, 23, 23
	};
	const unsigned char *		cp = (const unsigned char *)p;
	unsigned long			hval = 0L;
	int				shift = 0;
	DASSERT(MAXBITS <= 64);
	if ((nbits <= 0) | (nbits > MAXBITS))
		nbits = MAXBITS;
	const int			myShiftInc = shiftinc[nbits];
	while (len-- > 0) {
		unsigned long	v = (unsigned long)shuffle[*cp++];
		if (shift+8 > nbits)
			hval ^= v >> (nbits-shift);
		hval ^= v << shift;
		if ((shift += myShiftInc) >= nbits)
			shift -= nbits;
	}
	if (nbits < MAXBITS)
		hval &= ((unsigned long)1<<nbits)-1L;
	return hval;
}

// Compute nbits-wide hash value from string (caseless)
unsigned long
strhash_nc(const char *s, int nbits)
{
	char		buf[4096];
	int		n;
	unsigned long	hval = 0L;

	while (*s) {
		for (n = 0; *s && n < sizeof(buf); s++)
			buf[n++] = tolower(*s);

		hval ^= memhash(buf, n, nbits);
	}
	return hval;
}

struct SHead {
	SHead *		next;		// next in hash list
	int		nl;		// links count
};				// followed by the string itself

static SHead *		stab[S_NHASH];

const char		strNoMemory[] = "OutOfMemory";

#define	S_HASH(s)	int(strhash(s) % S_NHASH)
#define	S_STR(sp)	((char *)((sp)+1))

// Share a common copy of a read-only string
const char *
strlink(const char *str)
{
	if (str == NULL)
		return NULL;
	if (!*str)
		return "";
	const int	hval = S_HASH(str);
	SHead *		sp;
	for (sp = stab[hval]; sp != NULL; sp = sp->next)
		if (!istrcmp(str, S_STR(sp))) {
			sp->nl++;
			return(S_STR(sp));
		}
	sp = (SHead *)malloc(sizeof(SHead)+1+strlen(str));
	if (sp == NULL) {
		DMESG(DMCmemory, "malloc() failed in strlink");
		return strNoMemory;
	}
	strcpy(S_STR(sp), str);
	sp->nl = 1;
	sp->next = stab[hval];
	stab[hval] = sp;
	return S_STR(sp);
}

// Release link to read-only string
void
strunlink(const char *s)
{
	if ((s == NULL) | (s == strNoMemory) || !*s)
		return;
	const int	hval = S_HASH(s);
	SHead		*spl, *sp;
	for (spl = NULL, sp = stab[hval]; sp != NULL; spl = sp, sp = sp->next)
		if (s == S_STR(sp)) {
			if (--sp->nl > 0)
				return;
			if (spl != NULL)
				spl->next = sp->next;
			else
				stab[hval] = sp->next;
			free(sp);
			return;
		}
	DMESG(DMCwarning, "Bogus call to strunlink()");
}

// Print our string hash table
void
strprinttab(FILE *fp)
{
	for (int hval = 0; hval < S_NHASH; hval++) {
		if (stab[hval] == NULL)
			continue;
		fprintf(fp, "%d:", hval);
		for (SHead *sp = stab[hval]; sp != NULL; sp = sp->next)
			fprintf(fp, " %d(%s)", sp->nl, S_STR(sp));
		fputc('\n', fp);
	}
}

// Find a string in our set
int
PStringSet::Find(const char *s) const
{
	if (s == NULL)
		return -1;
	int	i;
	if (ns < 12) {				// linear search
		for (i = ns; i-- > 0; )
			if (!istrcmp(s, sl[i]))
				return i;
		return -1;
	}
	int     ilower = 0, iupper = ns;
	int	c = iupper;			// binary search
	while ((i = (iupper + ilower) >> 1) != c) {
		c = istrcmp(s, sl[i]);
		if (c > 0)
			ilower = i;
		else if (c < 0)
			iupper = i;
		else
			return i;
		c = i;
	}
	return -1;
}

// Insert a string in set -- return true if added
bool
PStringSet::Insert(const char *s)
{
	if (s == NULL || InSet(s))
		return false;
	if (ns >= maxs-1) {
		const char **   newsl;
		maxs += (maxs>>3) + 16;		// give us some space
		if (sl == NULL)
			newsl = (const char **)malloc(maxs*sizeof(const char *));
		else
			newsl = (const char **)realloc(sl,
						maxs*sizeof(const char *));
		if (newsl == NULL) {
			DMESG(DMCmemory, "No more space in PStringSet::Insert");
			return false;
		}
		sl = newsl;
	}
	s = strlink(s);				// insert in our list
	int	i;
	for (i = ns; i > 0; i--) {
		if (istrcmp(sl[i-1], s) < 0)
			break;
		sl[i] = sl[i-1];
	}
	sl[i] = s;
	sl[++ns] = NULL;			// add NULL terminator
	return true;
}

// Delete a string from set -- return true if deleted
bool
PStringSet::Delete(const char *s)
{
	int	i = Find(s);
	if (i < 0)
		return false;
	strunlink(sl[i]);			// delete from our list
	for ( ; i < ns; i++)
		sl[i] = sl[i+1];
	--ns;					// copied terminator
	return true;
}

// Copy operator
PStringSet &
PStringSet::operator=(const PStringSet &src)
{
	if (this == &src)
		return *this;
	Clear();
	if (!src.ns)
		return *this;
	sl = (const char **)malloc((src.ns+1)*sizeof(const char *));
	if (sl == NULL) {
		DMESG(DMCmemory, "malloc() failed in PStringSet::operator=()");
		return *this;
	}
	ns = maxs = src.ns;
	for (int i = ns; i--; )
		sl[i] = strlink(src.sl[i]);
	sl[ns] = NULL;				// add NULL terminator
	return *this;
}

// Find near match in string set (ignores case)
const char *
typo_match(const char *entry, const PStringSet &sset)
{
	if (entry == NULL || sset.GetSize() <= 0)
		return NULL;
	if (sset.Find(entry) >= 0)
		return NULL;				// has perfect match
	const int	slen = strlen(entry);
	int		nearest = slen/10 + 2;		// need better than this
	const char *	bestMatch = NULL;
	int		n = sset.GetSize();
	while (n--) {					// find closest match
		const char *	sstr = sset[n];
		if (abs((int)strlen(sstr) - slen) >= nearest)
			continue;
		int		nmiss = 0;
		const char *	estr = entry;
		while (*estr)				// count typos
			if (tolower(*estr++) != tolower(*sstr++)) {
				if (!sstr[-1]) {
					--sstr;
					nmiss += strlen(--estr);
					break;
				}
							// sync up
				if (tolower(*estr) != tolower(*sstr)) {
					bool	eback = (tolower(estr[-1]) == tolower(*sstr));
					bool	sback = (tolower(sstr[-1]) == tolower(*estr));
					if (eback ^ sback) {
						estr -= eback;
						sstr -= sback;
					} else {	// transposed letters count as 1 typo
						estr += eback;
						sstr += sback;
					}
				}
				if (++nmiss >= nearest)
					break;
			}
		nmiss += strlen(sstr);
		if (nmiss < nearest) {
			bestMatch = sset[n];
			nearest = nmiss;
		}
	}
	if ((nearest > 0) & (slen <= 3))
		return NULL;				// too short for typos
	return bestMatch;				// may be case mismatch
}
