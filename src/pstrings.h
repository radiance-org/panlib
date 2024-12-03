/*
 *  pstrings.h
 *  panlib
 *
 *  Include after <stdio.h>, <string.h> and <stdlib.h>
 *
 *  Declarations for strings and string sets
 *
 *  Created by Greg Ward on Wed Apr 14 2004.
 *  Copyright (c) 2004 Anyhere Software. All rights reserved.
 *
 */

#ifndef _PSTRINGS_H_
#define _PSTRINGS_H_

// Inline string comparison function with pointer equality test
inline int
istrcmp(const char *s1, const char *s2)
{
	if (s1 == s2) return 0;
	if (!s1) return -1;
	if (!s2) return 1;
	while (*s1++ == *s2)
		if (!*s2++) return 0;
	return int(*--s1) - int(*s2);
}

// Caseless string compare
extern int		strcmp_nc(const char *s1, const char *s2);

// Caseless find string in string
extern const char *	strstr_nc(const char *big, const char *little);

// Get link to read-only string
extern const char *	strlink(const char *s);

// strlink return string when memory is gone
extern const char	strNoMemory[];

// Release link to read-only string
extern void		strunlink(const char *s);

// Print our string hash table
extern void		strprinttab(FILE *fp);

// Compute nbits-wide hash value from memory
extern unsigned long	memhash(const void *p, int len, int nbits = 32);

// Compute nbits-wide hash value from string
inline unsigned long
strhash(const char *s, int nbits = 32)
{
	return memhash(s, strlen(s), nbits);
}

// Compute nbits-wide hash value from string (caseless)
extern unsigned long	strhash_nc(const char *s, int nbits = 32);

// Class for keeping a set of sorted strings
class PStringSet {
private:
	const char **	sl;			// NULL-terminated string list
	int		maxs;			// actual array length
protected:
	int		ns;			// number of strings
public:
			PStringSet() {
				sl = 0; ns = maxs = 0;
			}
			PStringSet(const PStringSet &orig) {
				sl = 0; ns = maxs = 0;
				*this = orig;
			}
			~PStringSet() {
				Clear();
			}
			// Find index for the given string, or -1 if none
	int		Find(const char *s) const;
			// Check to see if a string is in our set
	bool		InSet(const char *s) const {
				return (Find(s) >= 0);
			}
			// Insert a string in set -- return true if added
	bool		Insert(const char *s);
			// Add a list of strings -- return # added
	int		Add(const char * const sarr[], int n) {
				if (!sarr) return 0;
				int	nnew = 0;
				while (n > 0) nnew += (int)Insert(sarr[--n]);
				return nnew;
			}
			// Add a set of strings (union) -- return # added
	int		Add(const PStringSet &ss) {
				return Add(ss.sl, ss.ns);
			}
			// Delete a string from set -- return true if deleted
	bool		Delete(const char *s);
			// Get number of strings
	int		GetSize() const {
				return ns;
			}
			// Get nth string in sorted set
	const char *	operator[](int i) const {
				if ((i < 0) | (i >= ns)) return 0;
				return sl[i];
			}
			/*
			// Get sorted string array (NULL-terminated)
	const char * const *
			GetArray() const {
				static const char *	empty = 0;
				if (!sl) return &empty;
				return sl;
			}
			*/
			// Clear the set
	void		Clear() {
				while (ns > 0) strunlink(sl[--ns]);
				if (sl) {
					free(sl); sl = 0; maxs = 0;
				}
			}
			// Copy operator
	PStringSet &	operator=(const PStringSet &src);
};

// Find near match in string set (ignores case)
extern const char *	typo_match(const char *entry, const PStringSet &sset);

#endif  // ! _PSTRINGS_H_
