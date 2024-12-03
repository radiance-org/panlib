/*
 *  ablockarray.h
 *  pan
 *
 *  Template class for efficient allocation of objects in ABLKSIZ-length blocks.
 *  Faster than individual "new" operations, and allows for smaller address spaces.
 *  Overheads other than unused entries are minimal, and a free list is maintained
 *  at a cost of ~1 bit/object.  This code assumes that you are allocating tens of
 *  thousands to billions of objects, not just a handful.  The ABLKBITS macro
 *  allows some tuning, but block allocation doesn't make sense for small arrays.
 *  Also, sizeof(T) should probably be greater than 4 bytes, otherwise there is
 *  no savings from trading the contents for a 32-bit index to an array.
 *  (Of course, your copy of an index may be less than 32 bits if you know your array
 *  will fit in a smaller address space, but it's a game of diminishing returns.)
 *
 *  32-bit unsigned indexing starts at 1, reserving 0 as a NULL-equivalent.
 *  Adjacent entries are not necessarily adjacent in memory, or even ordered.
 *  Once an index is allocated, its memory location is fixed until recycled, with
 *  contents unaltered until FreeSpace() or RecycleA() is called or the array size
 *  is reduced.  Entries may be appended to the array by calling AllocA(nent).
 *  A maximum allowed index may be set with SetMaxIndex(), which defaults to 0xffffffff.
 *  FreeSpace() and RecycleA() free entire blocks at a time once they are empty.
 *  Tail blocks will also be freed as their entries are recycled, unless the
 *  "freeTail" argument passed to Recycle() is false.  At any time, the Last()
 *  function will report the largest allocated index.
 *
 *  Beware of pathological behavior when calls alternately allocate and recycle the
 *  last few entries, as this can result in excessive "new" and "delete" operations
 *  if Last() is near a block boundary (i.e., divisible by ABLKSIZ).  If you
 *  interpserse many calls to Alloc() and Recycle(), set "freeTail" to false
 *  to avoid performance penalties.  You can always call FreeSpace() at a later
 *  point to release any unused blocks.
 *
 *  Continuous allocation is supported, so (array[++i] = value) works indefinitely,
 *  but only if the requested index is not more than one past Last().  Recycling the
 *  last entry similarly reduces the array size by at least one, and possibly much
 *  more if the released entry is connected to a block of unassigned entries.
 *
 *  This is not as general as a container class, since constructors/destructors
 *  are only called en masse as blocks are allocated/freed.  Thus, recycling an
 *  entry does not immediately call the object's destructor or free its memory.
 *  The caller can and probably should manually free object data and return to
 *  its initialized state before calling Recycle(), or  assign a clearing function
 *  to the ABclear macro, which takes no arguments.
 *
 *  Most routines are O(1).  The notable exceptions are AllocA() and RecycleA(),
 *  which are O(M) for M allocations but much faster than serial calls, and Index(),
 *  which is O(log(N)) for N allocated entries.
 *
 *  Copyright (c) 2023 Anyhere Software. All rights reserved.
 *
 */

#ifndef _ABLOCKARRAY_H_
#define _ABLOCKARRAY_H_

#include <stdlib.h>
#include <string.h>
#include "abitmap.h"

#ifndef ABLKBITS
#define ABLKBITS	10		// controls block allocation size
#endif
#define ABLKSIZ		(1<<ABLKBITS)	// entries per allocated block

// A block array block (not meant to be used directly)
template <class T>
struct ABABlock {
	ABitMap		useMap;		// which entries are in use?
	mutable uint32	firstFree;	// first potentially available entry
	T		entry[ABLKSIZ];	// block of allocated class objects
			ABABlock() : useMap(ABLKSIZ) {
				firstFree = 0;
			}
			// Is the block empty?
	bool		IsEmpty() const {
				return (!firstFree && useMap.Find() >= ABLKSIZ);
			}
			// Is the block full?
	bool		IsFull() const {
				if (firstFree < ABLKSIZ)
					firstFree = useMap.Find(firstFree, false);
				return (firstFree >= ABLKSIZ);
			}
			// Get number allocated
	int		NAlloc() const {
				if (IsFull()) return ABLKSIZ;
				if (IsEmpty()) return 0;
				return useMap.SumTotal();
			}
			// Get first available slot
	int		Alloc() {
				if (IsFull()) return -1;
				useMap.Set(firstFree);
				return (int)firstFree++;
			}
			// Mark indicated slot as free
	bool		Free(int ei) {
				if (ei < 0 || !useMap.Check(ei)) return false;
#ifdef ABclear				// class clearing function
				entry[ei].ABclear();
#endif
				useMap.Reset(ei);
				if (firstFree > ei) firstFree = ei;
				return true;
			}
};

// A block array object allocator (not thread-safe if global)
template <class T>
class ABlockArray {
	ABABlock<T> **	block;		// pointers to allocation blocks
	int *		blksort;	// blocks sorted by address for Index()
	int		ablock, basiz;	// next allocation block, array size
	int		nballoc;	// number of blocks actually allocated
	uint32		lastE;		// largest allocated index
	uint32		maxIndex;	// caller-specified allocation limit
			// allows for more (or fewer) allocated blocks
	bool		ResizePointers(int nblks = -1);
			// require block allocation
	void		ReqBlock(int b);
			// delete an allocated block
	void		DelBlock(int b);
public:
			ABlockArray() {
				block = 0; blksort = 0;
				nballoc = ablock = basiz = 0;
				lastE = 0;
				maxIndex = 0xffffffff;
			}
			ABlockArray(uint32 siz, mx = 0) {
				block = 0; blksort = 0;
				nballoc = ablock = basiz = 0;
				lastE = 0;
				if (!mx) mx = 0xffffffff;
				else if (mx < siz) mx = siz;
				maxIndex = mx;
				AllocA(siz);
			}
			~ABlockArray() {
				Clear();
			}
			// Clear everything, return to empty state (maxIndex unchanged)
	void		Clear();
			// Find index for an allocated pointer
	uint32		Index(const T *p) const;
			// Allocate and return first unused entry
	uint32		Alloc();
			// Grow end of array, allocating nent new objects
	uint32		AllocA(uint32 nent);
			// How many entries are currently allocated?
	uint32		NAlloc() const;
			// What is the current fill factor (allocated/max.)?
	double		Occupancy() const {
				if (!lastE) return .0;
				return (double)NAlloc() / (double)lastE;
			}
			// Free (delete) empty blocks and unused pointers
	uint32		FreeSpace();
			// Query largest allocated index
	uint32		Last() const {
				return lastE;
			}
			// Query maximum index
	uint32		GetMaxIndex() const {
				return maxIndex;
			}
			// (Re)set maximum index
	void		SetMaxIndex(uint32 mx = 0) {
				if (!mx) { maxIndex = 0xffffffff; return; }
				if (lastE > mx) RecycleA(mx+1, lastE-mx);
				maxIndex = mx;
			}
			// Recycle the specified entry
	bool		Recycle(uint32 i, bool freeTail = true);
			// Same but based on pointer
	bool		Recycle(T *p, bool freeTail = true) {
				return Recycle(Index(p), freeTail);
			}
			// Recycle an array of entries (always frees what it can)
	uint32		RecycleA(uint32 i, uint32 n);
/*
			// Recycle an array of entries
	uint32		RecycleA(uint32 i, uint32 n) {
				if (!i | (i > lastE)) return 0;
				if (n > lastE+1-i) n = lastE+1-i;
				uint32	cnt = 0;
				while (n--) cnt += Recycle(i++);
				return cnt;
			}
*/
			// Allocate entry as pointer (return index if requested)
	T *		Alloc(uint32 *ip) {
				uint32	i = Alloc();
				if (ip) *ip = i;
				if (!i--) return 0;
				return &block[i >> ABLKBITS]->entry[i & ABLKSIZ-1];
			}
			// Get pointer to read-only entry if allocated
	const T *	Get(uint32 i) const {
				if (!i--) return 0;
				const int	b = i >> ABLKBITS;
				if (b >= basiz || !block[b]) return 0;
				i &= ABLKSIZ-1;
				if (!block[b]->useMap.Check(i)) return 0;
				return &block[b]->entry[i];
			}
			// Get pointer to writeable entry (allocates if alloc)
	T *		Entry(uint32 i, bool alloc = false);
			// Read-only array operator
	const T &	operator[](uint32 i) const {
				static const T	dummy;
				const T *	p = Get(i);
				if (!p) return dummy;
				return *p;
			}
			// Access writeable entry with array operator
	T &		operator[](uint32 i) {
				static T	dummy;
				T *	p = Entry(i,true);
				if (!p) return dummy;
				return *p;
			}
};

template <class T>
void
ABlockArray<T>::Clear()
 {
 	ablock = 0;
	lastE = 0;
	if (!block | !blksort) {
		nballoc = basiz = 0;
 		return;
	}
	while (nballoc > 0)
		delete block[blksort[--nballoc]];
	free(block); block = 0;
	free(blksort); blksort = 0;
	basiz = 0;
}

template <class T>
uint32
ABlockArray<T>::Index(const T *p) const
{
	int	tupper = nballoc, tlower = 0;
	int	b = -1;		// binary search for block
	while (tupper != tlower) {
		int		t = (tupper + tlower) >> 1;
		const T *	earr = block[blksort[t]]->entry;
		if (p < earr) {
			tupper = t;
		} else if (p >= earr + ABLKSIZ) {
			if (t == tlower) break;
			tlower = t;
		} else {	// found entry's block
			b = blksort[t];
			break;
		}
	}
	if (b < 0)
		return 0;	// failure to locate in any block

	int	ei = p - block[b]->entry;

	if (!block[b]->useMap.Check(ei))
		return 0;	// not on record as allocated

	return ((uint32)b << ABLKBITS) + ei + 1;
}

template <class T>
bool
ABlockArray<T>::ResizePointers(int nblks)
{
	if (nblks < 0)			// default is to grow some amount
		nblks = basiz + (basiz >> 1) + 32;
	if (nblks > maxIndex>>ABLKBITS) {
		nblks = maxIndex >> ABLKBITS;
		if (nblks <= basiz)
			return false;	// end of addressable range
	}
	if (nblks == basiz)
		return (block != 0);
	if (!nblks) {
		free(block); block = 0;
		free(blksort); blksort = 0;
	} else if (!block) {
		blksort = (int *)malloc(nblks*sizeof(int));
		basiz = 0;
		block = (ABABlock<T> **)malloc(nblks*sizeof(ABABlock<T> *));
	} else {
		blksort = (int *)realloc(blksort, nblks*sizeof(int));
		block = (ABABlock<T> **)realloc(block, nblks*sizeof(ABABlock<T> *));
	}
	if (!block | !blksort) {	// catastrophic unless requested
		if (block) { free(block); block = 0; }
		if (blksort) { free(blksort); blksort = 0; }
		nballoc = ablock = basiz = 0;
		lastE = 0;
		return false;
	}
	if (nblks > basiz)		// clear new pointers
		memset(block+basiz, 0, (nblks-basiz)*sizeof(ABABlock<T> *));
	else if (ablock > nblks)
		ablock = nblks;
	basiz = nblks;
	return true;
}

template <class T>
void
ABlockArray<T>::ReqBlock(int b)
{					// assumes b is in range
	if (block[b])
		return;
	block[b] = new ABABlock<T>;
	int	bi = ++nballoc;		// insert in sorted list
	while (--bi > 0)
		if (block[blksort[bi-1]] < block[b])
			break;
	memmove(blksort+bi+1, blksort+bi, (nballoc-1-bi)*sizeof(int));
	blksort[bi] = b;
}

template <class T>
void
ABlockArray<T>::DelBlock(int b)
{					// assumes b is in range
	delete block[b]; block[b] = 0;
	int	bi = nballoc--;		// find and remove from sorted list
	while (bi-- > 0)
		if (blksort[bi] == b)
			break;
	if (bi < 0) abort();		// code consistency error!
	memmove(blksort+bi, blksort+bi+1, (nballoc-bi)*sizeof(int));
}

template <class T>
uint32
ABlockArray<T>::FreeSpace()
{
	uint32	nfreed = 0;
	int	bmax = 0;
	int	bi = nballoc;		// free empty blocks from end
	while (bi-- > 0)
		if (block[blksort[bi]]->IsEmpty()) {
			DelBlock(blksort[bi]);
			nfreed += ABLKSIZ;
		} else if (blksort[bi] >= bmax)
			bmax = blksort[bi] + 1;

	ResizePointers(bmax);		// resize/free pointer array
	return nfreed;
}

template <class T>
uint32
ABlockArray<T>::Alloc()
{
	while (ablock < basiz || ResizePointers()) {
		ReqBlock(ablock);	// allocate block if necessary
		int	ei = block[ablock]->Alloc();
		if (ei >= 0) {		// allocated entry from this block?
			uint32	i = (uint32)ablock << ABLKBITS | ei;
			if (++i > maxIndex) {
				block[ablock]->Free(ei);
				break;	// hit hard stop
			}
			if (i > lastE) lastE = i;
			return i;	// got new entry
		}
		++ablock;		// else move on to next
	}
	return 0;
}

template <class T>
uint32
ABlockArray<T>::AllocA(uint32 nent)
{
	if (!nent)
		return 0;
	if (nent > maxIndex - lastE)
		return 0;		// past addressable range
	const uint32	ebeg = lastE+1;
	if (nent <= ABLKSIZ) {		// short allocation
		nent += lastE;
		while (lastE < nent)
			if (!Entry(lastE+1, true))
				return 0;
		return ebeg;
	}
	if (lastE)
		nent += lastE;
	else if (!Alloc())		// initialize if necessary
		return 0;
	const int	bi = (nent-1) >> ABLKBITS;
	int		b = (lastE-1) >> ABLKBITS;
	int		i = lastE-1 & ABLKSIZ-1;
					// flag entries in first block
	block[b]->useMap.ClearBits(i+1, ABLKBITS-1-i, true);
					// check pointer array capacity
	if (bi >= basiz && !ResizePointers(bi+1))
		return 0;
	while (++b <= bi) {		// initialize new blocks as needed
		ReqBlock(b);
		if (b < bi)
			block[b]->useMap.ClearBitMap(true);
	}
	lastE = nent;			// flag used up to new end
	block[bi]->firstFree = (lastE-1 & ABLKSIZ-1) + 1;
	block[bi]->useMap.ClearBits(0, block[bi]->firstFree, true);
	return ebeg;
}

template <class T>
uint32
ABlockArray<T>::RecycleA(uint32 i, uint32 n)
{
	if (!i | (i > lastE) | !n)
		return 0;
	if (i+n > lastE+1)
		n = lastE+1 - i;
	uint32		cnt = 0;
	if ((i == 1) & (n == lastE)) {	// free everything?
		cnt = NAlloc();
		Clear();
		return cnt;
	}
	if (n <= ABLKSIZ) {		// short array?
		while (n--)
			cnt += Recycle(i++,true);
		return cnt;
	}
	const int	b1 = (ABLKSIZ-2+i) >> ABLKBITS;
	const int	bi = (i+n-2) >> ABLKBITS;
	int		b;
	for (b = b1; b < bi; b++)	// free middle blocks
		if (block[b]) {
			cnt += block[b]->NAlloc();
			DelBlock(b);
		}
	if (ablock > b1)		// first available
		ablock = b1;
	if (i+n == lastE+1) {		// free to end?
		if (block[bi]) {
			cnt += block[bi]->NAlloc();
			DelBlock(bi);
		}
		lastE = b1 << ABLKBITS;
		n = lastE+1-i;
	}
	while (n--)			// free head/tail
		if ((n >= ABLKSIZ-1) & !(i-1 & ABLKSIZ-1)) {
			i += ABLKSIZ;	// freed middle block
			n -= ABLKSIZ-1;
		} else
			cnt += Recycle(i++,true);
	return cnt;
}

template <class T>
uint32
ABlockArray<T>::NAlloc() const
{
	uint32		totass = 0;
	int		bi = nballoc;
	while (bi-- > 0)
		totass += block[blksort[bi]]->NAlloc();
	return totass;
}

template <class T>
T *
ABlockArray<T>::Entry(uint32 i, bool alloc)
{
	if (!i-- || (i > lastE) | (i >= maxIndex))
		return 0;
	const int	b = i >> ABLKBITS;
	if (!alloc && (b >= basiz || !block[b]))
		return 0;
	if (b >= basiz && !ResizePointers())
		return 0;
	ReqBlock(b);			// allocates if not already
	const int	ei = i & ABLKSIZ-1;
	if (alloc)
		block[b]->useMap.Set(ei);
	else if (!block[b]->useMap.Check(ei))
		return 0;
	lastE += (i == lastE);
	return &block[b]->entry[ei];
}

template <class T>
bool
ABlockArray<T>::Recycle(uint32 i, bool freeTail)
{
	if (!i | (i > lastE))		// counting on lastE to be right
		return false;
	if (i-- == lastE) {
		do {			// find new last allocation
			if (!--lastE) {
				Clear();
				return true;
			}
		} while (!Get(lastE));
	}
	const int	bi = i >> ABLKBITS;
	if (bi >= basiz || !block[bi])
		return false;
	i &= ABLKSIZ-1;
	if (!block[bi]->Free(i))	// fails if already freed
		return false;
	if (ablock > bi)		// track first block with free space
		ablock = bi;
	int 	b = (lastE-1) >> ABLKBITS;
	if (ablock > b)
		ablock = b;
	if (!freeTail)			// hold off deletions?
		return true;
	while (++b <= bi)		// else delete empty tail blocks
		if (block[b] && block[b]->IsEmpty())
			DelBlock(b);
	return true;
}

#endif /* _ABLOCKARRAY_H_ */
