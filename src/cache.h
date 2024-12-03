/*
 *  cache.h
 *  panlib
 *
 *  Depends on <string.h>
 *
 *  Base class for caching objects in memory and on disk
 *
 *  Created by gward on Mon Apr 30 2001.
 *  Modified by plonghurst on Thursday Feb 15 2007.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#ifndef _CACHE_H_
#define _CACHE_H_

#ifndef _MEMOBJECT_H_
#include "memobject.h"
#endif

extern double	wallTime();		// return wall time (seconds from 1/1/1970)
extern double	cpuTime();		// return CPU time (seconds used by process)

class CacheObject;

// Table for managing cached object memory (should be only one such table)
class CacheTab {
private:
	int		nentries;	// number of cache entries
	struct CacheEntry {
		int		refCount;	// reference count
		CacheObject *	co;		// cache object
	} *		entry;		// cache entry table
	int		lockState;	// 0==unlocked, -1==write, >0==read
	static void	Snore() {
				(void)time(0);
			}
	void		ObtainWriteLock() {
				while (lockState != 0) Snore();
				lockState = -1;
			}
	void		ReleaseWriteLock() {
				lockState = 0;
			}
	void		ObtainReadLock() {
				while (lockState < 0) Snore();
				++lockState;
			}
	void		ReleaseReadLock() {
				if (lockState > 0) --lockState;
			}
public:
	size_t		memoryBudget;	// allowed memory budget (bytes)
			CacheTab(size_t maxbytes = 0) {
				nentries = 0;
				entry = NULL;
				lockState = 0;
				memoryBudget = maxbytes;
			}
			~CacheTab() {
				if (nentries > 0) free(entry);
			}
			// Free memory to get under budget, or say how short
	size_t		CheckBudget(size_t nbytes, bool purgeAnyway = false);
			// Make room for allocating nbytes
	void		MakeRoom(size_t nbytes) {
				CheckBudget(nbytes, true);
			}
			// Create a new cache table entry
	int		NewCacheEntry(CacheObject *co);
			// Get cache object by its ID
	CacheObject *	GetCacheObject(int id);
			// Get cache object by a pointer to its memory
	CacheObject *	GetCacheObject(const void *ptr);
			// Release cache object
	void		ReleaseCacheObject(CacheObject *co);
			// Delete a cache entry from table
	void		DeleteCacheEntry(int id);
};

// Abstract class for cached (recreatable) memory objects
// Note: MemObject::nref counts memory r.t. class object references
// WARNING: Do not incorporate as a member object -- call new & HolderRelease()
class CacheObject : protected MemObject {
private:
	int		id;		// identifier for this object
	int		nhold;		// # references to holder
	unsigned long	lastUse;	// tick value at last use
	static unsigned long
			nuses;		// LRU counter
protected:
	int		numRestores;	// # of times object restored
	double		restoreTime;	// time req'd for last restore
	size_t		objSize;	// memory block size (bytes)
	void *		objCache;	// memory block (NULL if free)
			// Our override for the (*freo)() call
	static void	FreeCacheObject(MemObject *p);
public:
	static CacheTab	tab;		// table of cached objects
			CacheObject() {
				id = tab.NewCacheEntry(this);
				freo = &CacheObject::FreeCacheObject;
				nref = 0; nhold = 1;
				objSize = 0; objCache = NULL;
				numRestores = 0; restoreTime = 0;
			}
	virtual		~CacheObject();
			// Get name of subclass
	virtual const char *
			ClassName() const {
				return "CacheObject";
			}
			// Someone wants to keep holder
	int		HolderRetain() {
				return ++nhold;
			}
			// Replacement for delete operator
	virtual void	HolderRelease() {
				// XXX only works if created using "new"
				nhold -= (nhold > 0);
				if (!nhold & !nref) delete this;
			}
			// How many people are retaining us?
	int		HolderCount() const {
				return nhold;
			}
			// Get cache table entry ID
	int		GetCacheID() const {
				return id;
			}
			// Get size of this object if known
	size_t		GetSize() const {
				return objSize;
			}
			// Check if object is memory-resident
	bool		IsResident() const {
				return (objCache != NULL);
			}
			// Is object currenty in use?
	int		InUse() const {
				return RetainCount();
			}
			// Can object be purged?
	bool		IsPurgeable() const {
				return (IsResident() && !InUse());
			}
			// What is the relative cost of purging this object?
	double		CostEstimate() const {
				return restoreTime/double(nuses - lastUse);
			}
			// Check out the object, regenerating if not resident
			// PL: Renamed to stop clash with MS Windows function
	const void *	GetCacheObject(size_t *siz = NULL);
			// Release pointer to this object
	void		ReleaseObject(const void *m);
/*
			// Release reference to any cache memory
	static bool	ReleaseReference(const void *ptr);
*/
			// (Re)create this object
	virtual bool	RestoreMemory() = 0;
			// Free this object
	virtual void	FreeMemory() = 0;
};

#define	gCacheSize	CacheObject::tab.memoryBudget

// Class for cached objects with separate file backing store
class DiskCacheObject : public CacheObject {
private:
	static bool	outOfTemp;		// out of temp space?
	bool		origMemory;		// using original's memory?
	bool		mappedMemory;		// using mmap()?
	bool		origFile;		// using original's file?
	double		origTime;		// time to (re)create original
	char		objFile[1024];		// backing store file name
	static size_t	readBytes;		// total disk bytes read
	static double	readTime;		// total disk read time
	void		Empty() {
				origMemory = false;
				mappedMemory = false;
				origFile = false;
				origTime = 0;
				objFile[0] = '\0';
				reproducible = true;
			}
protected:
	bool		reproducible;		// NewOriginal works
			// Free disk file if we allocated one
	void		FreeFile();
			// Free disk cache object after we are done with it
	void		Free();
			// Set up disk cache object using memory block
	bool		SetMemory(void *orig, size_t siz);
			// Set up disk cache object using file
	bool		SetFile(const char *fname);
			// (Re)create original object
	virtual bool	NewOriginal();
			// Free original object
	virtual void	FreeOriginal();
public:
			DiskCacheObject() {
				Empty();
			}
			DiskCacheObject(void *orig, size_t siz) {
				Empty();
				SetMemory(orig, siz);
			}
			DiskCacheObject(const char *fname) {
				Empty();
				SetFile(fname);
			}
	virtual		~DiskCacheObject() {
				Free();
			}
			// Get name of subclass
	virtual const char *
			ClassName() const {
				return "DiskCacheObject";
			}
	virtual bool	RestoreMemory();
	virtual void	FreeMemory();
};

// Reference to cached object (handy volatile reference using cache id)
class CacheReference {
private:
	CacheObject *	co;		// cache object pointer
public:
	const void *	cacheMemory;	// pointer to cache memory
	size_t		cacheSize;	// size of cache memory
			CacheReference(int id) {
				co = CacheObject::tab.GetCacheObject(id);
				if (co == NULL) {
					cacheMemory = NULL;
					cacheSize = 0;
					return;
				}
				cacheMemory = co->GetCacheObject(&cacheSize);
			}
			~CacheReference() {
				co->ReleaseObject(cacheMemory);
				CacheObject::tab.ReleaseCacheObject(co);
			}
};

extern "C" {				// C-callable cache allocation
extern void		CacheMakeRoom(size_t nbytes);
extern void *		CacheMalloc(size_t nbytes, const char *file, int line);
#define Cmalloc(n)	CacheMalloc(n, __FILE__, __LINE__)
#define Cfree		free
}

#endif	// !_CACHE_H_
