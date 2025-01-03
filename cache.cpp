/*
 *  cache.cpp
 *  panlib
 *
 *  Implementation of caching classes defined in cache.h
 *
 *  Created by gward on Mon Apr 30 2001.
 *  Modified by plonghurst on Thursday Feb 15 2007.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;
#include <string.h>
#include <fcntl.h>
#include "system.h"
#include "dmessage.h"
#include "cache.h"

#ifndef USE_MMAP
#if defined(_WIN32) || defined(_WIN64)
#define USE_MMAP		0
#else
#include <sys/mman.h>
#define USE_MMAP		32000			// map files larger than this
#endif
#endif

#ifndef PDC_WRITE_TIME
#define PDC_WRITE_TIME		1.2			// write this much slower
#endif
#ifndef PDC_RELOADS_OBJ
#define PDC_RELOADS_OBJ		0.75			// approx. future reloads/object
#endif

#ifndef timeMeas
#define timeMeas		wallTime		// time measurement routine
#endif

CacheTab	CacheObject::tab;			// cache management table
unsigned long	CacheObject::nuses = 1;			// LRU counter

			// 70 Mbytes/second (initial guess) for disk read
size_t		DiskCacheObject::readBytes = 70L*1024L*1024L;
double		DiskCacheObject::readTime = 1.0;
bool		DiskCacheObject::outOfTemp = false;

// Return wall time (seconds from 1/1/1970)
double
wallTime()
{
	struct timeval	time_now;

	gettimeofday(&time_now, NULL);

	return (double)time_now.tv_sec + 1e-6*(double)time_now.tv_usec;
}

// Return CPU time (seconds from process start)
double
cpuTime()
{
	return (double)clock()*(1.0/(double)CLOCKS_PER_SEC);
}

// Our override for the (*freo)(void *) call
void
CacheObject::FreeCacheObject(MemObject *p)
{
	CacheObject *	co = (CacheObject *)p;	// may apply pointer adjustment

	if (!co->HolderCount())
		delete co;
}

// Destroy cache object, removing from cache table
// Memory must be freed in derived class
CacheObject::~CacheObject()
{
	DASSERT(RetainCount() <= 0);
	DASSERT(HolderCount() <= 1);
	tab.DeleteCacheEntry(id);
}

// Get memory object
const void *
CacheObject::GetCacheObject(size_t *siz)
{
	if (objCache == NULL) {			// have we been purged?
		double	start = timeMeas();
		if (!RestoreMemory())		// restore memory
			return NULL;
						// record time required
		restoreTime = timeMeas() - start;
		sprintf(dmessage_buf, "%s %lu byte %s",
				numRestores ? "Restored" : "Loaded",
				(unsigned long)objSize, ClassName());
		DMESG(DMCtrace, dmessage_buf);
		++numRestores;
	}
	Retain();				// increment reference counter
	if (siz != NULL)			// return size
		*siz = objSize;
	lastUse = nuses++;			// update LRU counter
	return objCache;			// return memory object
}

// Release memory object
void
CacheObject::ReleaseObject(const void *m)
{
	if ((m == NULL) | (RetainCount() <= 0))
		return;				// XXX should report error?
	DASSERT(m == objCache);
	Release();
}

// C-callable, cache-limiting memory allocation
extern "C" void
CacheMakeRoom(size_t nbytes)
{
	CacheObject::tab.MakeRoom(nbytes);
}

extern "C" void *
CacheMalloc(size_t nbytes, const char *file, int line)
{
	size_t	ob = CacheObject::tab.CheckBudget(nbytes, true);
	void *	p = malloc(nbytes);
	if (p == NULL) {
		sprintf(dmessage_buf,
		"malloc failed on %lu byte request (%u Mbytes over cache budget)",
				(unsigned long)nbytes, (unsigned int)(ob>>20));
		dmessage(DMCmemory, dmessage_buf, file, line);
	}
	return p;
}

// Free disk cache file if we allocated one
void
DiskCacheObject::FreeFile()
{
	if (objFile[0] && !origFile) {
		remove(objFile);
		objFile[0] = '\0';
		outOfTemp = false;
	}
}

// Free disk cache object after we are done with it
void
DiskCacheObject::Free()
{
	FreeFile();
	DASSERT(!InUse());
	if (objCache == NULL) {
		objSize = 0;
		return;
	}
	DASSERT(objSize != 0);
	if (origMemory)
		FreeOriginal();
#if USE_MMAP
	else if (mappedMemory)
		munmap(objCache, objSize);
#endif
	else
		Cfree(objCache);
	objCache = NULL;
	objSize = 0;
	Empty();
}

// Set up disk cache object using memory block
bool
DiskCacheObject::SetMemory(void *orig, size_t siz)
{
	Free();
	if (orig == NULL)
		return false;
	objCache = (void *)orig;
	objSize = siz;
	restoreTime = (double)siz*readTime/readBytes;
	origMemory = true;
	objFile[0] = '\0';
	reproducible = false;
	return true;
}

// Set up disk cache object using file
bool
DiskCacheObject::SetFile(const char *fname)
{
	Free();
	if (fname == NULL || strlen(fname) >= sizeof(objFile)
			|| access(fname, 4) < 0)
		return false;
	strcpy(objFile, fname);
	origFile = true;
	reproducible = false;
	return true;
}

// Restore memory for disk cache object
bool
DiskCacheObject::RestoreMemory()
{
	if (IsResident())
		return true;
	double	start = timeMeas();
	if (!objFile[0]) {			// need to (re)create?
		if (!NewOriginal())
			return false;
		origTime = timeMeas() - start;
		origMemory = true;
		return true;
	}
	int	fd = open(objFile, O_RDONLY);
	if (fd < 0) {				// check file existence
		DMESG(DMCresource, "Could not reopen disk cache file");
		return false;
	}
	SET_FD_BINARY(fd);
	if (!objSize) {				// first time? (SetFile)
		off_t	flen = lseek(fd, 0, SEEK_END);
		if (flen < 0) {
			DMESGF(DMCsystem, "Cannot seek on file '%s'", objFile);
			return false;
		}
		if (!(objSize = (size_t)flen))
			DMESGF(DMCwarning, "Empty cache file '%s'", objFile);
		else
			lseek(fd, 0, SEEK_SET);
	}
					// restore from file
#if USE_MMAP
	if (objSize >= USE_MMAP) {
		objCache = mmap(NULL, objSize, PROT_READ|PROT_WRITE, MAP_PRIVATE, fd, 0);
		if (objCache != MAP_FAILED) {
			close(fd);
			origMemory = false;
			mappedMemory = true;
			return true;
		}
		objCache = NULL;	// failure fall-back
	}
#endif
	char *	objMem = (char *)Cmalloc(objSize);
	if (objMem == NULL) {
		close(fd);
		return false;
	}
	if (read(fd, objMem, objSize) != (int)objSize) {
		Cfree(objMem);
		close(fd);
		DMESGF(DMCsystem, "Disk cache read error on '%s'", objFile);
		return false;
	}
	close(fd);
						// add to avg. statistics
	if (readBytes + objSize > readBytes) {	// prevents integer wrap
		readBytes += objSize;
		readTime += timeMeas() - start;
	}
	objCache = (void *)objMem;
	origMemory = false;
	mappedMemory = false;
	return true;
}

// Free memory for disk cache object
void
DiskCacheObject::FreeMemory()
{
	if (!IsPurgeable())
		return;
	DASSERT(objSize != 0);
	double	diskFactor = double(PDC_WRITE_TIME+PDC_RELOADS_OBJ-1 + numRestores) /
				double(PDC_RELOADS_OBJ-1 + numRestores);
	if (!objFile[0] && (!reproducible || (!outOfTemp &&
			origTime > diskFactor*objSize*readTime/readBytes))) {
		tmpnam(objFile);		// create disk copy
		int	fd = open(objFile, O_CREAT|O_WRONLY|O_EXCL, 0600);
		if (fd < 0) {
			objFile[0] = '\0';
		} else {
			SET_FD_BINARY(fd);
			if (write(fd, objCache, objSize) != (size_t)objSize) {
				FreeFile();
			} else {
				sprintf(dmessage_buf, "Wrote %u byte %s to %s",
					(unsigned int)objSize, ClassName(),
						objFile);
				DMESG(DMCtrace, dmessage_buf);
			}
			close(fd);
		}
		outOfTemp = !(objFile[0]);
		if (outOfTemp & !reproducible) {
			DMESG(DMCwarning, "FreeMemory: Out of temp disk");
			return;			// we'll lose it if we free it
		}
	}
	if (origMemory)				// free memory
		FreeOriginal();
#if USE_MMAP
	else if (mappedMemory)
		munmap(objCache, objSize);
#endif
	else
		Cfree(objCache);
	objCache = NULL;
	origMemory = false;
	mappedMemory = false;
}

// Default virtual function to get original memory (never called)
bool
DiskCacheObject::NewOriginal()
{
	DMESG(DMCassert, "Missing NewOriginal call");
	return false;
}

// Default virtual function to free original memory (never called)
void
DiskCacheObject::FreeOriginal()
{
	DMESG(DMCassert, "Missing FreeOriginal call");
}

// Comparison function used by CacheTab::CheckBudget
static int
CostCmp(const void *cop1, const void *cop2)
{
	double	diff = (**(const CacheObject **)cop1).CostEstimate() -
			(**(const CacheObject **)cop2).CostEstimate();
	if (diff > 0) return 1;
	if (diff < 0) return -1;
	return 0;
}

// Check memory budget and purge cache as necessary
size_t
CacheTab::CheckBudget(size_t nbytes, bool purgeAnyway)
{
	if (memoryBudget <= 0)
		return 0;			// unlimited budget
	const size_t	MinPurge = memoryBudget/64 + 1024*1024;
	size_t		inUse = 0;
	size_t		toPurge = 0;
	size_t		overB = 0;
	int		npurgeable = 0;
	int		i;
	ObtainReadLock();
	for (i = nentries; i--; )
		if (entry[i].refCount > 0 && entry[i].co->IsResident())
			if (entry[i].co->InUse())
				inUse += entry[i].co->GetSize();
			else {
				toPurge += entry[i].co->GetSize();
				npurgeable++;
			}
				

	if (toPurge + inUse + nbytes <= memoryBudget) {
		ReleaseReadLock();
		return 0;			// plenty o'space
	}
	if (inUse + nbytes > memoryBudget) {	// over limit?
		overB = inUse + nbytes - memoryBudget;
		if (!purgeAnyway | (toPurge == 0)) {
			ReleaseReadLock();
			return overB;
		}
	} else if (toPurge > MinPurge) {	// adjust amount
		toPurge -= memoryBudget - inUse - nbytes;
		if (toPurge < MinPurge)
			toPurge = MinPurge;
	}
	DASSERT(toPurge > 0);			// assured by above
	/*
	 *  Perform actual purge...
	 *  The method below is designed to maximize overall response
	 *  time by freeing lowest cost items first until we satisfy
	 *  the toPurge value.
	 *  Response cost is roughly the restoreTime for the object
	 *  over the number of uses that have occurred since this
	 *  object was last used:  restoreTime/(nuses - lastUse).
	 *  We sort our purgeable items by this cost and take
	 *  items until we reach our toPurge goal.
	 *  We then back up the list to reprieve objects we
	 *  don't absolutely need to meet our goal.
	 */
	CacheObject **	purgeList = new CacheObject * [npurgeable];
	int		j = 0;
	for (i = nentries; i--; )		// purgeable entries
		if (entry[i].refCount > 0 &&
				entry[i].co->IsPurgeable())
			purgeList[j++] = entry[i].co;

	DASSERT(j == npurgeable);		// sort by cost
	qsort(purgeList, npurgeable, sizeof(CacheObject *), CostCmp);

	size_t	excess = 0;			// what to purge
	for (i = 0; i < npurgeable; i++)
		if ((excess += purgeList[i]->GetSize()) >= toPurge)
			break;
	DASSERT(i < npurgeable);
	npurgeable = i + 1;
	excess -= toPurge;			// grant reprieves
	while (excess > 0 && i-- > 0)
		if (excess >= purgeList[i]->GetSize()) {
			excess -= purgeList[i]->GetSize();
			purgeList[i] = NULL;
		}
	for (i = 0; i < npurgeable; i++)	// purge
		if (purgeList[i] != NULL) {
			purgeList[i]->FreeMemory();
			sprintf(dmessage_buf, "Purged %u byte %s",
					(unsigned int)purgeList[i]->GetSize(),
					purgeList[i]->ClassName());
			DMESG(DMCtrace, dmessage_buf);
		}
	delete [] purgeList;			// clean up
	ReleaseReadLock();
	return overB;
}

// Assign new cache table entry and return id
int
CacheTab::NewCacheEntry(CacheObject *co)
{
	const int	CTblockSize = 64;
	int		i;
	DASSERT(co != NULL);
	ObtainWriteLock();				// lock table
	for (i = nentries; i--; )
		if (entry[i].refCount <= 0)
			break;				// found free entry
	if (i < 0) {					// else grow table
		if (!nentries)
			entry = (CacheEntry *)malloc(
					CTblockSize*sizeof(CacheEntry));
		else
			entry = (CacheEntry *)realloc(entry,
					(nentries+CTblockSize)*sizeof(CacheEntry));
		DASSERT(entry != NULL);
		for (i = nentries+CTblockSize-1; i > nentries; i--) {
			entry[i].refCount = 0;
			entry[i].co = NULL;
		}
		nentries += CTblockSize;
	}
	entry[i].refCount = 1;				// assign new entry
	entry[i].co = co;
	ReleaseWriteLock();				// unlock table
	return i;
}

// Find cache object that contains this bit of memory
CacheObject *
CacheTab::GetCacheObject(const void *ptr)
{
	if (ptr == NULL)
		return NULL;
	ObtainReadLock();
	int	i = nentries;
	while (i--)
		if (entry[i].refCount > 0 && entry[i].co->IsResident()) {
			size_t		sz;
			const char *	op = (const char *)
						entry[i].co->GetCacheObject(&sz);
			entry[i].co->ReleaseObject(op);
			if ((ptr >= op) & (ptr < op + sz))
				break;
		}
	if (i < 0)
		return NULL;
	CacheObject *	co = entry[i].co;		// get pointer
	entry[i].refCount++;				// add reference
	ReleaseReadLock();
	return co;
}

// Return pointer to cache object
CacheObject *
CacheTab::GetCacheObject(int id)
{
	ObtainReadLock();
	DASSERT((id >= 0) & (id < nentries) && entry[id].refCount > 0);
	CacheObject *	co = entry[id].co;		// get pointer
	entry[id].refCount++;				// add reference
	ReleaseReadLock();
	return co;
}

// Release cache object pointer
void
CacheTab::ReleaseCacheObject(CacheObject *co)
{
	ObtainReadLock();
	int	id = co->GetCacheID();
	DASSERT((id >= 0) & (id < nentries) && entry[id].refCount > 1);
	entry[id].refCount--;				// remove reference
	ReleaseReadLock();
}

// Delete cache table entry
void
CacheTab::DeleteCacheEntry(int id)
{
	ObtainWriteLock();
	DASSERT((id >= 0) & (id < nentries) && entry[id].refCount <= 1);
	entry[id].co = NULL;
	entry[id].refCount = 0;
	ReleaseWriteLock();
}
