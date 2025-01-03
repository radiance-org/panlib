/*
 *  memobject.c
 *  panlib
 *
 *  Created by Gregory Ward on Thu Oct 31 2002.
 *  Copyright (c) 2002 Anyhere Software. All rights reserved.
 *
 *  Simple memory reference manager.
 */

#include <stdio.h>
#include <stdlib.h>
#include "dmessage.h"
#include "memobject.h"

/* Basic free routine */
void
SysFree(MemObject *p)
{
	free(p);
}

/* Basic memory object allocator */
MemObject *
MobjAlloc(size_t n)
{
	MemObject *	mo = (MemObject *)malloc(sizeof(MemObject) + n);

	if (mo != NULL) {
		mo->nref = 1;
		mo->freo = SysFree;
	}
	return mo;
}

/* Memory object allocator with error reporting */
MemObject *
MobjAllocReport(size_t n, const char *fn, int ln)
{
	MemObject *	mo = MobjAlloc(n);
	
	if (mo == NULL) {
		sprintf(dmessage_buf, "malloc() failed on %lu byte request",
				(unsigned long)n);
		dmessage(DMCmemory, dmessage_buf, fn, ln);
	}
	return mo;
}

/* Memory object allocator that returns pointer to user memory */
void *
MobjMemAlloc(size_t n)
{
	MemObject *	mo = MobjAlloc(n);
	
	if (mo == NULL)
		return NULL;
	return MobjMem(mo);
}

/* Memory object allocator with reporting that returns user memory */
void *
MobjMemAllocReport(size_t n, const char *fn, int ln)
{
	MemObject *	mo = MobjAllocReport(n, fn, ln);

	if (mo == NULL)
		return NULL;
	return MobjMem(mo);
}

/* Free user memory pointer attached to memory object */
void
MOMrelease(void *p)
{
	MemObject *	mo;

	if (p == NULL)
		return;
	mo = MobjObj(p);
	if (--(mo->nref) > 0)
		return;
	if (mo->nref < 0) {
		DMESGF(DMCparameter, "Negative reference count (%d)", mo->nref);
		return;
	}
	if (mo->freo != NULL)
		(*mo->freo)(mo);
}
