/*
 *  memobject.h
 *  panlib
 *
 *  Created by Gregory Ward on Thu Oct 31 2002.
 *  Copyright (c) 2002 Anyhere Software. All rights reserved.
 *
 *  A simple shared memory reference allocator.
 *  Use MemObject as a public base class to provide for
 *  shared object references.  Use the provided Retain()
 *  and Release() functions to link and delete objects.
 *  (The MobjRetain() and MobjRelease() macros also
 *  work from C or C++ code for derived objects.)
 */

#ifndef _MEMOBJECT_H_
#define _MEMOBJECT_H_

#include <sys/types.h>

#ifdef __cplusplus

/* Base object sharing class (C++ definiton) */
struct MemObject {
	static void
		FreeObj(MemObject *p) {
			delete p;
		}
	int	nref;			/* reference count -- free at 0 */
	void	(*freo)(MemObject *);	/* free allocated object */
		MemObject() {
			nref = 1; freo = &MemObject::FreeObj;
		}
	int	Retain() {
			return ++nref;
		}
	void	Release() {
			if (!--nref && freo) (*freo)(this);
		}
	int	RetainCount() const {
			return nref;
		}
};

#define MobjRetain(mo)		((mo)->Retain())
#define MobjRelease(mo)		((mo)->Release())

extern "C" {

#else	/* ! __cplusplus */

/* Struct for sharing memory references (C definiton) */
typedef struct _MemObject {
	int	nref;				/* reference count -- free when 0 */
	void	(*freo)(struct _MemObject *);	/* free our allocated memory */
} MemObject;

#define	MobjRetain(mo)		(++((mo)->nref))
#define	MobjRelease(mo)		if (!--((mo)->nref) && (mo)->freo) \
					(*(mo)->freo)(mo); else

/* Basic free routine */
extern void			SysFree(MemObject *);

#endif	/* ! __cplusplus */

#define	MobjMem(mo)		((mo) ? (void *)((mo)+1) : (void *)0)
#define	MobjObj(p)		((p) ? (MemObject *)(p)-1 : (MemObject *)0)
#define	MobjMemRetain(p)	MobjRetain(MobjObj(p))
#define MobjMemRelease(p)	MobjRelease(MobjObj(p))

#define	MOalloc(n)		MobjAllocReport(n, __FILE__, __LINE__)
#define MOMalloc(n)		MobjMemAllocReport(n, __FILE__, __LINE__)

extern MemObject *		MobjAlloc(size_t n);
extern MemObject *		MobjAllocReport(size_t n, const char *fn, int ln);

extern void *			MobjMemAlloc(size_t n);
extern void *			MobjMemAllocReport(size_t n, const char *fn, int ln);

extern void			MOMrelease(void *p);

#ifdef __cplusplus
}
#endif

#endif	/* ! _MEMOBJECT_H_ */
