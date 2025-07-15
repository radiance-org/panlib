/*
 *  astack.h
 *  pan
 *
 *  Template class for LIFO stack with block allocation (e.g., for flood fills)
 *  Random access and FIFO also supported, but not quite as efficiently
 *
 *  Created by Greg Ward on 2/21/18.
 *  Copyright 2018 Anyhere Software. All rights reserved.
 *
 */

#ifndef ASTACK_L
#define	ASTACK_L	512
#endif

template <class T> struct AStackBlock {
	AStackBlock<T> *	prev;		// previous stack block
	int			ne;		// number of assigned entries
	T			ent[ASTACK_L];	// block entries
				AStackBlock(AStackBlock<T> *pp = 0) {
					ne = 0; prev = pp;
				}
				~AStackBlock() {
					delete prev;
				}
	bool			Full() const {
					return (ne >= ASTACK_L);
				}
};

template <class T> class AStack {
	AStackBlock<T> *	top;		// top stack block
	int			start;		// starting position for FIFO
	int			len;		// active stack length
	static T		dummy;
public:
			AStack() {
				top = 0; start = len = 0;
			}
			~AStack() {
				delete top;
			}
			// Clear the stack
	void		Clear() {
				delete top; top = 0; start = len = 0;
			}
			// Push object onto LIFO top, return new stack length
	int		Push(const T &t) {
				if (!top || top->Full())
					top = new AStackBlock<T>(top);
				top->ent[top->ne++] = t;
				return ++len;
			}
			// Copy LIFO top if any
	bool		Peek(T *tp = 0) const {
				if (!top) return false;
				if (tp) *tp = top->ent[top->ne-1];
				return true;
			}
			// Pop last entry off LIFO top unless exhausted
	bool		Pop(T *tp = 0) {
				if (!top) return false;
				top->ne--;
				if (tp) *tp = top->ent[top->ne];
				if (--len <= 0) Clear();
				else if (!top->ne) {
					AStackBlock<T> *	tdone = top;
					top = top->prev; tdone->prev = 0;
					delete tdone;
				}
				return true;
			}
			// Pull first entry from bottom (FIFO mode)
	bool		Pull(T *tp = 0);
			// Get active stack length
	int		Length() const {
				return len;
			}
			// Look up ith entry from bottom (first==0)
	const T &	operator[](int i) const {
				if ((i < 0) | (i >= len)) return dummy;
				i = len - i;
				const AStackBlock<T> *	ptop = top;
				while (i > ptop->ne) {
					i -= ptop->ne; ptop = ptop->prev;
				}
				return ptop->ent[ptop->ne - i];
			}
			// Writable stack entry
	T &		operator[](int i) {
				if ((i < 0) | (i >= len)) return dummy;
				i = len - i;
				AStackBlock<T> *	ptop = top;
				while (i > ptop->ne) {
					i -= ptop->ne; ptop = ptop->prev;
				}
				return ptop->ent[ptop->ne - i];
			}
};

// Pull first entry from bottom (FIFO mode)
template <class T>
bool
AStack<T>::Pull(T *tp)
{
	if (!top)
		return false;
	if (start < ASTACK_L-1) {		// normal case
		if (tp) *tp = (*this)[0];
		++start;
	} else if (top->prev) {			// losing first block
		AStackBlock<T> *	pend = top;
		while (pend->prev->prev)
			pend = pend->prev;
		if (tp) *tp = pend->prev->ent[ASTACK_L-1];
		delete pend->prev; pend->prev = 0;
		start = 0;
	} else if (tp) {			// corner case on last block
		*tp = top->ent[ASTACK_L-1];
	}
	if (--len <= 0) Clear();		// last entry?
	return true;
}

// Commonly used stack types

struct IxyPoint {
	int		x, y;
			IxyPoint() {}
			IxyPoint(int ix, int iy) {
				x = ix; y = iy;
			}
};
typedef AStack<IxyPoint>	IxyStack;

struct IxyzPoint {
	int		x, y, z;
			IxyzPoint() {}
			IxyzPoint(int ix, int iy, int iz) {
				x = ix; y = iy; z = iz;
			}
};
typedef AStack<IxyzPoint>	IxyzStack;
