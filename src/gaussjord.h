/*
 *  gaussjord.h
 *
 *  Solve set of linear equations using Gauss-Jordan elimination.
 */

#ifndef _GAUSSJORD_H_
#define _GAUSSJORD_H_

template <class R>
inline R
absval(R v)
{
	return (v >= R(0)) ? v : -v;
}

// Use Gauss-Jordan elimination to solve matrix equation Ax==b
template <class R>
bool
GJsolveLinearEq(
			const int dim,
			const R *a,	// a[dim][dim] (a[row][col])
			const R *b,	// b[dim]
			R *x)		// returned solution x[dim]
{
	const int n1=dim+1;
	int i,j,k,pivot=0;
	
	if (dim <= 0)
		return false;

	R *m = new R[dim*n1];

	for(i=0;i<dim;i++)
	{
		int in0,in1;

		in0=i*dim;
		in1=i*n1;
		for(j=0;j<dim;j++)
			m[in1+j]=a[in0+j];
		m[in1+dim]=b[i];
	}

	for(k=0;k<dim;k++)
	{
		int kn1=k*n1;

		R pvval=0;
		for(i=k;i<dim;i++)
			if(absval(m[i*n1+k])>absval(pvval))
			{
				pivot=i;
				pvval=m[i*n1+k];
			}
		
		if (pvval == 0) {
			delete [] m;
			return false;
		}

		int pn1=pivot*n1;
		for(j=0;j<n1;j++)
		{
			R t=m[pn1+j];
			m[pn1+j]=m[kn1+j];
			m[kn1+j]=t;
		}

		for(j=k;j<n1;j++)
			m[kn1+j]/=pvval;
	
		for(i=0;i<dim;i++)
			if(i!=k)
			{
				int in1=i*n1;
				pvval=m[in1+k];
				for(j=k;j<n1;j++)
					m[in1+j]-=m[kn1+j]*pvval;
			}
	}

	for(i=0;i<dim;i++)
		x[i]=m[i*n1+dim];

	delete [] m;
	return true;
}

// Minimize overdetermined system Ax==b using least squares
template <class R>
bool
GJsolveLeastSq(
			const int nc, const int nr,
			const R *a,	// a[nr][nc] (a[row][col])
			const R *b,	// b[nr]
			R *x,		// returned solution x[nc]
			R *sa = 0, R *sb = 0) // optional returned square eq.
{
	int	i, j, k;
	if (nr == nc) {
		if (sa) for (i = nr*nc; i--; ) sa[i] = a[i];
		if (sb) for (i = nc; i--; ) sb[i] = b[i];
		
		return GJsolveLinearEq(nc, a, b, x);
	}
	if (nr < nc)
		return false;
	if (nc <= 0)
		return false;

	bool	our_sa = !sa;
	bool	our_sb = !sb;

	if (our_sa) sa = new R[nc*nc];
	if (our_sb) sb = new R[nc];

	for (i = nc; i--; )		// for each A^T row
		for (j = nc; j--; ) {		// for each A column
			sa[i*nc+j] = 0;
			for (k = nr; k--; )		// accum. dot prod.
				sa[i*nc+j] += a[k*nc+i] * a[k*nc+j];
		}
	for (i = nc; i--; ) {		// for each A^T row
		sb[i] = 0;
		for (k = nr; k--; )		// accum. dot prod. with b
			sb[i] += a[k*nc+i] * b[k];
	}
					// solve reduced linear system
	bool	ok = GJsolveLinearEq(nc, sa, sb, x);

	if (our_sa) delete [] sa;
	if (our_sb) delete [] sb;
	return ok;
}

#endif	// ! _GAUSSJORD_H_
