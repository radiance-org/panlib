/*
 *  system.c
 *  panlib
 *
 *  System-specific routines
 *
 *  Created by gward on Fri Jan 25 2002.
 *  Modified by plonghurst on Thursday Feb 15 2007.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include "system.h"

#if defined(_WIN32) || defined(_WIN64)
#include <sys/utime.h>
#endif

				/* Check if the given path is a directory */
int
isdirectory(const char *path)
{
	struct stat	sbuf;
				/* stat file */
	if (path == NULL || stat(path, &sbuf) < 0)
		return -1;

	return ((sbuf.st_mode & S_IFMT) == S_IFDIR);
}
	
				/* Get file modified time (UNIX seconds) */
time_t
getfiletime(const char *path)
{
	struct stat	sbuf;
				/* stat file */
	if (path == NULL || stat(path, &sbuf) < 0)
		return 0;
	
	return sbuf.st_mtime;
}

				/* Set file modified time (UNIX seconds) */
int
setfiletime(const char *path, time_t mtime)
{
/* pl: utimes() ~= utime() in Windows */
#if defined(_WIN32) || defined(_WIN64)
	struct _utimbuf tval[2];
	
	tval[0].modtime = time(0);		/* set access time to now */
	tval[0].actime = 0;
	if (!mtime)
		mtime = tval[0].modtime; /* special case 0==now */
	else if (mtime > tval[0].modtime)
		return -1;		/* modified time in future?? */
	tval[1].modtime = mtime;
	tval[1].actime = 0;

	return _utime(path, tval);
#else
	struct timeval  tval[2];

	tval[0].tv_sec = time(0);       /* set access time to now */
	tval[0].tv_usec = 0;
	if (!mtime)
		mtime = tval[0].tv_sec; /* special case 0==now */
	else if (mtime > tval[0].tv_sec)
		return -1;		/* modified time in future?? */
	tval[1].tv_sec = mtime;
	tval[1].tv_usec = 0;

	return utimes(path, tval);
#endif
}

				/* Convert UNIX seconds to "YYYY:MM:DD HH:MM:SS" */
char *
cvtdate(char dt[20], time_t tval)
{
	struct tm	*tm;

	if (tval == 0)
		return NULL;
	tm = localtime(&tval);
	if (tm == NULL)
		return NULL;
	sprintf(dt, "%04d:%02d:%02d %02d:%02d:%02d",
			tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday,
			tm->tm_hour, tm->tm_min, tm->tm_sec);
	return dt;
}


#if defined(_WIN32) || defined(_WIN64)

int
gettimeofday(struct timeval *tp, void *dummy)
{
    FILETIME        ft;
    LARGE_INTEGER   li;
    __int64         t;

	SYSTEMTIME		st;
	FILETIME		ft2;
	LARGE_INTEGER   li2;
	__int64			t2;

	st.wYear = 1970;
	st.wHour = 0;
	st.wMinute = 0;
	st.wSecond = 0;
	st.wMilliseconds = 1;

	SystemTimeToFileTime(&st, &ft2);
	li2.LowPart = ft2.dwLowDateTime;
	li2.HighPart = ft2.dwHighDateTime;
	t2 = li2.QuadPart;

    GetSystemTimeAsFileTime(&ft);
    li.LowPart  = ft.dwLowDateTime;
    li.HighPart = ft.dwHighDateTime;
    t  = li.QuadPart;      
    t -= t2; // From 1970
    t /= 10; // In microseconds
    tp->tv_sec  = (long)(t / 1000000);
    tp->tv_usec = (long)(t % 1000000);
    return 0;
}

#endif
