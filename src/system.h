/*
 *  system.h
 *  panlib
 *
 *  System-related dependencies
 *
 *  Created by gward on Mon Apr 30 2001.
 *  Modified by plonghurst on Thursday Feb 15 2007.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#ifndef _SYSTEM_H_
#define _SYSTEM_H_
				/* Windows */
#if defined(_WIN32) || defined(_WIN64)

	#ifdef _MSC_VER
		#include <BaseTsd.h>
		typedef SSIZE_T	ssize_t;
		#define HAVE_BOOLEAN
	#else
		#define isnan(x)	_isnan(x)
		#define isinf(x)	(!_finite(x))
		#ifndef snprintf
				#define snprintf	_snprintf
		#endif
	#endif

	#include <io.h>
	#include <time.h>
	#include <Winsock2.h>
	#include <sys/types.h>
	#include <sys/stat.h>
	#include <float.h>

	#define _Ios_Fmtflags	ios_base::fmtflags
	#define DIRSEP		'/'
	#define DIRSEP2		'\\'
	#define IS_DIRSEP(c)	(((c)==DIRSEP) | ((c)==DIRSEP2))

	#ifndef SET_DEFAULT_BINARY
		#define SET_DEFAULT_BINARY()	(_fmode = _O_BINARY)
		#define SET_FILE_BINARY(fp)	_setmode(fileno(fp),_O_BINARY)
		#define SET_FD_BINARY(fd)	_setmode(fd,_O_BINARY)
	#endif

	/*#ifndef snprintf*/
	/*	#define snprintf	_snprintf*/
	/*#endif*/

	#define strncasecmp	_strnicmp
	#ifndef strcasecmp
		#define strcasecmp	_stricmp
		#define strncasecmp	_strnicmp
	#endif

	#ifndef gethostname
		#define	gethostname	unix_gethostname
	#endif

	#ifndef kill
		#define	kill		unix_kill
	#endif

#endif
				/* Unix */
#ifndef DIRSEP

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#define DIRSEP		'/'

#endif

#include <string.h>
				/* Generic */
#ifndef MEM_CACHE_SIZE
#define MEM_CACHE_SIZE		(64L*1024*1024)
#endif
#ifndef TMP_FNAME_DIR
#define TMP_FNAME_DIR		NULL
#endif
#ifndef TMP_FNAME_PREFIX
#define TMP_FNAME_PREFIX	"PAN_"
#endif
#ifndef IS_DIRSEP
#define IS_DIRSEP(c)	((c)==DIRSEP)
#endif
#ifndef SET_DEFAULT_BINARY
#define SET_DEFAULT_BINARY()
#define SET_FILE_BINARY(fp)
#define SET_FD_BINARY(fd)
#endif

#ifdef __cplusplus
extern "C" {
#endif
				/* Recursively remove a directory + contents */
extern int	removedirectory(const char *path);
				/* Check whether the given file is a directory */
extern int	isdirectory(const char *path);
				/* Get file modified time (UNIX seconds) */
extern time_t	getfiletime(const char *path);
				/* Set file modified time (UNIX seconds) */
extern int	setfiletime(const char *path, time_t mtime);
				/* Convert UNIX seconds to "YYYY:MM:DD HH:MM:SS" */
extern char *	cvtdate(char dt[20], time_t tval);
				/* Get current date */
#define		getcurrentdate(dt)	cvtdate(dt, time(0))
				/* Get file modified date */
#define		getfiledate(dt,path)	cvtdate(dt, getfiletime(path))

#if defined(_WIN32) || defined(_WIN64)
extern int	unix_gethostname(char *hostname, size_t name_len);
extern int	unix_kill(long pid, int sig);
extern int	gettimeofday(struct timeval *tp, void *);
#endif

#ifdef __cplusplus
}
#endif

#endif /* _SYSTEM_H_ */
