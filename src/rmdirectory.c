/*
 *  rmdirectory.c
 *
 *  Remove a directory and all its subdirectories and files
 *
 *  Created by Greg Ward on Wed Apr 28 2004.
 *  Copyright (c) 2004 Anyhere Software. All rights reserved.
 *
 */

#include <stdlib.h>

#if !defined(_WIN32) && !defined(_WIN64)

#include <errno.h>
#include <string.h>
#include "system.h"
#include <dirent.h>
#include <sys/param.h>

#ifndef _D_ALLOC_NAMLEN
#define _D_ALLOC_NAMLEN(d)	((d)->d_namlen + 1)
#endif

				/* recursive directory removal (private) */
static int
_rmdirectory(char *pbuf)
{
	int		status = 0;
	size_t		dblen = 4*sizeof(struct dirent);
	char		*dbuf = (char *)malloc(dblen);
	int		dbend = 0;
	int		pblen;
	DIR		*dirp;
	struct dirent   *dp;
					/* read & memorize directory */
	dirp = opendir(pbuf);
	if (dirp == NULL)
		return -1;
	pblen = strlen(pbuf);
	pbuf[pblen++] = '/';
	while ((dp = readdir(dirp)) != NULL) {
		if (dp->d_name[0] == '.' && (!dp->d_name[1] ||
				(dp->d_name[1] == '.' && !dp->d_name[2])))
			continue;
		if (pblen + _D_ALLOC_NAMLEN(dp) > MAXPATHLEN) {
			errno = ENAMETOOLONG;
			status = -1;
			continue;
		}
		if (dbend + dp->d_reclen > dblen)
			dbuf = (char *)realloc(dbuf, dblen *= 2);
		if (dbuf == NULL) {
			errno = ENOMEM;
			closedir(dirp);
			return -1;
		}
		memcpy((dbuf + dbend), dp, dp->d_reclen);
		dbend += dp->d_reclen;
	}
	if (closedir(dirp) < 0)		/* close directory */
		status = -1;
					/* remove files & subdirectories */
	for (dp = (struct dirent *)dbuf; (char *)dp < dbuf+dbend;
			dp = (struct dirent *)((char *)dp + dp->d_reclen)) {
		strcpy(pbuf+pblen, dp->d_name);
		if ((dp->d_type == DT_DIR ? _rmdirectory(pbuf)
					: unlink(pbuf)) < 0)
			status = -1;
	}
	pbuf[--pblen] = '\0';
	free(dbuf);			/* release memory */
	if (status < 0)
		return -1;
	return rmdir(pbuf);		/* remove empty directory */
}

				/* Recursively remove a directory + contents */
int
removedirectory(const char *path)
{
	char		pbuf[MAXPATHLEN];

	if (!path || !*path)
		goto badarg;
	if (path[0] == '.') {
		if (!path[1])
			goto badarg;
		if (path[1] == '.' && (!path[2] || path[2] == '/'))
			goto badarg;
	}
	if (strlen(path) >= MAXPATHLEN) {
		errno = ENAMETOOLONG;
		return -1;
	}
	return _rmdirectory(strcpy(pbuf, path));
badarg:
	errno = EINVAL;
	return -1;
}

#else	/* Windows version */

int
removedirectory(const char *path)
{
	char	command[2048];

	strcpy(command, "rmdir /s /q ");
	strcat(command, path);

	return (system(command) == 0) ? 0 : -1;
}

#endif
