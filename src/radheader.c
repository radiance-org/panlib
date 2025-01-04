/*
 *  radheader.c
 *  pan
 *
 *  Common routines for handling Radiance header i/o
 *
 *  Created by Greg Ward on 7/27/19.
 *  Copyright 2019 Anyhere Software. All rights reserved.
 *
 */

#include <stdlib.h>
#include <ctype.h>
#include "imgio.h"
#include "radheader.h"
#include "color.h"
#include "resolu.h"

const short	RHortab[8] = {	/* orientation conversion table */
	YMAJOR|YDECR,
	YMAJOR|YDECR|XDECR,
	YMAJOR|XDECR,
	YMAJOR,
	YDECR,
	XDECR|YDECR,
	XDECR,
	0
};

/* Extract horizontal view angle from view parameter string */
int
RHgetAngle(ImgInfo *info)
{
	char	*cp;
	if (!(info->flags & IIFview))
		return 0;
	cp = strstr(info->view, "-vt");
	if (cp != NULL) {
		if (cp[3] == 'l') {
			info->hvangle = 0;	/* parallel view */
			info->flags |= IIFhvangle;
			return 1;
		}
		if (cp[3] != 'v')
			return 0;		/* cannot handle these */
	}
	cp = strstr(info->view, "-vh ");
	if (cp == NULL)			/* should use default?? */
		return 0;
	info->hvangle = atof(cp+4);
	info->flags |= IIFhvangle;
	return 1;
}

/* Check if header line fits "VAR= value" format */
int
RHisParameter(const char *hl)
{
	if (!isalpha(*hl++))
		return 0;
	while (isalnum(*hl) | (*hl == '_'))
		hl++;
	if (*hl++ != '=')
		return 0;
	while (isspace(*hl))
		hl++;

	return (*hl != '\0');
}

/* Scan header line, modifying ImgInfo struct appropriately */
int
RHscanHeadline(char *hl, void *p)
{
	ImgInfo *       info = (ImgInfo *)p;
	int		i;

	if (isheadid(hl))		/* ignore header ID */
		return 0;
	if (isformat(hl))		/* presume caller knows format */
		return 0;
	if (isaspect(hl))		/* ditto */
		return 0;
	if (isprims(hl))		/* ditto */
		return 0;
	if (isexpos(hl))		/* ditto */
		return 0;
	if (isbigendian(hl) >= 0)	/* not a parameter */
		return 0;
	if (isncomp(hl))		/* ditto */
		return 0;
	if (!strncmp(hl, "NCOLS=", 6))	/* ditto */
		return 0;
	if (!strncmp(hl, "NROWS=", 6))	/* ditto */
		return 0;
	if (iswlsplit(hl))		/* spectral data hazard! */
		return 0;
					/* record version if Radiance */
	if (!strncmp(hl, "SOFTWARE= RADIANCE ", 19)) {
		info->source.stype = ISTrender;
		strcpy(info->source.make, "LBNL");
		strcpy(info->source.model, "RADIANCE");
		for (i = 19; hl[i] && !isdigit(hl[i]); i++)
			;
		strlcpy(info->source.vers, hl+i, sizeof(info->source.vers));
		for (i = 0; info->source.vers[i] &&
				!isspace(info->source.vers[i]); i++)
			;
		info->source.vers[i] = '\0';
		info->flags |= IIFsource;
		return 1;
	}
					/* add to Radiance view */
	if (!strncmp(hl, "VIEW= ", 6) && hl[6]) {
		if (info->flags & IIFview)
			i = strlen(info->view);
		else
			i = 0;
		strlcpy(info->view+i, hl+5, sizeof(info->view)-i);
		i += strlen(info->view+i) - 1;
		if (info->view[i] == '\n')
			info->view[i] = '\0';
		info->flags |= IIFview;
		RHgetAngle(info);	/* get horizontal angle */
		return 1;
	}
					/* save capture date */
	if (isdate(hl)) {
		hl += 8; while (isspace(*hl)) hl++;
		strncpy(info->capdate, hl, 19);
		info->capdate[19] = '\0';
		info->flags |= IIFcapdate;
		return 1;
	}
					/* save GMT */
	if (isgmt(hl)) {
		hl += 4; while (isspace(*hl)) hl++;
		strncpy(info->gmtdate, hl, 19);
		info->gmtdate[19] = '\0';
		info->flags |= IIFgmtdate;
		return 1;
	}
					/* save latitude/longitude */
	if (islatlon(hl)) {
		latlonval(info->latlong, hl);
		info->flags |= IIFlatlong;
		return 1;
	}
					/* save owner */
	if (!strncmp(hl, "OWNER= ", 7)) {
		strlcpy(info->owner, hl+7, sizeof(info->owner));
		i = strlen(info->owner) - 1;
		if (info->owner[i] == '\n')
			info->owner[i] = '\0';
		info->flags |= IIFowner;
		return 1;
	}
					/* other parameter? */
	if (RHisParameter(hl)) {
		i = strlen(info->params);
		if (i + strlen(hl) < sizeof(info->params)) {
			strcpy(info->params+i, hl);
			info->flags |= IIFparams;
		}
		return 1;
	}
					/* save misc. as comment */
	i = strlen(info->comments);
	if (i + strlen(hl) < sizeof(info->comments)) {
		strcpy(info->comments+i, hl);
		info->flags |= IIFcomments;
		return 1;
	}
	return 0;
}

/* Get header information and bundle for calling application */
int
RHgetInfo(FILE *fp, ImgInfo *info)
{
	if (!fp | !info)
		return -1;
					/* scan info. header from current */
	return getheader(fp, &RHscanHeadline, info);
}

/* Create Radiance information header */
void
RHstartHeader(const ImgInfo *info, FILE *pout)
{
	if (!info | !pout)
		return;
	newheader((char *)"RADIANCE", pout);
						/* write software version */
	if (info->flags & IIFsource) {
		if ((info->source.stype == ISTrender) |
				(info->source.stype == ISTeditor))
			fprintf(pout, "SOFTWARE= %s %s by %s\n",
					info->source.model,
					info->source.vers,
					info->source.make);
		else if (info->source.stype == ISTdigicam)
			fprintf(pout, "CAMERA= %s %s version %s\n",
					info->source.make,
					info->source.model,
					info->source.vers);
		else if ((info->source.stype == ISTfilm) |
				(info->source.stype == ISTflatbed))
			fprintf(pout, "SCANNER= %s %s version %s\n",
					info->source.make,
					info->source.model,
					info->source.vers);
	}
						/* comments and parameters */
	if (info->flags & IIFcomments && info->comments[0] != '\n' &&
			!strstr(info->comments, "\n\n"))
		fputs(info->comments, pout);
	if (info->flags & IIFparams)
		fputs(info->params, pout);
						/* write view parameters */
	if (info->flags & IIFview) {
		fputs("VIEW=", pout);
		if (info->view[0] != ' ')
			fputc(' ', pout);
		fputs(info->view, pout);
		if (info->view[strlen(info->view)-1] != '\n')
			fputc('\n', pout);
	}
						/* write owner */
	if (info->flags & IIFowner)
		fprintf(pout, "OWNER= %s\n", info->owner);
						/* write capture date */
	if (info->flags & IIFcapdate)
		fprintf(pout, "CAPDATE= %s\n", info->capdate);
	if (info->flags & IIFgmtdate)
		fprintf(pout, "GMT= %s\n", info->gmtdate);
						/* write latitude/longitude */
	if (info->flags & IIFlatlong)
		fputlatlon(info->latlong[0], info->latlong[1], pout);
}
