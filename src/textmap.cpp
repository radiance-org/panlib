/*
 *  textmap.cpp
 *  panlib
 *
 *  Created by Gregory Ward on Fri November 19, 2021.
 *  Copyright (c) 2021 Anyhere Software. All rights reserved.
 *
 *  Class to render text onto a 2-D bitmap
 */

#include "rtio.h"
#include "pdraw.h"
#include "textmap.h"
#include "font.h"

#define	INTERSPC	30		// intercharacter spacing
#define	FWASPECT	0.6f		// glyph width/height

					// standard font file suffixes
static const char	fontSuffix[4][16] = {
				".fnt",
				"Italic.fnt",
				"Bold.fnt",
				"BoldItalic.fnt"
			};

// Free any allocated pointers
TextMap::~TextMap()
{
	for (int i = 4; i--; )
		if (fref[i]) freefont(fref[i]);
}

// Change font flags and return what is set
int
TextMap::ChangeFFlags(int flgs)
{
	fndx = flgs & (TFFplain|TFFitalic|TFFbold);
	if (!(fpresent & 1<<fndx))		// no such variant?
		return fndx = 0;
	if (fref[fndx])				// already loaded?
		return fndx;
						// else load variant
	char	fullName[sizeof(baseFont)+sizeof(fontSuffix[0])];
	strcpy(fullName, baseFont);
	strcat(fullName, fontSuffix[fndx]);
	if (!(fref[fndx] = getfont(fullName))) {
		fpresent &= ~(1<<fndx);		// problem with file
		fndx = 0;
	}
	return fndx;
}

// Update our font height (new family name optional)
bool
TextMap::SetFont(float fht, const char *fn)
{
	EndLine();

	if (fht > 1)				// set height if legal
		fheight = fht;

	if (fn && *fn) {
		if (!strcmp(fn, baseFont))	// font already set?
			return true;
		int	i;
		for (i = 4; i--; )
			if (fref[i]) freefont(fref[i]);
		Clear();
		const int	fnlen = strlen(fn);
		if (fnlen >= sizeof(baseFont))
			return false;
		char	fullName[sizeof(baseFont)+sizeof(fontSuffix[0])];
		strcpy(fullName, fn);		// load plain font
		strcpy(fullName+fnlen, fontSuffix[0]);
		if (!(fref[0] = getfont(fullName)))
			return false;
		strcpy(baseFont, fn);		// got base font at least
		fpresent = 1;
		for (i = 4; --i; ) {		// check for others in family
			strcpy(fullName+fnlen, fontSuffix[i]);
			fpresent |= (getpath(fullName, getrlibpath(), 0) != NULL) << i;
		}
	}
	return (fref[fndx] != NULL);
}

// Retain fonts when not in use
void
TextMap::KeepFonts(bool set2)
{
	retainfonts = set2;
}

// Get word width (in points)
float
TextMap::GetWordWidth(const char *word) const
{
	if (!fref[fndx] | (fheight <= 1) || !word || !*word)
		return 0;

	short *		sp = new short [strlen(word)+1];
	int		wlen = squeeztext(sp, (char *)word, fref[fndx], INTERSPC);
	delete [] sp;

	return fheight*(FWASPECT/256.f) * wlen;
}

// Render a word at the given map coordinates
bool
TextMap::RenderWord(const char *word, float x, float y)
{
	if (!Ready())
		return false;
	if ((x > map.Width()) | (y-fheight > map.Height()))
		return true;
					// get character spacing
	const float	hf = fheight/256.f;
	const float	wf = hf*FWASPECT;
	const int	len = strlen(word);
	short *		sp = new short [len+1];
	int		wlen = squeeztext(sp, (char *)word, fref[fndx], INTERSPC);

	if (x + wf*wlen <= 0) {
		delete [] sp;		// word's end is left of margin
		return true;
	}
	int		(*vxy)[2] = new int [fref[fndx]->maxgv][2];
	for (int i = 0; i < len; i++) {	// draw and fill each character glyph
		if ((x += wf*sp[i]) >= map.Width())
			break;		// walked past right margin
		GLYPH *		gp = fref[fndx]->fg[word[i]&0xff];
		if (!gp) continue;	// undefined character?
		if (x + wf*gp->right < 0)
			continue;	// letter left of left margin
		GORD *		gv = gvlist(gp);
		for (int j = 0; j < gp->nverts; j++, gv += 2) {
			vxy[j][0] = int(x + wf*gv[0]);
			vxy[j][1] = int(y - hf*gv[1]);
		}
		BfillPolygon(&map, vxy, gp->nverts);
			// this letter done; advance to next
	}
	delete [] vxy;
	delete [] sp;
	return true;
}
