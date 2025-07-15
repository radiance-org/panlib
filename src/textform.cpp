/*
 *  textform.cpp
 *  panlib
 *
 *  Created by Gregory Ward on Wed Oct 30 2002.
 *  Copyright (c) 2002 Anyhere Software. All rights reserved.
 *
 *  Abstract text rendering class implementation.
 */

#include <stdio.h>
#include <string.h>
#include "textform.h"

// Allocate and assign a word string
void
TFWordList::SetWord(const char *w)
{
	static char	emt[1];

	if (word && *word)
		delete [] word;
	if (!w || !*w) {
		word = emt;
		width = 0;
		addspc = false;
		return;
	}
	int	len = strlen(w);
	word = strcpy(new char [len+1], w);
	width = len;
	addspc = false;
}

// Put a word at the current cursor position (space if !*word)
bool
TextForm::PutWord(const char *word, int flgs, bool spaceAfter)
{
	if (trgn.top >= trgn.bottom)
		return false;
	if (!word)
		return true;
	if (!*word & !spaceAfter)
		return true;
	const float	spc = spaceAfter*spaceScale*GetFontSize();
	TFWordList *	tfw = new TFWordList(word);
	bool		ok = true;
	if (*word) {					// non-empty word?
		if ((tfw->wflags = ChangeFFlags(flgs)) == TFFerror) {
			tfw->wflags = TFFplain;
			ok = false;
		}
		tfw->width = GetWordWidth(word);	// check for full line
		if (wordList && wordList->addspc &
				(lwcur + tfw->width > GetWidth()))
			if (!EndLine()) {		// new line
				delete tfw;
				return false;
			}
	}
	tfw->addspc = spaceAfter;			// add word
	tfw->next = wordList;
	wordList = tfw;
	lwcur += tfw->width;
	if (spaceAfter)
		lwcur += spc;
	return ok;
}

// Put out an integer value
bool
TextForm::PutWord(long iv, int flgs, bool spaceAfter)
{
	char	wrd[64];
	sprintf(wrd, "%ld", iv);
	return PutWord(wrd, flgs, spaceAfter);
}

// Put out a float value
bool
TextForm::PutWord(double fv, int flgs, bool spaceAfter)
{
	char	wrd[64];
	sprintf(wrd, "%g", fv);
	return PutWord(wrd, flgs, spaceAfter);
}

// Break text into words and render from the current position
bool
TextForm::PutText(const char *text, int flgs)
{
	int		maxw;
	char		wtmp[1024];
	char *		dp;
							// set max. word length
	maxw = int(GetWidth()/(spaceScale*GetFontSize()));
	if (maxw >= (int)sizeof(wtmp))
		maxw = sizeof(wtmp)-1;
							// break into words
	for (dp = wtmp; ; text++) {
		switch (*text) {
		case '\r':
			if (text[1] == '\n')
				continue;		// cr-lf
		case '\n':
		case '\f':
			*dp = '\0';
			if (wtmp[0] && !PutWord(wtmp, flgs, false))
				return false;
			if (!EndLine(true))
				return false;
			if (!text[1])
				break;			// all done
			continue;
		case ' ':
		case '\t':
		case '\0':
			*dp = '\0';
			if (!PutWord(wtmp, flgs, (*text != '\0')))
				return false;
			if (!*text)
				break;			// all done
			if (*text == '\t') {
				PutSpace(); PutSpace(); PutSpace();
			}
			dp = wtmp;
			continue;
		default:				// normal char
			if (dp - wtmp >= maxw) {
				*dp = '\0';
				if (!PutWord(wtmp, flgs, false))
					return false;
				dp = wtmp;
			}
			*dp++ = *text;
			continue;
		}
		break;
	}
	return true;
}

// Finish the current line, advancing unless empty or advEmpty
bool
TextForm::EndLine(bool advEmpty)
{
	float		spc = spaceScale*GetFontSize();
	TFWordList *	tfw;
							// elide trailing spaces
	while (wordList && wordList->width <= 0) {
		tfw = wordList;
		wordList = tfw->next;
		tfw->next = NULL;
		if (tfw->addspc)
			lwcur -= spc;
		delete tfw;
	}
	if (!wordList) {				// check for empty line
		lwcur = 0;
		if (advEmpty)
			trgn.top += GetLineHeight();
		return (trgn.top <= trgn.bottom);
	}
	if (wordList->addspc) {				// remove final space
		wordList->addspc = false;
		lwcur -= spc;
	}
	if (justify == TJFfull) {			// adjust spacing
		int	spccnt = 0;
		for (tfw = wordList; tfw; tfw = tfw->next)
			spccnt += tfw->addspc;
		if (spccnt) {
			float	adjspc = spc + (GetWidth() - lwcur)/spccnt;
			if (maxAdj < 1 || adjspc <= maxAdj*spc) {
				spc = adjspc;
				lwcur = GetWidth();
			}
		}
	}
	const float	baseline = trgn.top + GetFontSize();
	float		xcur = trgn.left;
	bool		ok = true;
	if (justify & TJFleft)				// write from end
		xcur += lwcur;
	else if (justify & TJFright)
		xcur += GetWidth();
	else
		xcur += 0.5f*(GetWidth() + lwcur);
	for (tfw = wordList; ok & (tfw != NULL); tfw = tfw->next) {
		if (tfw->addspc)
			xcur -= spc;
		if (tfw->width > 0) {
			xcur -= tfw->width;
			ok &= (ChangeFFlags(tfw->wflags) != TFFerror);
			ok &= RenderWord(tfw->word, xcur, baseline);
		}
	}
	delete wordList; wordList = NULL;		// advance line
	lwcur = 0;
	trgn.top += GetLineHeight();
	ok &= (trgn.top <= trgn.bottom);
	return ok;
}
