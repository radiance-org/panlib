/*
 *  textmap.h
 *  panlib
 *
 *  Created by Gregory Ward on Fri November 19, 2021.
 *  Copyright (c) 2021 Anyhere Software. All rights reserved.
 *
 *  Class to render text onto a 2-D bitmap
 */

#ifndef _TEXTMAP_H_
#define _TEXTMAP_H_

#include "textform.h"
#ifndef _ABITMAP_H_
#include "abitmap.h"
#endif

struct font;				// forward declaration

// Derived class using a polygonal font to render text into 2-D bitmaps
class TextMap : public TextForm {
	float		fheight;	// font height (adjustable)
	char		baseFont[64];	// base font name
	int		fpresent;	// flags of known font subfamilies
	struct font *	fref[4];	// loaded poly font reference(s)
	int		fndx;		// current font index
protected:
			// Clear class
	void		Clear() {
				baseFont[0] = '\0'; fndx = 0; fpresent = 0;
				fref[0] = fref[1] = fref[2] = fref[3] = 0;
			}
			// Change font flags and return what is set
	int		ChangeFFlags(int flgs);
			// Get word width (in points)
	float		GetWordWidth(const char *word) const;
			// Render a word at the given global coordinates
	bool		RenderWord(const char *word, float x, float y);
public:
	ABitMap2	map;		// output 2-D bitmap (canvas)
			TextMap() {
				Clear();
				fheight = 12.f;
			}
			TextMap(int w, int h, const char *fn=0, float fht=0) : map(w,h) {
				Clear();
				fheight = 12.f;
				SetFont(fht, fn);
				SetRect(0, 0, w, h);
			}
	virtual		~TextMap();
			// Erase/resize canvas
	void		Erase(int w=0, int h=0) {
				if (w <= 0) w = map.Width();
				if (h <= 0) h = map.Height() ? map.Height() : int(fheight+.999f);
				map.NewBitMap(w, h);
				SetRect(0, 0, w, h);
			}
			// Update our font height (new family name optional)
	bool		SetFont(float fht, const char *fn=0);
			// Retain fonts when not in use
	static void	KeepFonts(bool set2 = true);
			// Return current font height
	float		GetFontSize() const {
				return fheight;
			}
			// Are we ready to render some text?
	bool		Ready() const {
				return (fref[fndx] && fheight > 1 && map.Height() > 0);
			}
};

#endif	// !_TEXTMAP_H_
