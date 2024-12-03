/*
 *  textform.h
 *  panlib
 *
 *  Created by Gregory Ward on Wed Oct 30 2002.
 *  Copyright (c) 2002 Anyhere Software. All rights reserved.
 *
 *  Abstract text rendering class declarations.
 */
 
#ifndef _TEXTFORM_H_
#define _TEXTFORM_H_
					// flags for font variations
enum TextFormFlags {
	TFFplain	= 0x00,
	TFFitalic	= 0x01,
	TFFbold		= 0x02,
	TFFunderscore	= 0x04,
	TFFstrikeout	= 0x08,
	TFFsuperscript	= 0x10,
	TFFsubscript	= 0x20,
	TFFerror	= -1
};
					// flags for margin justification
enum TextJustFlags {
	TJFcenter	= 0x0,
	TJFleft		= 0x1,
	TJFright	= 0x2,
	TJFfull		= 0x3
};

// A renderable text area
struct TFRect {
	float		left, right, top, bottom;
};

// A list of words and widths
struct TFWordList {
	TFWordList *	next;			// next in list
	char *		word;			// allocated word
	int		wflags;			// associated font flags
	float		width;			// word width (points)
	bool		addspc;			// add space after word
			TFWordList() {
				static char	emt[1];
				next = 0; word = emt; width = 0; addspc = false;
			}
			TFWordList(const char *w, TFWordList *n = 0) {
				word = 0; next = n;
				SetWord(w);
			}
			~TFWordList() {
				delete next;
				if (word && *word) delete [] word;
			}
	void		SetWord(const char *w);
};

// Abstract base class for formatting text to a graphics device
class TextForm {
	float		lwcur;			// width of line in progress
	TFWordList *	wordList;		// our word list (stack order)
protected:
	TFRect		orgn;			// renderable text region
	TFRect		trgn;			// region remaining
	float		lineSpacing;		// current line spacing
			// Get remaining renderable height (in points)
	float		RemainingHeight() const {
				return trgn.bottom - trgn.top;
			}
			// Change font flags and return what is set
	virtual int	ChangeFFlags(int flgs) {
				return TFFplain;
			}
			// Get word width (in points)
	virtual float	GetWordWidth(const char *word) const = 0;
			// Render a word at the given global coordinates
	virtual bool	RenderWord(const char *word, float x, float y) = 0;
public:
	float		spaceScale;		// inter-word spacing multiplier
	float		lineScale;		// font spacing multiplier
	float		maxAdj;			// maximum space scale for TJFfull
	TextJustFlags	justify;		// left/right justification
			TextForm() {
				lwcur = 0; wordList = 0;
				trgn.left = trgn.right = trgn.top = trgn.bottom = 0;
				orgn = trgn;
				lineSpacing = 1;
				spaceScale = 0.33f;
				lineScale = 1.25f;
				maxAdj = 2.2f;
				justify = TJFleft;
			}
	virtual		~TextForm() {
				delete wordList;
			}
			// Get renderable rectangle
	const TFRect &	GetRect() const {
				return orgn;
			}
			// Set new renderable rectangle (call before output)
	void		SetRect(const TFRect &rct) {
				EndLine();
				trgn = orgn = rct;
			}
	void		SetRect(float x1, float y1, float x2, float y2) {
				TFRect	rct;
				if (x1 <= x2) { rct.left = x1; rct.right = x2; }
				else { rct.left = x2; rct.right = x1; }
				if (y1 <= y2) { rct.top = y1; rct.bottom = y2; }
				else { rct.top = y2; rct.bottom = y1; }
				SetRect(rct);
			}
			// Get renderable width (in points)
	float		GetWidth() const {
				return orgn.right - orgn.left;
			}
			// Get renderable height (in points)
	float		GetHeight() const {
				return orgn.bottom - orgn.top;
			}
			// Get the current number of lines
	int		CurrentLine() const {
				return int((trgn.top - orgn.top)/GetLineHeight()+0.5f);
			}
			// Get number of remaining lines
	int		RemainingLines() const {
				return int(RemainingHeight()/GetLineHeight());
			}
			// Get current line height
	float		GetLineHeight() const {
				return lineSpacing*lineScale*GetFontSize();
			}
			// Set new line height by changing spacing
	bool		SetLineHeight(float lh) {
				return SetLineSpacing(lh/(lineScale*GetFontSize()));
			}
			// Get current line spacing (relative to single)
	float		GetLineSpacing() const {
				return lineSpacing;
			}
			// Set new line spacing (relative to single)
	bool		SetLineSpacing(float mul = 1) {
				if (!EndLine()) return false;
				lineSpacing = mul;
				return true;
			}
			// Set left/right justification
	bool		SetJustify(TextJustFlags jf = TJFleft) {
				if (!EndLine()) return false;
				justify = jf;
				return true;
			}
			// Put a word out on the current line
	bool		PutWord(const char *word, int flgs = TFFplain,
					bool spaceAfter = true);
			// Put out an integer value
	bool		PutWord(long iv, int flgs = TFFplain,
					bool spaceAfter = true);
			// Put out a float value
	bool		PutWord(double fv, int flgs = TFFplain,
					bool spaceAfter = true);
			// Add an extra space
	bool		PutSpace() {
				return PutWord("");
			}
			// Put a required (fixed, non-breaking) space
	bool		PutRSpace(int flgs = TFFplain) {
				return PutWord(" ", flgs, false);
			}
			// Put out a single character
	bool		PutChar(char ch, int flgs = TFFplain) {
				if (ch == ' ') return PutSpace();
				if (ch == '\n') return EndLine(true);
				char	str[2];
				str[0] = ch; str[1] = '\0';
				return PutWord(str, flgs, false);
			}
			// Break text into words and render from the current position
	bool		PutText(const char *text, int flgs = TFFplain);
			// Finish the current line (and advance)
	bool		EndLine(bool advEmpty = false);
			// Get the current font height (points)
	virtual float	GetFontSize() const = 0;
};

inline TextForm &
operator<<(TextForm &tf, char ch)
{
	tf.PutChar(ch);
	return tf;
}

inline TextForm &
operator<<(TextForm &tf, const char *str)
{
	tf.PutText(str);
	return tf;
}

inline TextForm &
operator<<(TextForm &tf, long iv)
{
	tf.PutWord(iv, TFFplain, false);
	return tf;
}

inline TextForm &
operator<<(TextForm &tf, int iv)
{
	tf.PutWord((long)iv, TFFplain, false);
	return tf;
}

inline TextForm &
operator<<(TextForm &tf, double fv)
{
	tf.PutWord(fv, TFFplain, false);
	return tf;
}

#endif	// ! _TEXTFORM_H_
