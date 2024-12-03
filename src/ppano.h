/*
 *  ppano.h
 *  Photosphere
 *
 *  Created by Greg Ward on 5/4/06.
 *  Copyright 2006 Anyhere Software. All rights reserved.
 *
 *  Include after pancine.h
 *
 */

#ifndef _PPANO_H_
#define	_PPANO_H_

// Class definition for building panorama from positioned input images
class PPanorama {
private:
	mutable ImgColorSpace
			cspc;			// best target color space
	DBFieldInfo	fieldInfo;		// our field information
	ImgRect		anchRect;		// corresponding anchor box
	int		anchor;			// open anchor record index
	int		linkID;			// link field ID
	int		FindImage(const DBRecord &irec) const;
	int		InsertImage(const DBRecord &rec);
	bool		GetImageSize(int siz[2], int i) const;
	bool		Linked(int i1, int i2, ABitMap *bm,
				int offset[2] = NULL) const;
	bool		AlignAnchor(int matchpt[4], const DBRecord &irec2,
				const ImgRect &anchr2) const;
	bool		RenderLinked(ImgStruct *pimg, const float isca,
				int i1, const ImgStruct *ipar, ABitMap *bm,
				ABitMap2 *cvg = NULL);
public:
	DBRecordList	irl;			// our panoramic image set
						// progress report function
	int		(*reportProgress)(const char *, int);
			PPanorama() {
				cspc.cstatic = -1;
				anchor = -1;
				reportProgress = NULL;
			}
			// Clear the current image set
	void		ClearPano() {
				cspc.cstatic = -1;
				anchor = -1;
				irl.Init();
				fieldInfo.Clear();
			}
			// Can we set the given anchor origin?
	bool		CheckOrig(const DBRecord &irec1, const ImgRect &anchr1) const;
			// Assign image with anchor origin (opening)
	bool		LinkOrig(const DBRecord &irec1, const ImgRect &anchr1);
			// Verify anchor destination and check for loops
	bool		CheckDest(const DBRecord &irec2, const ImgRect &anchr2) const;
			// Assign image with anchor destination (closing)
	bool		LinkDest(const DBRecord &irec2, const ImgRect &anchr2);
			// Check to make sure we're ready & get size
	const char *	CheckPano(ImgColorSpace *dcsp = NULL,
					int siz[2] = NULL,
					int org0[2] = NULL) const;
			// Get corrected stonits (returns -1.f if unset)
	float		GetStoNits() const {
				float	stf;
				if (irl.GetSize() > 0 &&
				  PDBgetField(irl.Get(0),PDBFstonits).Get(&stf))
					return stf;
				return -1.f;
			}
			// Render panorama from current image set
	bool		RenderPano(ImgStruct *pimg, ABitMap2 *coverage = NULL);
};

#endif	// ! _PPANO_H_
