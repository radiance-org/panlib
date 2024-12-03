/*
 * TIFF reader implementation
 */

#include <iostream>
#include <fstream>
using namespace std;
#include <string>
#include <sstream>
#include <ctype.h>
#include <string.h>
#include "system.h"
#include "tiffin.h"
#include "tiff.h"			/* needed for int32, etc. */

using namespace std;

const short TIFFbytes_comp[] = {0, 1, 1, 2, 4, 8, 1, 1, 2, 4, 8, 4, 8};

// Check machine byte order
static bool
IsBigendian()
{
	union { int32 i; char c[4]; } u;
	u.i = 1;
	return (u.c[0] == 0);
}

const bool	TIFFin::nativeLE = !IsBigendian();

// Return the nominal data size if the type matches one of the float types
static inline size_t
getRealSize(TIFFdatatype dt)
{
	switch (dt) {
	case TIFFsigned_rational:
	case TIFFunsigned_rational:
	case TIFFsingle_float:
		return sizeof(float);
	case TIFFdouble_float:
		return sizeof(double);
	}
	return 0;
}

// Remove leading and trailing white space from string
bool
cleanString(char *str)
{
	char *firstc = str;
	while (isspace(*firstc)) firstc++;		// find first non-white
	char *lastc;
	for (lastc = firstc; *lastc; lastc++)
		if ((*lastc < ' ') | (*lastc > '~'))	// stomp non-printing characters
			*lastc = '?';
	while (--lastc >= firstc && isspace(*lastc)) ;	// find last non-white
	if (lastc < firstc) {
		*str = '\0';				// empty string
		return false;
	}
	*++lastc = '\0';				// eliminate trailing white
	if (firstc > str)				// eliminate leading white
		while (firstc <= lastc)
			*str++ = *firstc++;
	return true;					// finito
}

// Check TIFF magic number and get size and byte order
bool TIFFin::check()
{
	if (instr == NULL)
		return false;
	if (state != TIFFbeforeMagic)		// position at beginning
		instr->seekg(0);
	switch (instr->get()) {			// get byte order
	case 'I':				// Intel
		if (instr->get() != 'I')
			goto fail;
		littlendian = true;
		break;
	case 'M':				// Motorola
		if (instr->get() != 'M')
			goto fail;
		littlendian = false;
		break;
	default:				// EOF or bad magic
		goto fail;
	}
	if (getWord(2) != 0x2a)			// 0x2a is tag mark
		goto fail;
	instr->seekg(0, ios::end);		// find file length
	tifflen = instr->tellg();
	state = TIFFunknown;
	firstpos = 0L;
	return true;
fail:
	Close();				// close input on failure
	return false;
}

// Read in a word from TIFF file, swapping bytes as necessary
unsigned long TIFFin::getWord(int wlen)
{
	unsigned char	buf[4];
	unsigned long	val = 0;
	unsigned char	*dp;
	instr->read((char *)buf, wlen);
	if (littlendian)
		for (dp = buf+wlen; dp > buf; ) {
			val <<= 8;
			val |= *--dp;
		}
	else
		for (dp = buf; dp < buf+wlen; ) {
			val <<= 8;
			val |= *dp++;
		}
	return val;
}

// Open a TIFF file for reading
bool TIFFin::OpenFile(const char *fname)
{
	Close();				// close old file
	if (fname == NULL)
		return false;
	instr = new ifstream(fname, ios::in|ios::binary);
	state = TIFFbeforeMagic;
	return check();
}

// Register TIFF buffer for reading
bool TIFFin::UseBuffer(const char *buf, int len)
{
	Close();				// close old file
	if ((buf == NULL) | (len <= 8))
		return false;
	instr = new istringstream(string(buf, len));
	state = TIFFbeforeMagic;
	return check();
}

// Set current IFD
bool TIFFin::SetIFD(unsigned int ifd)
{
	if (!Ready())
		return false;
	unsigned long	thisIFD;
	seek(4L);				// start at first tag mark
	state = TIFFunknown;
	for ( ; ; ifd--) {
		thisIFD = getWord(4);		// get offset to IFD
		if (thisIFD == 0L || !seek(thisIFD))
			return false;			// requested IFD doesn't exist
		if (!ifd)			// got it!
			break;
		if (!seek(thisIFD + 12*getWord(2) + 2))	// skip to next IFD pointer
			return false;
	}
	ntags = getWord(2);			// initialize IFD
	firstpos = thisIFD + 2;
	tagnum = 0;
	state = TIFFbeforeTag;
	return true;
}

// Seek to new IFD at user-specified location
bool TIFFin::SeekIFD(unsigned long pos)
{
	if (!Ready())
		return false;
	if (pos == 0L || !seek(pos))
		return false;
	ntags = getWord(2);
	firstpos = pos + 2;
	tagnum = 0;
	state = TIFFbeforeTag;
	return true;
}

// Seek to position and return stream pointer on success
istream * TIFFin::GrabStream(unsigned long pos)
{
	if (instr == NULL)
		return NULL;
	if (!seek(pos))
		return NULL;
	state = TIFFgrabbed;
	return instr;
}

// Release input stream after grab
void TIFFin::ReleaseStream(istream *istr)
{
	if (istr != instr)
		delete istr;			// we closed while grabbed
	else
		state = TIFFunknown;
}

// Read the next tag in this IFD, returning 0 at end
int TIFFin::NextTag()
{
	if (!Ready() || tagnum >= ntags)
		return -1;				// tags exhausted
	if (state == TIFFafterTag) {		// didn't use last tag -- skip it
		instr->ignore(10);
	} else if (state != TIFFbeforeTag) {	// need to reposition pointer
		if (firstpos == 0L)
			return -1;			// no IFD set!
		if (!seek(firstpos + 12*tagnum))
			return -1;
	}
	state = TIFFafterTag;			// prepare to return next tag
	tagnum++;
	return getWord(2);			// read it and return it
}

// Find the given tag in the current IFD
bool TIFFin::FindTag(int tag)
{
	int	thistag;
	for (thistag = FirstTag(); thistag >= 0; thistag = NextTag())
		if (thistag == tag)
			return true;
	return false;
}

// Get single integer value (won't convert real)
long TIFFin::GetInteger(int tag)
{
	if (!Ready())
		return 0;
	if (tag != 0 && !FindTag(tag))		// find tag if specified
		return 0;
	if (state != TIFFafterTag)
		return 0;
	TIFFdatatype	dt = (TIFFdatatype)getWord(2);
	bool		fitsword = (getWord(4)*TIFFbytes_comp[dt] <= 4);
	int32		lval;
	state = TIFFunknown;			// in case of premature return
	if (!fitsword && !seek(getWord(4)))
		return 0;
	switch (dt) {
	case TIFFsigned_byte:
	case TIFFunsigned_byte:
		lval = instr->get();
		if (fitsword) instr->ignore(3);
		if (dt == TIFFsigned_byte && lval & 0x80)
			lval |= ~0xff;
		break;
	case TIFFsigned_short:
	case TIFFunsigned_short:
		lval = getWord(2);
		if (fitsword) instr->ignore(2);
		if (dt == TIFFsigned_short && lval & 0x8000)
			lval |= ~0xffff;
		break;
	case TIFFsigned_long:
	case TIFFunsigned_long:
		lval = getWord(4);
		break;
	default:
		return 0;
	}
	if (fitsword)
		state = TIFFbeforeTag;
	return lval;
}

// Get single real value (will convert integer)
double TIFFin::GetReal(int tag)
{
	if (!Ready())
		return 0;
	if (tag != 0 && !FindTag(tag))		// find tag if specified
		return 0;
	if (state != TIFFafterTag)
		return 0;
	TIFFdatatype	dt = (TIFFdatatype)getWord(2);
	bool		fitsword = (getWord(4)*TIFFbytes_comp[dt] <= 4);
	uint32		num, den;
	double		rval;
	state = TIFFunknown;			// in case of premature return
	if (!fitsword && !seek(getWord(4)))
		return 0;
	switch (dt) {
	case TIFFsigned_byte:
	case TIFFunsigned_byte:
		num = instr->get();
		if (fitsword) instr->ignore(3);
		if (dt == TIFFsigned_byte && num & 0x80)
			rval = (int32)(num |= ~0xff);
		else
			rval = num;
		break;
	case TIFFsigned_short:
	case TIFFunsigned_short:
		num = getWord(2);
		if (fitsword) instr->ignore(2);
		if (dt == TIFFsigned_short && num & 0x8000)
			rval = (int32)(num |= ~0xffff);
		else
			rval = num;
		break;
	case TIFFsigned_long:
	case TIFFunsigned_long:
		num = getWord(4);
		if (dt == TIFFunsigned_long && num & 0x80000000)
			rval = (int32)num;
		else
			rval = num;
		break;
	case TIFFsigned_rational:
	case TIFFunsigned_rational:
		num = getWord(4);
		den = getWord(4);
		if (den == 0)
			return 0;		// don't deal infinity
		if (dt == TIFFsigned_rational)
			rval = (int32)num / (double)den;
		else
			rval = (double)num / (double)den;
		break;
	case TIFFsingle_float:			// assume IEEE float
		if (sizeof(float) != 4)
			return 0;
		num = getWord(4);
		rval = *(float *)&num;
		break;
	case TIFFdouble_float:			// assume IEEE double
		if (sizeof(double) != 8)
			return 0;
		num = getWord(4);
		den = getWord(4);
		if (littlendian == nativeLE) {
			((uint32 *)&rval)[0] = num;
			((uint32 *)&rval)[1] = den;
		} else {
			((uint32 *)&rval)[0] = den;
			((uint32 *)&rval)[1] = num;
		}
		break;
	default:
		return 0;
	}
	if (fitsword)
		state = TIFFbeforeTag;
	return rval;
}

// Read TIFF tag data
int TIFFin::GetData(void *dp, TIFFdatatype dt, int maxcomp, int elemsz)
{
	if (!Ready() || state != TIFFafterTag)
		return 0;
	state = TIFFunknown;			// in case of premature return
	if (!elemsz)
		elemsz = getRealSize(dt);
						// get/check type
	TIFFdatatype	indt = (TIFFdatatype)getWord(2);
	if (dt == TIFFany || (getRealSize(dt) && getRealSize(indt)))
		dt = indt;
	else if (dt != indt)
		return 0;
	if (!elemsz)				// get/check element size
		elemsz = TIFFbytes_comp[dt];
	else if (elemsz != TIFFbytes_comp[dt] && elemsz != getRealSize(dt))
		return 0;
	if (dt == TIFFascii_string) maxcomp--;	// for nul terminator
	int	ncomp = getWord(4);
	bool	fitsword = (ncomp*TIFFbytes_comp[dt] <= 4);
	if (ncomp > maxcomp) ncomp = maxcomp;
	if (ncomp <= 0)
		return 0;
	if (!fitsword && !seek(getWord(4)))
		return 0;
	int i = ncomp;
	switch (dt) {				// load data appropriately
	case TIFFascii_string:
		((char *)dp)[i] = '\0';		// assure nul termination
	    // fall through...
	case TIFFsigned_byte:
	case TIFFunsigned_byte:
	case TIFFundefined:
		instr->read((char *)dp, i);
		break;
	case TIFFsigned_long:
	case TIFFunsigned_long: {
			uint32 *lp = (uint32 *)dp;
			while (i--) *lp++ = getWord(4);
		} break;
	case TIFFsigned_short:
	case TIFFunsigned_short: {
			unsigned short *sp = (unsigned short *)dp;
			while (i--) *sp++ = getWord(2);
		} break;
	case TIFFsigned_rational:
	case TIFFunsigned_rational:
		if (elemsz == sizeof(float)) {
			float *rp = (float *)dp;
			uint32 num, den;
			while (i--) {
				num = getWord(4); den = getWord(4);
				if (den == 0)
					*rp++ = 0;	// don't deal infinity
				else if (dt == TIFFsigned_rational)
					*rp++ = (signed long)num / (float)den;
				else
					*rp++ = (float)num / (float)den;
			}
		} else if (elemsz == sizeof(double)) {
			double *rp = (double *)dp;
			uint32 num, den;
			while (i--) {
				num = getWord(4); den = getWord(4);
				if (den == 0)
					*rp++ = .0;	// don't deal infinity
				else if (dt == TIFFsigned_rational)
					*rp++ = (signed long)num / (double)den;
				else
					*rp++ = (double)num / (double)den;
			}
		} else
			return 0;
		break;
	case TIFFsingle_float: {		// assume IEEE fp
			if (sizeof(float) != 4)
				return 0;
			uint32 *lp = (uint32 *)dp;
			while (i--) *lp++ = getWord(4);
		} break;
	case TIFFdouble_float: {		// assume IEEE fp
			if (sizeof(double) != 2*sizeof(uint32))
				return 0;
			uint32	*lp = (uint32 *)dp;
			const int	swap = (littlendian != nativeLE);
			while (i--) {
				lp[swap] = getWord(sizeof(uint32));
				lp[1-swap] = getWord(sizeof(uint32));
				lp += 2;
			}
		} break;
	case TIFFany:				// avoided at routine initialization
		return 0;
	}
	if (fitsword) {
		if ((i = 4 - ncomp*TIFFbytes_comp[dt]) > 0)
			instr->ignore(i);
		state = TIFFbeforeTag;
	}
	return ncomp;
}

// Read tag and get value at the same time
bool TIFFin::GetTagValue(int tag, void *dp, TIFFdatatype dt)
{
	if (!Ready())
		return false;
	if (!FindTag(tag))
		return false;				// tag not found
	int elemsz = getRealSize(dt);
	if (!elemsz)
		elemsz = TIFFbytes_comp[dt];
	int len = 1;
	if (dt == TIFFascii_string)
		len = TIFFSTRLEN;			// big assumption -- take care!
	return (GetData(dp, dt, len, elemsz) > 0);
} 
