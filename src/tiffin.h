/*
 * Header file for reading TIFF data.
 *
 * Requires <iostream.h>
 */

#ifndef _TIFFIN_H_
#define _TIFFIN_H_

// TIFF data types
enum TIFFdatatype {
	TIFFany=0,
	TIFFunsigned_byte=1,
	TIFFascii_string=2,
	TIFFunsigned_short=3,
	TIFFunsigned_long=4,
	TIFFunsigned_rational=5,
	TIFFsigned_byte=6,
	TIFFundefined=7,
	TIFFsigned_short=8,
	TIFFsigned_long=9,
	TIFFsigned_rational=10,
	TIFFsingle_float=11,
	TIFFdouble_float=12
};

// TIFF data sizes
extern const short TIFFbytes_comp[];

// TIFF file position state
enum TIFFstate {
	TIFFgrabbed,
	TIFFunknown,
	TIFFbeforeMagic,
	TIFFbeforeTag,
	TIFFafterTag
};

#ifndef TIFFSTRLEN
#define TIFFSTRLEN		512		// default string buffer length
#endif

// TIFF input class
class TIFFin {
private:
	istream		*instr;			// pointer to input stream
	TIFFstate	state;			// current file position state
	bool		check();		// check TIFF header & set byte order
	unsigned long	tifflen;		// total file length in bytes
	bool		littlendian;		// file byte order Intel?
	bool		seek(off_t p) {		// go to a specific file postion
				if (p > tifflen) return false;
				instr->clear(); instr->seekg(p);
				return instr->good();
			}
	unsigned long	getWord(int);		// get word, swapping bytes as needed
	unsigned long	firstpos;		// offset of first tag in current IFD
	int		ntags;			// number of tags in this IFD
	int		tagnum;			// next tag to read
public:
	static const bool
			nativeLE;		// native little endian?
			TIFFin()
				{ instr=NULL; state=TIFFunknown; }
			TIFFin(const char *fname)
				{ instr=NULL; OpenFile(fname); }
			TIFFin(const char *buf, int len)
				{ instr=NULL; UseBuffer(buf, len); }
			~TIFFin()
				{ Close(); }
	void		Close()
				{ if (state != TIFFgrabbed) delete instr;
					instr=NULL; }
	bool		Ready() const
				{ return (instr != NULL & state != TIFFgrabbed); }
	bool		OpenFile(const char *fname);
	bool		UseBuffer(const char *buf, int len);
	bool		SetIFD(unsigned int ifd);
	bool		SeekIFD(unsigned long pos);
	istream *	GrabStream(unsigned long pos=0);
	void		ReleaseStream(istream *istr=NULL);
	int		NextTag();
	int		FirstTag()
				{ tagnum=0; state=TIFFunknown; return NextTag(); }
	bool		FindTag(int);
	int		GetData(void *dp, TIFFdatatype dt=TIFFany,
					int maxcomp=1, int elemsz=0);
	bool		GetTagValue(int tag, void *dp,
					TIFFdatatype dt=TIFFany);
	long		GetInteger(int tag=0);
	double		GetReal(int tag=0);
};

bool	cleanString(char *str);		// generic utility for cleaning up text string

#endif	// ! _TIFFIN_H_