/*
 * Header for JPEG stream and error handling helper functions
 */

#ifndef _JSTREAMSRC_H_
#define	_JSTREAMSRC_H_

#include "jpeglib.h"
#include "jerror.h"
#include <iostream>
using namespace std;
#include <setjmp.h>

struct jumper_struct {
	jmp_buf	b;
	int	w;
};

extern void	jpeg_stream_src(j_decompress_ptr cinfo, istream *instream);

extern void	jpeg_stream_dest(j_compress_ptr cinfo, ostream *outstream);

extern void	jpeg_error_jump(j_common_ptr jp);

extern void	jpeg_emit_message(j_common_ptr jp, int msg_level);

extern void	jpeg_error_output(j_common_ptr jp);

#define JPEG_WARNTHRESH		12		// warning threshold for abort

// Errors that occur during JPEG loading are stored in the following buffer
extern char	jpeg_error_buffer[JMSG_LENGTH_MAX];

extern void	jpeg_sf_idct_float(j_decompress_ptr cinfo, jpeg_component_info *compptr,
		 JCOEFPTR coef_block, JSAMPARRAY output_buf, JDIMENSION output_col);

extern bool	jpeg_sf_idct(j_decompress_ptr cinfo);

#endif	// ! _JSTREAMSRC_H_
