/*
 * Include after "tiffio.h"
 *
 * TIFF error and warning handlers that save messages to buffer
 * rather than writing to stderr.  Enable by calling:
 *
 *	TIFFSetErrorHandler(BufferTiffError);
 *	TIFFSetWarningHandler(BufferTiffWarning);
 */

#ifndef _TIFFMSG_H_
#define _TIFFMSG_H_

#ifdef __cplusplus
extern "C" {
#endif

extern void	BufferTiffWarning(const char *module, const char *fmt, va_list ap);
extern void	BufferTiffError(const char *module, const char *fmt, va_list ap);

extern char	TiffMessageBuffer[];	/* most recent TIFF error or warning */

#ifdef __cplusplus
}
#endif

#endif	/* _TIFFMSG_H_ */
