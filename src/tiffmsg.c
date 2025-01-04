/*
 * Implementation of TIFF error and warning message buffering.
 */

#include "tiffio.h"
#include "tiffmsg.h"

char		TiffMessageBuffer[256];	/* most recent error or warning */

/* Save warning message to TiffMessageBuffer */
void
BufferTiffWarning(const char *module, const char *fmt, va_list ap)
{
	char	*bp = TiffMessageBuffer;
	if (module != NULL) {
		sprintf(bp, "%s: ", module);
		while (*bp) bp++;
	}
	sprintf(bp, "Warning, ");
	while (*bp) bp++;
	vsprintf(bp, fmt, ap);
}

/* Save error message to TiffMessageBuffer */
void
BufferTiffError(const char *module, const char *fmt, va_list ap)
{
	char	*bp = TiffMessageBuffer;
	if (module != NULL) {
		sprintf(bp, "%s: ", module);
		while (*bp) bp++;
	}
	vsprintf(bp, fmt, ap);
}
