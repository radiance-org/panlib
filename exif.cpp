/*
 *  exif.cpp
 *  panlib
 *
 *  Extraction routines for Exif headers
 *
 *  Created by gward on Tue May 29 2001.
 *  Copyright (c) 2001 Anyhere Software. All rights reserved.
 *
 */

#include <math.h>
#include <ctype.h>
#include <iostream>
using namespace std;
#include "rtio.h"
#include "color.h"
#include "tiffin.h"
#include "pimage.h"
#include "exif.h"
#include "tiff.h"			/* needed for int32, etc. */

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

// Read & interpret Exif color space
bool
GetExifCS(ImgColorSpace *ics, TIFFin *ti)
{
	if ((ics == NULL) | (ti == NULL))
		return false;
	if (!ti->SetIFD(0))			// main TIFF directory
		return false;
	const int	bits_sample = ti->GetInteger(0x0102);
	const int	samps_pixel = ti->GetInteger(0x0115);
						// start making assumptions
	switch (bits_sample) {
	case 8:
		if (samps_pixel == 3) {
			PcopyCS(ics, &ICS_sRGB);
			break;
		}
		if (samps_pixel == 1) {
			PcopyCS(ics, &ICS_Y8);
			return true;
		}
		return false;
	case 0:
		PcopyCS(ics, &ICS_sRGB);	// must be JPEG-compressed
		break;
	default:
		return false;			// punt other cases for now
	}
	bool	gotColorPrims = false;		// has primary chromaticities?
	if (samps_pixel != 1 && ti->FindTag(0x013f)) {
		if (ti->GetData(ics->chroma, TIFFsigned_rational, 6) != 6)
			return false;
		if (ti->FindTag(0x013e))	// and white point?
			ti->GetData(ics->chroma[3], TIFFsigned_rational, 2);
		int	ok = colorprimsOK((RGBPRIMP)ics->chroma);
		if (!ok)
			return false;
		if (ok < 0)
			ics->format = IPFxyz;
		gotColorPrims = true;
	}
	uint32		Offset;			// go to Exif portion
	if (!ti->GetTagValue(0x8769, &Offset, TIFFunsigned_long) ||
			!ti->SeekIFD(Offset))
		return false;
	float	gamval = ti->GetReal(0xa500);
	if (gamval > .5f)
		ics->gamma = gamval;
	if (gotColorPrims)			// good enough?
		return true;
						// else check other flags
	switch (ti->GetInteger(0xa001)) {	// check color space tag
	case 1:					// sRGB space for sure
		return true;
	case 2:					// *probably* Adobe RGB
		PcopyCS(ics, &ICS_AdobeRGB);
		return true;
	default:				// assume "uncalibrated"
		break;
	}
						// check interoperability
	if (!ti->GetTagValue(0xa005, &Offset, TIFFunsigned_long) ||
			!ti->SeekIFD(Offset))
		return false;
	char	io_string[TIFFSTRLEN];
	if (!ti->GetTagValue(0x1, io_string, TIFFascii_string))
		return false;
	if (!strcmp(io_string, "R03")) {
		PcopyCS(ics, &ICS_AdobeRGB);	// Adobe RGB identifier
		return true;
	}
	return false;				// else "uncalibrated"
}

// Check that DateTime value is reasonable
static bool
dateOK(const char *dt)
{
	int		yr, mo, da, hr, mn, sc;
	if (dt == NULL)
		return false;
	if (strlen(dt) != 19)
		return false;
	if (sscanf(dt, "%d:%d:%d %d:%d:%d", &yr, &mo, &da, &hr, &mn, &sc) != 6)
		return false;
	return !((yr < 1850) | (mo < 1) | (mo > 12) | (da < 1) | (da > 31) |
			(hr < 0) | (hr > 23) | (mn < 0) | (mn > 59) |
				(sc < 0) | (sc > 59));
}

// Parse and distribute image description to appropriate info field(s)
static void
parseDescription(ImgInfo *info, char *descr)
{
	if (!cleanString(descr))
		return;
	char *	cp;
	for (cp = descr; *cp; cp++)		// put back line feeds
		if (*cp == ';')
			*cp = '\n';
	if (cp[-1] != '\n') {
		*cp++ = '\n'; *cp = '\0';
	}
	char *	start = descr;
	bool	isParam = false;
	for (cp = descr; *cp; cp++) {		// get each line
		isParam |= (*cp == '=');
		if (*cp != '\n') continue;
		isParam &= isalpha(*start);
		if (isParam) {			// parameter?
			strncat(info->params, start, cp-start+1);
			info->flags |= IIFparams;
		} else if (cp > start) {
			strncat(info->comments, start, cp-start+1);
			info->flags |= IIFcomments;
		}
		start = cp+1;
		isParam = false;
	}
}

// Set camera properties from Exif header
bool
GetExifInfo(ImgInfo *info, TIFFin *ti)
{
#define make	(info->source.make)
#define model	(info->source.model)
#define rev	(info->source.vers)
#define artist	(info->owner)
	if ((info == NULL) | (ti == NULL))
		return false;
	if (!ti->SetIFD(0))			// get relevant IFD0 tags
		return false;
	int		tag;
	float		rval;
	unsigned short	sval;
	char		date[20];
	char		descr[4096];
	uint32		ExifOffset = 0;
	uint32		GPSOffset = 0;
	uint32		horiz_size = 0;
	float		horiz_dens = -1;
	float		horiz_unit_mm = 25.4f;
	info->source.stype = ISTdigicam;	// Exif assumption
	make[0] = model[0] = '\0'; rev[0] = '\0';
	artist[0] = '\0'; descr[0] = '\0';
	for (tag = ti->FirstTag(); tag >= 0; tag = ti->NextTag())
		switch (tag) {
		case 0x0100:			// Image Width
			horiz_size = ti->GetInteger();
			break;
		case 0x0112:			// Image orientation
			info->orientation = (ImgOrientation)ti->GetInteger();
			if (info->orientation != IOtopleft)
				info->flags |= IIForientation;
			break;
		case 0x011a:			// Horizontal Density
			horiz_dens = ti->GetReal();
			break;
		case 0x0128:			// Density Units
			switch (ti->GetInteger()) {
			case 2:			// inches
				horiz_unit_mm = 25.4f;
				break;
			case 3:			// centimeters
				horiz_unit_mm = 10.f;
			}
			break;
		case 0x8769:			// Exif Offset
			ti->GetData(&ExifOffset, TIFFunsigned_long);
			break;
		case 0x8825:			// GPS Offset
			ti->GetData(&GPSOffset, TIFFunsigned_long);
			break;
		case 0x882a:			// Time Zone Offset
			// signed short value -- meaning??
			break;
		case 0x010e:			// Image Description
			if (ti->GetData(descr, TIFFascii_string, sizeof(descr)))
				parseDescription(info, descr);
			break;
		case 0x010f:			// Camera Make
			ti->GetData(make, TIFFascii_string, sizeof(make));
			break;
		case 0x0110:			// Camera Model
			ti->GetData(model, TIFFascii_string, sizeof(model));
			break;
		case 0x0131:			// Firmware Revision
			ti->GetData(rev, TIFFascii_string, sizeof(rev));
			break;
		case 0x013b:			// Artist
			if (ti->GetData(artist, TIFFascii_string, sizeof(artist))
					&& cleanString(artist))
				info->flags |= IIFowner;
			break;
		case 0x923f:			// Sample-to-nits conversion
			info->stonits = ti->GetReal();
			if (info->stonits > 0)
				info->flags |= IIFstonits;
			break;
		case 0x0132:			// DateTime modified
			if (ti->GetData(date, TIFFascii_string, sizeof(date))
					&& dateOK(date)) {
				strcpy(info->capdate, date);
				info->flags |= IIFcapdate;
			}
			break;
		}
	if (!cleanString(rev))
		strcpy(rev, "v.0");
	if (cleanString(make) && cleanString(model))
		info->flags |= IIFsource;
	if (horiz_dens > 0) {
		info->hdensity = horiz_dens * 25.4 / horiz_unit_mm;
		info->flags |= IIFhdensity;
	}
						// get relevant Exif tags
	if (ExifOffset == 0 || !ti->SeekIFD(ExifOffset))
		return false;
	int		focal_35mm = -1;
	horiz_dens = -1.f;			// now means image plane dens.
	horiz_unit_mm = 25.4f;
	for (tag = ti->FirstTag(); tag >= 0; tag = ti->NextTag())
		switch (tag) {
		case 0x829a:			// ExposureTime
			if (!ti->GetData(&rval, TIFFunsigned_rational))
				break;
			if ((rval < 1e-6f) | (rval > 3600.f))
				break;
			info->exptime = rval;
			info->flags |= IIFexptime;
			break;
		case 0x829d:			// FNumber
			if (!ti->GetData(&rval, TIFFunsigned_rational))
				break;
			if ((rval < 1.f) | (rval > 128))
				break;
			info->aperture = rval;
			info->flags |= IIFaperture;
			break;
		case 0x8827:			// ISO SpeedRatings
			if (!ti->GetData(&sval, TIFFunsigned_short))
				break;
			if (sval < 1)
				break;
			info->asa = sval;
			info->flags |= IIFasa;
			break;
		case 0x0132:			// DateTime modified
			if (info->flags & IIFcapdate)
				break;
			// fall through
		case 0x9004:			// Date TimeDigitized
		case 0x9003:			// Date TimeOriginal
			if (ti->GetData(date, TIFFascii_string, sizeof(date))
					&& dateOK(date)) {
				strcpy(info->capdate, date);
				info->flags |= IIFcapdate;
			}
			break;
		case 0x9201:			// ShutterSpeedValue
			if (info->flags & IIFexptime)
				break;
			if (!ti->GetData(&rval, TIFFsigned_rational))
				break;
			if ((rval < -11.8) | (rval > 15.6))
				break;
			info->exptime = pow(.5f, rval);
			info->flags |= IIFexptime;
			break;
		case 0x9202:			// ApertureValue
			if (info->flags & IIFaperture)
				break;
			if (!ti->GetData(&rval, TIFFunsigned_rational))
				break;
			if ((rval == 0) | (rval > 14))
				break;
			info->aperture = pow(2.f, .5f*rval);
			info->flags |= IIFaperture;
			break;
		case 0x9206:			// SubjectDistance
			rval = ti->GetReal();
			if (rval <= 0)
				break;
			info->focus = rval;
			info->flags |= IIFfocus;
			break;
		case 0x9208:			// LightSource
			if (!ti->GetData(&sval, TIFFunsigned_short))
				break;
			info->whitebal = (ImgWhiteBal)sval;
			info->flags |= IIFwhitebal;
			break;
		case 0x9209:			// Flash
			if (!ti->GetData(&sval, TIFFunsigned_short))
				break;
			info->flash = (ImgFlashMode)sval;
			info->flags |= IIFflash;
			break;
#if 0					// seems to be unused
		case 0x9211:			// Image Number
			ti->GetData(&lval, TIFFunsigned_long);
			printf("ImageNumber: %ld\n", lval);
			break;
#endif
		case 0x920a:			// Focal Length (mm)
			info->flength = ti->GetReal();
			info->flags |= IIFflength*(info->flength > 1.f);
			break;
		case 0x9286: {			// User comment
			char	cdata[2048];
			int	nc = ti->GetData(cdata, TIFFany, sizeof(cdata));
			if (nc > 8 && !memcmp(cdata, "ASCII\0\0", 8)) {
				int	i = nc;
				while (i-- > 8)
					if (isspace(cdata[i]))
						--nc;
					else
						break;
				if (i >= 8) {
					cdata[nc++] = '\n';
					cdata[nc] = '\0';
					for (i = 8; isspace(cdata[i]); i++)
						;
					strlcat(info->comments, cdata+i,
							sizeof(info->comments));
					info->flags |= IIFcomments;
				}
			}
			} break;
		case 0xa405:			// 35mm equiv. Focal Length
			focal_35mm = ti->GetInteger();
			break;
		case 0xa002:			// Image Width
			horiz_size = ti->GetInteger();
			break;
		case 0xa20e:			// Focal Plane X Resolution
			horiz_dens = ti->GetReal();
			break;
		case 0xa210:			// Focal Plane Resolution Unit
			if (!ti->GetData(&sval, TIFFunsigned_short))
				break;
#if 1					// XXX which is correct??
			switch (sval) {
			case 2:			// inches
				horiz_unit_mm = 25.4f; break;
			case 3:			// centimeters
				horiz_unit_mm = 10.f; break;
			}
#else
			switch (sval) {
			case 1:			// inches
				horiz_unit_mm = 25.4f; break;
			case 2:			// meters
				horiz_unit_mm = 1000.f; break;
			case 3:			// centimeters
				horiz_unit_mm = 10.f; break;
			case 4:			// millimeters
				horiz_unit_mm = 1.f; break;
			case 5:			// microns
				horiz_unit_mm = .001f; break;
			}
#endif
			break;
		case 0xa215:			// Exposure Index (ASA)
			if (info->flags & IIFasa)
				break;
			if (!ti->GetData(&rval, TIFFunsigned_rational))
				break;
			if ((rval < 1) | (rval > 100000))
				break;
			info->asa = rval;
			info->flags |= IIFasa;
			break;
#if 0
		default:
			fprintf(stderr, "Exif tag 0x%x ignored\n", tag);
			break;
#endif
		}
						// compute horizontal angle
	if (horiz_dens <= 0 && info->flags&IIFsource) {
		horiz_unit_mm = 1.f;		// cheat to get density
		if (!strncmp(make, "OLYMPUS", 7)) {
			if (!strncmp(model, "C20", 3))
				horiz_dens = horiz_size / 6.14f;
			else if (!strncmp(model, "C30", 3))
				horiz_dens = horiz_size / 7.02f;
			else if (!strncmp(model, "C40", 3))
				horiz_dens = horiz_size / 7.02f;
			else if (!strncmp(model, "E-10", 4))
				horiz_dens = horiz_size / 8.75f;
			else if (!strncmp(model, "E-20", 4))
				horiz_dens = horiz_size / 8.75f;
		}
	}
	if (info->flags&IIFflength && (horiz_size > 0) & (horiz_dens > 0)) {
		info->hvangle = 360./M_PI * atan(.5*horiz_unit_mm/horiz_dens*
						horiz_size/info->flength);
		info->flags |= IIFhvangle;
	} else if (focal_35mm > 0) {
		info->hvangle = 360./M_PI * atan(17.5/(double)focal_35mm);
		info->flags |= IIFhvangle;
		if (!(info->flags&IIFflength) & (horiz_size > 0) &
						(horiz_dens > 0)) {
			info->flength = focal_35mm / 35.f *
					horiz_unit_mm/horiz_dens*horiz_size;
			info->flags |= IIFflength;
		}
	}
						// compute calibration factor
	if ((info->flags & (IIFexptime|IIFaperture|IIFasa|IIFstonits)) ==
					(IIFexptime|IIFaperture|IIFasa)) {
		static const float	Kfactor = 87.f;	// empirical
		info->stonits = Kfactor * info->aperture*info->aperture /
				(info->exptime * info->asa);
		info->flags |= IIFstonits;
	}
						// get relevant GPS tags
	if (GPSOffset == 0 || !ti->SeekIFD(GPSOffset))
		return false;
	int	latSign = 0, longSign = 0, altSign = 0;
	float	utc_hms[3] = {-1.f};
	float	dms_buf[3];
	char	rbuf[16];
	info->gmtdate[0] = '\0';
	for (tag = ti->FirstTag(); tag >= 0; tag = ti->NextTag())
		switch (tag) {
		case 0x01:			// Latitude reference
			if (ti->GetData(rbuf, TIFFascii_string, sizeof(rbuf)))
				latSign = (rbuf[0] == 'S') ? -1 : 1;
			break;
		case 0x02:			// Latitude
			if (!ti->GetData(dms_buf, TIFFunsigned_rational, 3))
				break;
			info->latlong[0] = dms_buf[0] + dms_buf[1]*(1.f/60.f) +
						dms_buf[2]*(1.f/3600.f);
			if (!latSign) latSign = 1;
			break;
		case 0x03:			// Longitude reference
			if (ti->GetData(rbuf, TIFFascii_string, sizeof(rbuf)))
				longSign = (rbuf[0] == 'W') ? -1 : 1;
			break;
		case 0x04:			// Longitude
			if (!ti->GetData(dms_buf, TIFFunsigned_rational, 3))
				break;
			info->latlong[1] = dms_buf[0] + dms_buf[1]*(1.f/60.f) +
						dms_buf[2]*(1.f/3600.f);
			if (!longSign) longSign = 1;
			break;
		case 0x05:			// Altitude reference
			if (ti->GetData(rbuf, TIFFunsigned_byte, 1))
				altSign = (rbuf[0] == '1') ? -1 : 1;
			break;
		case 0x06:			// Altitude
			info->altitude = ti->GetReal();
			if (!altSign) altSign = 1;
			break;
		case 0x07:			// UTC hours, minutes, seconds
			if (!ti->GetData(utc_hms, TIFFunsigned_rational, 3))
				utc_hms[0] = -1.f;
			break;
		case 0x1d:			// GMT year, month, day
			ti->GetData(info->gmtdate, TIFFascii_string,
					sizeof(info->gmtdate));
			break;
		case 0x18:			// Sight bearing
			info->bearing = ti->GetReal();
			info->flags |= IIFbearing;
			break;
		}
	if ((latSign != 0) & (longSign != 0)) {
		info->latlong[0] *= float(latSign);
		info->latlong[1] *= float(longSign);
		info->flags |= IIFlatlong;
	}
	if (altSign != 0) {
		info->altitude *= float(altSign);
		info->flags |= IIFaltitude;
	}
	if (info->gmtdate[0] && utc_hms[0] >= 0) {
		sprintf(info->gmtdate+10, " %02d:%02d:%02d",
				int(utc_hms[0]+.5f),
				int(utc_hms[1]+.5f),
				int(utc_hms[2]+.5f));
		if (dateOK(info->gmtdate))
			info->flags |= IIFgmtdate;
	}
	return true;
#undef make
#undef model
#undef rev
#undef artist
}
