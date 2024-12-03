/*
 *  panopen.cpp -- replacement for open/close/lseek using wxWidgets
 *
 *  Created by Greg Ward on 1/28/10.
 *  Copyright 2010 Anyhere Software. All rights reserved.
 *
 */

#ifdef _WIN32

#include <wx/utils.h>
#include "system.h"

extern "C" {			// permits calls from C modules

int
unix_gethostname(char* hostname, size_t name_len)
{
	wxString ourhost = wxGetFullHostName();
	if (ourhost.IsEmpty())
		return -1;
	strncpy(hostname, ourhost.c_str(), name_len);
	return 0;
}

int
unix_kill(long pid, int sig)
{
	return wxKill(pid, (wxSignal) sig);
}

}	// end of extern "C" function definitions


#endif
