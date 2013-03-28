
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
#include <float.h>
#include <limits.h>
#include <mpfr.h>
#include <mpf2mpfr.h>



#define BERTINI_REAL_VERSION_STRING "1.0"



#ifndef _PROGRAM_STARTUP_H
#define _PROGRAM_STARTUP_H

#include "polysolve.h"
#include "fileops.h"
#include "data_type.h"


int startup(int argC, char *args[], char **inputName, char **startName);
void get_tracker_config(tracker_config_t *T,int MPType);

#endif


