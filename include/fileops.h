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

#include <sys/stat.h>         /* declare the 'stat' structure */
#include <sys/types.h>


#include <dirent.h>




#ifndef _FILEOPS_H
#define _FILEOPS_H


#include "polysolve.h"


void purge_previous_directory(char *directoryName);

FILE *safe_fopen_read(char * filename);

FILE *safe_fopen_write(char * filename);

FILE *safe_fopen_append(char * filename);


#endif

