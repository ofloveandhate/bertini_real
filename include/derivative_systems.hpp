#ifndef DERIVATIVE_SYSTEMS_H
#define DERIVATIVE_SYSTEMS_H


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

#include <vector>
#include <map>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <sys/types.h>






#include "bertini_headers.hpp"



#include "programConfiguration.hpp"
#include "data_type.hpp"

#include "fileops.hpp"




void addItems(int *numItems, char ***itemNames, int **itemLines, FILE *IN, int lineNumber);

void parse_names(int *numItems, char ***itemNames, int **itemLines, FILE *IN, char *name, int num_declarations);

std::string just_constants(boost::filesystem::path filename,
																int numConstants, char **consts, int *lineConstants);


void write_matrix_as_constants(mat_mp M, std::string prefix, FILE *OUT);

void write_vector_as_constants(vec_mp V, std::string prefix, FILE *OUT);


#endif
