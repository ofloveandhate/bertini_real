
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





#ifndef SURFACE_H
#define SURFACE_H


#include "fileops.h"
#include "polysolve.h"
#include "checkSelfConjugate.h"
#include "curveNotSelfConj.h"
#include "curveSelfConj.h"
#include "data_type.h"
#include "isosingular.h"
#include "programConfiguration.h"
#include "witnessSet.h"
#include "output.h"
#include "missing_bertini_headers.h"

void surface_main(witness_set W,
									int self_conjugate,
									program_configuration *program_options,
									solver_configuration *solve_options);



#endif


