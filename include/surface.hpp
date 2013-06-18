
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

extern "C" {
#include "polysolve.h"
}

#include "fileops.hpp"

#include "checkSelfConjugate.hpp"
#include "curveNotSelfConj.hpp"
#include "curveSelfConj.hpp"
#include "data_type.hpp"
#include "isosingular.hpp"
#include "programConfiguration.hpp"
#include "witnessSet.hpp"
#include "output.hpp"
#include "missing_bertini_headers.hpp"

void surface_main(witness_set W,
									int self_conjugate,
									program_configuration *program_options,
									solver_configuration *solve_options);



#endif


