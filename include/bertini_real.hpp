

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
#include <sstream>




#ifndef _BERTINI_REAL_H
#define _BERTINI_REAL_H

extern "C" {
#include "cascade.h"
}
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
#include "surface.hpp"
#include "curve.hpp"
#include "parallelism.hpp"

#endif



