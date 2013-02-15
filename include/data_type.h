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


#include "polysolve.h"

#ifndef _DATA_TYPE_H
#define _DATA_TYPE_H

/*** low-level data types. ***/




typedef struct
{
  int num_pts;
  point_d *pts;
  int dim;
  int deg;
} witness_point_set_d;
typedef struct
{
  char *func;  //symbolic representation of function (straight from input file).
} function_d;

typedef struct
{
  function_d *funcs; //probably not used for now.
  //prog_t slp; //SLP -- we'll use this primarily, at least at first.
  char *fName; //system can be given in a file -- this is its name.
}system_d;

typedef struct
{
  system_d sys;
  vec_d L;
  witness_point_set_d W;
}witness_set_d;  //For a single irreducible component!  Need num. irred. decomp. type later.


// CURVE CELL DECOMP DATA TYPES
typedef struct
{
  point_d pt;
  comp_d projVal; //Value of projection pi applied to pt; used for easy comparison of projected values later....
  int type;  //See enum above.
  int num_left;  //this and next line track how often this vertex appears in
  int num_right; //edges;  good for sanity check AND determining V0.
}vertex_d;

typedef struct
{
  int left;  //index from V1
  int right; //index from V1
  int midpt; //index from midPts
  witness_set_d W; //contains functions; IS THIS OVERKILL????  Could just be a system_d....
  vec_d pi;  //projection
}edge_d;

typedef struct
{
  vertex_d *V0;  //Isolated real points.
  vertex_d *V1;  //Critical points AND new non-critical endpoints of edges.
  vertex_d *midPts;  //Midpoints of edges.
  edge_d *E;
  int      num_V0;
  int      num_V1;
  int      num_midPts;
  int      num_E;
}curveDecomp_d;

//The following lets us use words instead of numbers to indicate vertex types.
enum {CRITICAL=0, NEW=1, MIDPOINT=2};




#endif

