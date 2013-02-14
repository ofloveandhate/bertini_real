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

#ifndef _DATA_TYPE_H
#define _DATA_TYPE_H

/*** low-level data types. ***/
typedef struct
{
  double r, i;
} _comp_d;

typedef struct
{
  mpf_t r, i;
} _comp_mp;

typedef struct
{
  _comp_d *coord;
  int alloc_size; // allocated size
  int size;       // size of the point
} _point_d;

typedef struct
{
  _comp_mp *coord;
  int alloc_size; // allocated size
  int curr_prec;  // current precision
  int size;       // size of the point
} _point_mp;

typedef struct
{
  _comp_d **entry;
  int alloc_rows; // allocated number of rows
  int alloc_cols; // allocated number of cols
  int rows;       // number of rows for the matrix
  int cols;       // number of cols for the matrix
} _mat_d;

typedef struct
{
  _comp_mp **entry;
  int alloc_rows; // allocated number of rows
  int alloc_cols; // allocated number of cols
  int curr_prec;  // current precision
  int rows;       // number of rows for the matrix
  int cols;       // number of cols for the matrix
} _mat_mp;

typedef _comp_d   comp_d[1];   /* complex number */
typedef _point_d  point_d[1];  /* complex point */
typedef  point_d  vec_d;       /* complex vector (same structure as a point - useful to have two names) */
typedef _mat_d    mat_d[1];    /* complex matrix */

typedef _comp_mp  comp_mp[1];
typedef _point_mp point_mp[1];
typedef  point_mp vec_mp;
typedef _mat_mp   mat_mp[1];

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

//init
#define init_point_d(_a, _s) { (_a)->coord = (_comp_d *)bmalloc((_s) * sizeof(_comp_d)); (_a)->alloc_size = _s; (_a)->size = 0; }
#define init_vec_d   init_point_d
#define init_mat_d(_a, _r, _c) { int _i; (_a)->entry = (_comp_d **)bmalloc((_r) * sizeof(_comp_d *)); \
for (_i = 0; _i < _r; _i++) (_a)->entry[_i] = (_comp_d *)bmalloc((_c) * sizeof(_comp_d)); (_a)->alloc_rows = _r; (_a)->alloc_cols = _c; (_a)->rows = (_a)->cols = 0; }
// set real & imag values
#define set_double_d(_r, _x, _y) { (_r)->r = _x; (_r)->i = _y; }

// set zero
#define set_zero_d(_r) { (_r)->r = (_r)->i = 0; }

#define d_abs_d(_a)  (sqrt((_a)->r*(_a)->r + (_a)->i*(_a)->i))
#define add_d(_r,_a,_b)  { (_r)->r = (_a)->r + (_b)->r; (_r)->i = (_a)->i + (_b)->i; }
#define sub_d(_r,_a,_b)  { (_r)->r = (_a)->r - (_b)->r; (_r)->i = (_a)->i - (_b)->i; }
#define mul_d(_r,_a,_b) { double _s = (_a)->r * (_b)->r - (_a)->i * (_b)->i; \
(_r)->i = (_a)->r * (_b)->i + (_a)->i * (_b)->r; (_r)->r = _s;}
#define mul_d2(_r,_a,_b,_s) { _s = (_a)->r * (_b)->r - (_a)->i * (_b)->i; (_r)->i = (_a)->r * (_b)->i + (_a)->i * (_b)->r; (_r)->r = _s;}
#define vec_sub_d(_r, _u, _v) { int _i, _size = (_u)->size; increase_size_vec_d(_r, _size); \
(_r)->size = _size; for (_i = 0; _i < _size; _i++) sub_d(&(_r)->coord[_i], &(_u)->coord[_i], &(_v)->coord[_i]); }
// increase size
#define increase_size_point_d(_a, _new_size) { if ((_a)->alloc_size < _new_size) { (_a)->coord = (_comp_d *)brealloc((_a)->coord, (_new_size) * sizeof(_comp_d)); (_a)->alloc_size = _new_size; }}

// copy value
#define set_d(_r,_a)   { (_r)->r = (_a)->r; (_r)->i = (_a)->i; }
#define point_cp_d(_d, _s) { int _i, _size = (_s)->size; increase_size_point_d(_d, _size); \
(_d)->size = _size; for (_i = 0; _i < _size; _i++) set_d(&(_d)->coord[_i], &(_s)->coord[_i]); }

// clear
#define clear_point_d(_a)  { free((_a)->coord); (_a)->coord = NULL; (_a)->alloc_size = (_a)->size = 0; }
#define clear_vec_d   clear_point_d
#define clear_mat_d(_a)  { int _i; for (_i = (_a)->alloc_rows - 1; _i >= 0; _i--) { free((_a)->entry[_i]); } free((_a)->entry); (_a)->entry = NULL; \
(_a)->alloc_rows = (_a)->alloc_cols = (_a)->rows = (_a)->cols = 0; }






// divide
#define div_d(_r,_a,_b) { double _s = 1 / ((_b)->r*(_b)->r + (_b)->i*(_b)->i);\
double _rt = _s * ((_a)->r*(_b)->r + (_a)->i*(_b)->i);\
(_r)->i = _s * ((_a)->i*(_b)->r - (_a)->r*(_b)->i); (_r)->r = _rt; }
#define div_d2(_r,_a,_b, _s, _rt) { _s = 1 / ((_b)->r*(_b)->r + (_b)->i*(_b)->i); _rt = _s * ((_a)->r*(_b)->r + (_a)->i*(_b)->i);\
(_r)->i = _s * ((_a)->i*(_b)->r - (_a)->r*(_b)->i); (_r)->r = _rt; }

#define div_mp(_r, _a, _b) { int _oid = thread_num(); \
mpf_mul(_tempMPF1[_oid], (_b)->r, (_b)->r); mpf_mul(_tempMPF2[_oid], (_b)->i, (_b)->i); mpf_add(_tempMPF1[_oid], _tempMPF1[_oid], _tempMPF2[_oid]); \
mpf_ui_div(_tempMPF1[_oid], 1, _tempMPF1[_oid]); \
mpf_mul(_tempMPF2[_oid], (_a)->r, (_b)->r); mpf_mul(_tempMPF3[_oid], (_a)->i, (_b)->i); \
mpf_add(_tempMPF2[_oid], _tempMPF2[_oid], _tempMPF3[_oid]); mpf_mul(_tempMPF2[_oid], _tempMPF2[_oid], _tempMPF1[_oid]); \
mpf_mul(_tempMPF3[_oid], (_a)->i, (_b)->r); mpf_mul((_r)->i, (_a)->r, (_b)->i); mpf_sub((_r)->i, _tempMPF3[_oid], (_r)->i); mpf_mul((_r)->i, (_r)->i, _tempMPF1[_oid]); \
mpf_set((_r)->r, _tempMPF2[_oid]); }
#define div_omp_mp(_r, _a, _b, _oid) { \
mpf_mul(_tempMPF1[_oid], (_b)->r, (_b)->r); mpf_mul(_tempMPF2[_oid], (_b)->i, (_b)->i); mpf_add(_tempMPF1[_oid], _tempMPF1[_oid], _tempMPF2[_oid]); \
mpf_ui_div(_tempMPF1[_oid], 1, _tempMPF1[_oid]); \
mpf_mul(_tempMPF2[_oid], (_a)->r, (_b)->r); mpf_mul(_tempMPF3[_oid], (_a)->i, (_b)->i); \
mpf_add(_tempMPF2[_oid], _tempMPF2[_oid], _tempMPF3[_oid]); mpf_mul(_tempMPF2[_oid], _tempMPF2[_oid], _tempMPF1[_oid]); \
mpf_mul(_tempMPF3[_oid], (_a)->i, (_b)->r); mpf_mul((_r)->i, (_a)->r, (_b)->i); mpf_sub((_r)->i, _tempMPF3[_oid], (_r)->i); mpf_mul((_r)->i, (_r)->i, _tempMPF1[_oid]); \
mpf_set((_r)->r, _tempMPF2[_oid]); }
#define div_mp2(_r, _a, _b) { int _oid = thread_num(), _prec = mpf_get_prec((_a)->r), _prec2 = mpf_get_prec((_b)->r), _curr_prec[3]; \
if (_prec > _prec2) _prec = _prec2;\
_curr_prec[0] = mpf_get_prec(_tempMPF1[_oid]); _curr_prec[1] = mpf_get_prec(_tempMPF2[_oid]); _curr_prec[2] = mpf_get_prec(_tempMPF3[_oid]);\
if (_prec > _curr_prec[0]) mpf_set_prec(_tempMPF1[_oid], _prec);\
if (_prec > _curr_prec[1]) mpf_set_prec(_tempMPF2[_oid], _prec);\
if (_prec > _curr_prec[2]) mpf_set_prec(_tempMPF3[_oid], _prec);\
mpf_mul(_tempMPF1[_oid], (_b)->r, (_b)->r); mpf_mul(_tempMPF2[_oid], (_b)->i, (_b)->i); mpf_add(_tempMPF1[_oid], _tempMPF1[_oid], _tempMPF2[_oid]); \
mpf_ui_div(_tempMPF1[_oid], 1, _tempMPF1[_oid]); \
mpf_mul(_tempMPF2[_oid], (_a)->r, (_b)->r); mpf_mul(_tempMPF3[_oid], (_a)->i, (_b)->i); \
mpf_add(_tempMPF2[_oid], _tempMPF2[_oid], _tempMPF3[_oid]); mpf_mul(_tempMPF2[_oid], _tempMPF2[_oid], _tempMPF1[_oid]); \
mpf_mul(_tempMPF3[_oid], (_a)->i, (_b)->r); mpf_mul((_r)->i, (_a)->r, (_b)->i); mpf_sub((_r)->i, _tempMPF3[_oid], (_r)->i); mpf_mul((_r)->i, (_r)->i, _tempMPF1[_oid]); \
mpf_set((_r)->r, _tempMPF2[_oid]); \
if (_prec > _curr_prec[0]) mpf_set_prec(_tempMPF1[_oid], _curr_prec[0]);\
if (_prec > _curr_prec[1]) mpf_set_prec(_tempMPF1[_oid], _curr_prec[1]);\
if (_prec > _curr_prec[2]) mpf_set_prec(_tempMPF1[_oid], _curr_prec[2]); }

// reciprocate
#define recip_d(_r, _a) { double _s = 1 / ((_a)->r*(_a)->r + (_a)->i*(_a)->i); (_r)->r = (_a)->r * _s; (_r)->i = - (_a)->i * _s; }
#define recip_d2(_r, _a, _s) { _s = 1 / ((_a)->r*(_a)->r + (_a)->i*(_a)->i); (_r)->r = (_a)->r * _s; (_r)->i = - (_a)->i * _s; }
#define recip_mp(_r, _a) { int _oid = thread_num(); \
mpf_mul(_tempMPF1[_oid], (_a)->r, (_a)->r); mpf_mul(_tempMPF2[_oid], (_a)->i, (_a)->i); \
mpf_add(_tempMPF1[_oid], _tempMPF1[_oid], _tempMPF2[_oid]); mpf_ui_div(_tempMPF1[_oid], 1, _tempMPF1[_oid]); \
mpf_mul((_r)->r, (_a)->r, _tempMPF1[_oid]); mpf_neg(_tempMPF1[_oid], _tempMPF1[_oid]); mpf_mul((_r)->i, (_a)->i, _tempMPF1[_oid]); }
#define recip_rat(_r, _a) { mpq_t _t1, _t2; mpq_init(_t1); mpq_init(_t2); \
mpq_mul(_t1, (_a)[0], (_a)[0]); mpq_mul(_t1, (_a)[1], (_a)[1]); mpq_add(_t1, _t1, _t1); mpq_inv(_t1, _t1); \
mpq_mul((_r)[0], (_r)[0], _t1); mpq_neg(_t1, _t1); mpq_mul((_r)[1], (_a)[1], _t1); mpq_clear(_t1); mpq_clear(_t2); }


#endif

