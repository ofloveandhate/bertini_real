

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


#include "bertini_headers.hpp"


#include "fileops.hpp"
#include "data_type.hpp"
#include "programConfiguration.hpp"
#include "checkSelfConjugate.hpp"

#include "solver_multilintolin.hpp"
#include "surface.hpp"
#include "curve.hpp"

#ifndef _SAMPLER_H
#define _SAMPLER_H












void common_sampler_startup(const decomposition & D,
							sampler_configuration & sampler_options,
							solver_configuration & solve_options);





int get_dir_mptype_dimen(boost::filesystem::path & Dir_Name, int & MPType, int & dimension);





















void  set_witness_set_d(witness_set *W, vec_d L,vec_d pts,int num_vars);
void set_witness_set_mp(witness_set *W, vec_mp L,vec_mp pts,int num_vars);








void estimate_new_projection_value(comp_mp result, vec_mp left, vec_mp right, vec_mp pi);

void triangulate_two_ribs(const std::vector< int > & rib1, const std::vector< int > & rib2,
						  const vertex_set & V,
						  std::vector< triangle> & current_samples);

void triangulate_two_ribs_by_binning(const std::vector< int > & rib1, const std::vector< int > & rib2,
									 const vertex_set & V,
									 std::vector< triangle> & current_samples);


//int set_initial_sample_data(curve_sample_data & S, const curve_decomposition & C, vertex_set V,
//							int num_vars);
//
//void set_initial_refinement_flags(int & num_refinements, std::vector<bool> & refine_flags, std::vector<int> & current_indices,
//                                  vertex_set &V,
//                                  int current_edge, sampler_configuration & sampler_options);




#endif




