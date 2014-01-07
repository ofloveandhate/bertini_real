

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


#include "missing_bertini_headers.hpp"


#include "fileops.hpp"
#include "data_type.hpp"
#include "programConfiguration.hpp"
#include "checkSelfConjugate.hpp"

#include "solver_lintolin.hpp"
#include "solver_multilintolin.hpp"
#include "surface.hpp"
#include "curve.hpp"

#ifndef _SAMPLER_H
#define _SAMPLER_H





class curve_sample_data : public curve_decomposition
{
public:
	
	int num_variables;
	int num_edges;
	
	int *num_samples_each_edge;
	
	int **sample_indices;
};

void clear_sample(curve_sample_data *S, int MPType);





int  curve_sampler_startup(boost::filesystem::path directoryName,
						   boost::filesystem::path &inputName,
						   boost::filesystem::path &witnessSetName,
						   boost::filesystem::path &RandMatName,
						   boost::filesystem::path &samplingNamenew,
						   curve_decomposition &C,
						   vertex_set &V);





int get_dir_mptype(boost::filesystem::path & Dir_Name, int * MPType);





int  set_initial_sample_data(curve_sample_data *S, curve_decomposition C, vertex_set V,
							 int num_vars);

void  output_sampling_data(curve_sample_data S,  vertex_set V,
						   boost::filesystem::path samplingName,int num_vars, int MPType);
void  set_witness_set_d(witness_set *W, vec_d L,vec_d pts,int num_vars);
void set_witness_set_mp(witness_set *W, vec_mp L,vec_mp pts,int num_vars);

void generate_new_sampling_pts(curve_sample_data *S_new,
							   mat_mp randomizer_matrix,
							   curve_sample_data S_old,
							   curve_decomposition C,
							   vertex_set &V,
							   witness_set & W,
							   int  MPType,
							   sampler_configuration *sampler_options,
							   solver_configuration & solve_options);






void estimate_new_projection_value(comp_mp result, vec_mp left, vec_mp right, vec_mp pi);
#endif


void set_initial_refinement_flags(int *num_refinements, int **refine_flags, int **current_indices,
								  curve_sample_data S, vertex_set &V,
								  int current_edge, sampler_configuration *sampler_options);


