

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
	
	std::vector<int> num_samples_each_edge;
	
	std::vector< std::vector<int >> sample_indices;
	
	~curve_sample_data()
	{
		
		clear();
	}
	
	
private:
	void clear()
	{
		
	}
};


//void clear_sample(curve_sample_data *S, int MPType);




void common_sampler_startup(const decomposition & D);





int get_dir_mptype_dimen(boost::filesystem::path & Dir_Name, int & MPType, int & dimension);




















void  output_sampling_data(curve_sample_data S,  vertex_set V,
						   boost::filesystem::path samplingName,int num_vars, int MPType);
void  set_witness_set_d(witness_set *W, vec_d L,vec_d pts,int num_vars);
void set_witness_set_mp(witness_set *W, vec_mp L,vec_mp pts,int num_vars);

void curve_generate_new_sampling_pts_adaptive(curve_sample_data *S_new,
							   mat_mp randomizer_matrix,
							   curve_sample_data & S_old,
							   curve_decomposition & C,
							   vertex_set &V,
							   witness_set & W,
							   int  MPType,
							   sampler_configuration & sampler_options,
							   solver_configuration & solve_options);






void estimate_new_projection_value(comp_mp result, vec_mp left, vec_mp right, vec_mp pi);



int set_initial_sample_data(curve_sample_data & S, const curve_decomposition & C, vertex_set V,
							 int num_vars);

void set_initial_refinement_flags(int & num_refinements, std::vector<bool> & refine_flags, std::vector<int> & current_indices,
                                  const curve_sample_data & S, vertex_set &V,
                                  int current_edge, sampler_configuration & sampler_options);




#endif




