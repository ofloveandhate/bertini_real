#ifndef _SAMPLER_H
#define _SAMPLER_H


/** \file sampler.hpp */



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









/**
 \defgroup samplermethods Sampler Methods
	
 Methods common to sampling Bertini_real decompositions.
*/



/**
 \brief light startup processing common to all dimensions, for the sampler.
 
 parse the input file, get the tracker_config_t, set up the PPD, adjust a few sampler settings.
 
 \ingroup samplermethods
 
 \param D the base-class-type decomposition to refine.
 \param sampler_options The current state of sampler.
 \param solve_options the current state of the solver configuration.
 */
void common_sampler_startup(const decomposition & D,
							sampler_configuration & sampler_options,
							solver_configuration & solve_options);




/**
 \brief get the MPType, name of the directory to sample, and the dimension of the decomposition.
 
 \todo replace this function which reads decomposition metadata from a file in the directory.
 
 \return the MPType
 \param Dir_Name the name read in by this function
 \param MPType apparently this function returns this in two ways.
 \param dimension The retrieved dimension.
 */
int get_dir_mptype_dimen(boost::filesystem::path & Dir_Name, int & MPType, int & dimension);




















/**
 \brief Set the linear and point in a witness set to the input L and pts.
 
 \todo remove this function, or make a method of the witness_set class.
 
 \param W the witness set to modify.
 \param L the linear to set
 \param pts the input point
 \param num_vars The number of vars to change the W to.
 
 */
void set_witness_set_mp(witness_set *W, vec_mp L,vec_mp pts,int num_vars);







/**
 \brief given two points and a projection, estimate a projection value for the point halfway between.
 
 \param result the computed estimated projection value 
 \param left Input one.
 \param right Input two.
 \param pi the linear projection to use.
 */
void estimate_new_projection_value(comp_mp result, vec_mp left, vec_mp right, vec_mp pi);


/**
 \brief construct the triangulation of two ribs of a face
 
 \todo add a diagram illustrating this method
 
 \param rib1 input vector of indices into V
 \param rib2 the other input vector of indices into V
 \param V the vertex set into which the ribs index.  it holds the points and their projection values.
 \param current_samples The constructed triangles from this method.
 */
void triangulate_two_ribs(const std::vector< int > & rib1, const std::vector< int > & rib2,
						  const vertex_set & V,
						  std::vector< triangle> & current_samples);


/**
 \brief construct the triangulation of two ribs of a face, when the ribs have an unequal number of points on them
 
 \todo add a diagram illustrating this method
 
 \param rib1 input vector of indices into V
 \param rib2 the other input vector of indices into V
 \param V the vertex set into which the ribs index.  it holds the points and their projection values.
 \param current_samples The constructed triangles from this method.
 */
void triangulate_two_ribs_by_binning(const std::vector< int > & rib1, const std::vector< int > & rib2,
									 const vertex_set & V,
									 std::vector< triangle> & current_samples);





#endif




