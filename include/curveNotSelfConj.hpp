



#ifndef _CURVE_NOTSELFCONJ_H
#define _CURVE_NOTSELFCONJ_H

#include <gmp.h>
extern "C" {
#include "polysolve.h"
}



#include "programConfiguration.hpp"
#include "witnessSet.hpp"
#include "data_type.hpp"
extern "C" {
#include "partitionParse.h"
}

#include "fileops.hpp"
#include "missing_bertini_headers.hpp"


/**
 the main function for computing the self-intersection of a non-self-conjugate dimension 1 component.
 
 \param W		the witness set 
 \param pi	the projection for this decomposition
 \param C		one of the structures to hold the generated data.  indexes into V
 \param V		the vertex set being passed around.  C indexes into here.
 \param num_vars		the total number of variables in the problem.
 \param input_file		the name of the input file to use.
 \param program_options		main structure holding configuration
 \param	solve_options			structure holding options to pass to a solver.
 
 */
void 	computeCurveNotSelfConj(witness_set		W,
															vec_mp				pi,
															curveDecomp_d *C,
															vertex_set		*V,
															int						num_vars,
															char					*input_file,
															program_configuration *program_options,
															solver_configuration * solve_options);



void  get_random_mat_d(mat_d, int,int);

/** 
 writes the input file for the diagonal homotopy used in curve case to find the intersection points.
 
 \param outputFile the name of the file to create
 \param funcInputx	the name of the first input file
 \param funcInputy	the name of the second input file
 \param configInput	the name of the config file to use
 \param L			the linears for the problem
 \param num_vars		the number of variables in the problem, including the homogeneous ones.
 */
void 	diag_homotopy_input_file(char  *outputFile,
															 char  *funcInputx,
															 char  *funcInputy,
															 char  *configInput,
															 vec_d L,
															 int   num_vars);




/**
 writes the start file for the diagonal homotopy used in curve case to find the intersection points.
 
 \param startFile the name of the start file to write
 \param W the input witness set
 */
void 	diag_homotopy_start_file(char					*startFile,
															 witness_set  W);






#endif
