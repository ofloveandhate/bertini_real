
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





#ifndef CURVE_H
#define CURVE_H



#include "missing_bertini_headers.hpp"

#include "checkSelfConjugate.hpp"

#include "fileops.hpp"

#include "data_type.hpp"
#include "isosingular.hpp"
#include "programConfiguration.hpp"
#include "output.hpp"


#include "solver_lintolin.hpp"
#include "solver_linprodtodetjac.hpp"
#include "solver_detjactodetjac.hpp"

#include "nullspace_left.hpp"


/**
 a curve decomposition.
 
 includes methods to add vertices, look up vertices, etc
 */
class curve_decomposition : public decomposition
{
public:
	
	std::vector<edge> edges;
	int      num_edges;
	
	
	
	
	
	
	
	void add_edge(edge new_edge);
	
	int setup_edges(boost::filesystem::path INfile);
	
	void print_edges(boost::filesystem::path outputfile);
	
	void print(boost::filesystem::path base);
	
	curve_decomposition() : decomposition()
	{
		init();
	}
	
	curve_decomposition & operator=(const curve_decomposition& other){
		init();
		copy(other);
		return *this;
	}
	
	curve_decomposition(const curve_decomposition & other){
		init();
		copy(other);
	}
	
	void clear()
	{
		decomposition::clear();
		
		edges.clear();
		num_edges = 0;
	}
	
	void init(){
		decomposition::init();
		num_edges = 0;
		dimension = 1;
	}
	
	
	void copy(const curve_decomposition & other)
	{
		decomposition::copy(other);
		this->edges = other.edges;
		this->num_edges = other.num_edges;
	}
	
	
	
	
	bool is_degenerate()
	{
		if (num_edges==1) {
			if (edges[0].left == edges[0].right) // can it be otherwise?
				return true;
			else
				return false;
		}
		else
			return false;
	}
	
	int edge_w_midpt(int ind)
	{
		
		for (int ii=0; ii<num_edges; ii++){
			if (this->edges[ii].midpt == ind){
				return ii;
			}
		}

		return -1;
	}
	
	
	int edge_w_left(int ind)
	{
		
		for (int ii=0; ii<num_edges; ii++){
			if (this->edges[ii].left == ind){
				return ii;
			}
		}
		
		return -1;
	}
	int edge_w_right(int ind)
	{
		
		for (int ii=0; ii<num_edges; ii++){
			if (this->edges[ii].right == ind){
				return ii;
			}
		}
		
		return -1;
	}
	
	
	void main(vertex_set & V,
						witness_set & W,
						vec_mp *pi,
						BR_configuration & program_options,
						solver_configuration & solve_options);
	
	
	
	
	
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
	void 	computeCurveNotSelfConj(witness_set		& W,
																vec_mp				pi,
																vertex_set		&V,
																int						num_vars,
																BR_configuration & program_options,
																solver_configuration & solve_options);
	
	
	
	/**
	 //the main function for computing cell decom for a curve.  only for use on a self-conjugate component.
	 
	 \param inputFile the name of the input file.
	 \param W	witness_set containing linears, patches, points, etc.  much info and calculation performed from this little guy.
	 \param pi the set of projections to use.  in this case, there should be only 1.
	 \param C		curve decomposition structure into which to place computed data.
	 \param V		vertex set structure into which to place collected data.
	 \param num_vars		the total number of variables for the problem.
	 \param options program configuration.
	 \param solve_options solver configuration.
	 */
	void computeCurveSelfConj(witness_set & W,
														vec_mp *pi,
														vertex_set &V,
														int num_vars,
														BR_configuration & options,
														solver_configuration & solve_options);
	
	
	
	int interslice(witness_set & W_curve,
								 witness_set & W_crit_real,
								 mat_mp randomizer_matrix,
								 vec_mp *pi,
								 BR_configuration & program_options,
								 solver_configuration solve_options,
								 vertex_set & V);
}; // end curve_decomposition










//void  get_random_mat_d(mat_d, int,int);

/**
 writes the input file for the diagonal homotopy used in curve case to find the intersection points.
 
 \param outputFile the name of the file to create
 \param funcInputx	the name of the first input file
 \param funcInputy	the name of the second input file
 \param configInput	the name of the config file to use
 \param L			the linears for the problem
 \param num_vars		the number of variables in the problem, including the homogeneous ones.
 */
void 	diag_homotopy_input_file(boost::filesystem::path outputFile,
															 boost::filesystem::path funcInputx,
															 boost::filesystem::path funcInputy,
															 boost::filesystem::path configInput,
															 vec_mp L,
															 int   num_vars);




/**
 writes the start file for the diagonal homotopy used in curve case to find the intersection points.
 
 \param startFile the name of the start file to write
 \param W the input witness set
 */
void 	diag_homotopy_start_file(boost::filesystem::path startFile,
															 witness_set & W);









int curve_compute_critical_points(witness_set & W_curve,
																	mat_mp randomizer_matrix,
																	int *randomized_degrees,
																	vec_mp *pi,
																	BR_configuration & program_options,
																	solver_configuration solve_options,
																	witness_set & W_crit_real);
/**
 the linprodtodetjac method for getting the critical points
 */
int compute_crit_linprodtodetjac(witness_set *W_crit_real, // the returned value
																 witness_set & W,
																 mat_mp n_minusone_randomizer_matrix,
																 vec_mp pi,
																 int num_new_linears,
																 BR_configuration & program_options,
																 solver_configuration & solve_options);




int curve_get_additional_critpts(witness_set *W_crit_real,
																 witness_set & W,
																 mat_mp randomizer_matrix,
																 vec_mp pi,
																 int *randomized_degrees,
																 BR_configuration & program_options,
																 solver_configuration & solve_options);

//
///**
//
// */
//void check_patch_values(witness_set W);

///**
// checks to see what the determinant of the jacobian is at the points in W.
// */
//void check_detjac(witness_set W, prog_t SLP, tracker_config_t T, mat_d n_minusone_randomizer_matrix, vec_d projection);
//
///**
//// gets the jacobian (homogeneous) of the functions at the point current_values.  returns mat_d jacobian.  primarily for testing.
// */
//void get_jacobian(point_d current_values,
//									int MPType,
//									int num_var_gps,
//									prog_t SLP,
//									tracker_config_t T,
//									mat_d jacobian);

/**
 read the file "deg.out" and takes the sum of the numbers appearing there. used for determining the number of lintolin solves to perform to get the critical points WRT the projection and coordinate axes (bounding box).
 */
int get_sum_degrees(char filename[], int num_funcs);


void sort_for_membership(char * input_file,
												 witness_set *W_out,
												 witness_set & W_in,
												 char *stifle_text);









int verify_projection_ok(witness_set & W,
												 vec_mp projection,
												 solver_configuration & solve_options);

int verify_projection_ok(witness_set & W,
												 mat_mp randomizer_matrix,
												 vec_mp projection,
												 solver_configuration & solve_options);





#endif
