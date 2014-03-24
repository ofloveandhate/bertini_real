
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


#include <set>


#ifndef CURVE_H
#define CURVE_H



#include "bertini_headers.hpp"

#include "checkSelfConjugate.hpp"

#include "fileops.hpp"

#include "data_type.hpp"
#include "isosingular.hpp"
#include "programConfiguration.hpp"
#include "output.hpp"


#include "solver_lintolin.hpp"
#include "nullspace_left.hpp"
#include "solver_sphere.hpp"

/**
 a curve decomposition.
 
 includes methods to add vertices, look up vertices, etc
 */
class curve_decomposition : public decomposition
{
	
	std::vector<int> num_samples_each_edge;
	
	std::vector< std::vector<int >> sample_indices;
	
	
	friend class surface_decomposition;
	
public:
	
	std::vector<edge> edges;
	int      num_edges;
	
	
	
	
	
	
	std::set< int > all_edge_indices()
    {
        std::set< int > ind;
        
        for (int ii=0; ii<num_edges; ii++) {
            ind.insert(edges[ii].left);
            ind.insert(edges[ii].midpt);
            ind.insert(edges[ii].right);
        }
        
        return ind;
    }
	
	
	
	
	
	int add_edge(edge new_edge)
	{
		num_edges++;
		edges.push_back(new_edge);
		return num_edges-1; // -1 to correct for the fencepost problem
	}
	
	
	int setup(boost::filesystem::path containing_folder);
	
	int setup_edges(boost::filesystem::path INfile);
	
	void print_edges(boost::filesystem::path outputfile);
	
	void print(boost::filesystem::path base);
	
	
	
	
	
	
	
	
	int nondegenerate_edge_w_midpt(int ind)
	{
		
		for (int ii=0; ii<num_edges; ii++){
			if (this->edges[ii].midpt == ind){
				
				if (edges[ii].is_degenerate()) {
					continue;
				}
				
				return ii;
			}
		}
		
		return -10;
	}
	
	
	int nondegenerate_edge_w_left(int ind)
	{
		
		for (int ii=0; ii<num_edges; ii++){
			if (this->edges[ii].left == ind){
				
				if (edges[ii].is_degenerate()) {
					continue;
				}
				
				return ii;
			}
		}
		
		return -11;
	}
	int nondegenerate_edge_w_right(int ind)
	{
		
		for (int ii=0; ii<num_edges; ii++){
			if (this->edges[ii].right == ind){
				
				if (edges[ii].is_degenerate()) {
					continue;
				}
				
				return ii;
			}
		}
		
		return -12;
	}
	
	
	
	
	int edge_w_midpt(int ind)
	{
		
		for (int ii=0; ii<num_edges; ii++){
			if (this->edges[ii].midpt == ind){
				return ii;
			}
		}
		
		return -10;
	}
	
	
	
	int edge_w_left(int ind)
	{
		
		for (int ii=0; ii<num_edges; ii++){
			if (this->edges[ii].left == ind){
				return ii;
			}
		}
		
		return -11;
	}
	int edge_w_right(int ind)
	{
		
		for (int ii=0; ii<num_edges; ii++){
			if (this->edges[ii].right == ind){
				return ii;
			}
		}
		
		return -12;
	}
	
	
	int edge_w_removed(int ind)
	{
		
		for (int ii=0; ii<num_edges; ii++){
			
			for (std::vector<int>::iterator iter=edges[ii].removed_points.begin();
				 iter!=edges[ii].removed_points.end(); ++iter) {
				if (*iter == ind){
					return ii;
				}
			}
			
			
		}
		
		return -13;
	}
	
	
	std::vector<int> get_merge_candidate(const vertex_set & V);
	
	void merge(witness_set & W_midpt, vertex_set & V,
			   vec_mp * pi_in,
			   solver_configuration & solve_options);
	
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
	void 	computeCurveNotSelfConj(const witness_set		& W,
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
	void computeCurveSelfConj(const witness_set & W,
							  vec_mp *pi,
							  vertex_set &V,
							  BR_configuration & options,
							  solver_configuration & solve_options);
	
	
	
	int interslice(const witness_set & W_curve,
				   const witness_set & W_crit_real,
				   vec_mp *pi,
				   BR_configuration & program_options,
				   solver_configuration & solve_options,
				   vertex_set & V);
	
	
	
	
	int compute_critical_points(const witness_set & W_curve,
								vec_mp *pi,
								BR_configuration & program_options,
								solver_configuration & solve_options,
								witness_set & W_crit_real);
	
	
	
	
	
	
	int get_additional_critpts(witness_set *W_crit_real,
							   const witness_set & W,
							   BR_configuration & program_options,
							   solver_configuration & solve_options);
	
	
	void send(int target, parallelism_config & mpi_config);
	
	void receive(int source, parallelism_config & mpi_config);
	
	
	
	
	
	
	
	int adaptive_set_initial_sample_data();
	void adaptive_sampler(vertex_set &V,
						  sampler_configuration & sampler_options,
						  solver_configuration & solve_options);
	
	void adaptive_set_initial_refinement_flags(int & num_refinements, std::vector<bool> & refine_flags, std::vector<int> & current_indices,
											   vertex_set &V,
											   int current_edge, sampler_configuration & sampler_options);
	
	int fixed_set_initial_sample_data(int target_num_samples);
	
	void fixed_sampler(vertex_set &V,
					   sampler_configuration & sampler_options,
					   solver_configuration & solve_options,
					   int target_num_samples);
	
	void  output_sampling_data(boost::filesystem::path samplingName);
	
	
	
	
	
	
	
	
	void reset()
	{
		
		
		decomposition::reset();
		
		this->clear();
		this->init();
		
		
	}
	
	curve_decomposition() : decomposition()
	{
		this->init();
	}
	
	curve_decomposition & operator=(const curve_decomposition& other){
		this->init();
		this->copy(other);
		return *this;
	}
	
	curve_decomposition(const curve_decomposition & other){
		this->init();
		this->copy(other);
	}
	
	~curve_decomposition()
	{
		this->clear();
	}
	
	
	
	
protected:
	
	void clear()
	{
		edges.clear();
		num_edges = 0;
	}
	
	void init(){
		
		num_edges = 0;
		dimension = 1;
	}
	
	
	void copy(const curve_decomposition & other)
	{
		decomposition::copy(other);
		this->edges = other.edges;
		this->num_edges = other.num_edges;
	}
	
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
								 const witness_set & W);












int curve_get_additional_critpts(witness_set *W_crit_real,
								 const witness_set & W,
								 BR_configuration & program_options,
								 solver_configuration & solve_options);



/**
 read the file "deg.out" and takes the sum of the numbers appearing there. used for determining the number of lintolin solves to perform to get the critical points WRT the projection and coordinate axes (bounding box).
 */
int get_sum_degrees(char filename[], int num_funcs);











int verify_projection_ok(const witness_set & W,
						 vec_mp * projection,
						 solver_configuration & solve_options);

int verify_projection_ok(const witness_set & W,
						 system_randomizer * randomizer,
						 vec_mp * projection,
						 solver_configuration & solve_options);





#endif
