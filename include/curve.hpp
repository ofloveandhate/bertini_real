#ifndef CURVE_H
#define CURVE_H

/** \file curve.hpp */

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





#include "bertini_headers.hpp"

#include "checkSelfConjugate.hpp"

#include "fileops.hpp"

#include "data_type.hpp"
#include "isosingular.hpp"
#include "programConfiguration.hpp"
#include "output.hpp"


#include "nullspace_left.hpp"
#include "solver_sphere.hpp"

/**
 \brief A bertini_real cell curve decomposition.
 
 includes methods to add vertices, look up vertices, etc
 */
class curve_decomposition : public decomposition
{
	
	std::vector<int> num_samples_each_edge;
	
	std::vector< std::vector<int >> sample_indices;
	
	
	friend class surface_decomposition;
	
public:
	
	std::vector<edge> edges; ///< The edges (1-cells) computed by Bertini_real
	int      num_edges;  ///< How many edges this decomposition currently has.  This could also be inferred from edges.size()
	
	
	
	
	
	/**
	 \brief Get the set of all vertex indices to which this decomposition refers.
	 
	 This set is computed by reading all the edges.
	 
	 \return a std::set of ALL the vertex indices occuring in this decomposition.
	 */
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
	
	
	
	
	/**
	 \brief Commit an edge to the curve.
	 
	 Increments the edge counter, and pushes the edge to the pack of the vector of edges.
	 
	 \return the index of the edge just added.
	 \param new_edge the edge to add.
	 */
	int add_edge(edge new_edge)
	{
		num_edges++;
		edges.push_back(new_edge);
		return num_edges-1; // -1 to correct for the fencepost problem
	}
	
	/**
	 \brief Read a decomposition from text file.
	 
	 This function takes in the name of a folder containing a decomposition to set up, and parses the two folders "containing_folder/decomp" and "containing_folder/E.edge".  Note that the vertex_set is set up separately from decompositions.
	 
	 \return the number 1.  Seems stupid.
	 
	 \param containing_folder The name of the folder containing the decomposition.
	 */
	int setup(boost::filesystem::path containing_folder);
	
	
	
	/**
	 \brief open a folder as an edges file, and parse it into the vector of edges in this decomposition.
	 
	 \todo Add code ensuring the file is valid.  Probably would be best in xml.
	 
	 \return the number of edges added
	 \param INfile the name of the file to open and read as an edge file.
	 */
	int setup_edges(boost::filesystem::path INfile);
	
	
	/**
	 
	 \brief Write all the edges in this decomposition to a file.
	 
	 the format is 
	 
	 [
	 num_edges
	 
	 left mid right
	 ]
	 
	 \param outputfile The name of the file to write to.
	 */
	void print_edges(boost::filesystem::path outputfile);
	
	
	/**
	 \brief print the complete curve decomposition, including the base decomposition class.
	 
	 \param base The name of the base folder to print the decomposition to.
	 */
	void print(boost::filesystem::path base);
	
	
	
	
	
	
	
	/**
	 \brief Search this decomposition for an nondegenerate edge with input point index as midpoint.
	 
	 \param ind The index to search for.
	 \return The index of the found edge, or -10 if no edge had the point as midpoint.
	 */
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
	
	
	
	/**
	 \brief Search this decomposition for an nondegenerate edge with input point index as left point.
	 
	 \param ind The index to search for.
	 \return The index of the found edge, or -11 if no edge had the point as left point.
	 */
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
	
	
	/**
	 \brief Search this decomposition for an nondegenerate edge with input point index as right point.
	 
	 \param ind The index to search for.
	 \return The index of the found edge, or -12 if no edge had the point as right point.
	 */
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
	
	
	
	/**
	 \brief Search this decomposition for any edge (degenerate or not) with input point index as mid point.
	 
	 \param ind The index to search for.
	 \return The index of the found edge, or -10 if no edge had the point as mid point.
	 */
	int edge_w_midpt(int ind)
	{
		
		for (int ii=0; ii<num_edges; ii++){
			if (this->edges[ii].midpt == ind){
				return ii;
			}
		}
		
		return -10;
	}
	
	
	/**
	 \brief Search this decomposition for any edge (degenerate or not) with input point index as left point.
	 
	 \param ind The index to search for.
	 \return The index of the found edge, or -11 if no edge had the point as left point.
	 */
	int edge_w_left(int ind)
	{
		
		for (int ii=0; ii<num_edges; ii++){
			if (this->edges[ii].left == ind){
				return ii;
			}
		}
		
		return -11;
	}
	
	
	/**
	 \brief Search this decomposition for any edge (degenerate or not) with input point index as right point.
	 
	 \param ind The index to search for.
	 \return The index of the found edge, or -12 if no edge had the point as right point.
	 */
	int edge_w_right(int ind)
	{
		
		for (int ii=0; ii<num_edges; ii++){
			if (this->edges[ii].right == ind){
				return ii;
			}
		}
		
		return -12;
	}
	
	
	
	/**
	 \brief Search this decomposition for any edge (degenerate or not) with input point index as a removed point.
	 
	 \param ind The index to search for.
	 \return The index of the found edge, or -13 if no edge had the point as a removed point.
	 */
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
	
	
	/**
	 \brief  Find candidates for merging, by finding vertices with type NEW.
	 
	 \return a vector of integers containing the indices of edges to merge together into a single new edge.
	 \param V The vertex set to which the decomposition refers.
	 */
	std::vector<int> get_merge_candidate(const vertex_set & V);
	
	
	/**
	 \brief Merge method for curves, to remove edges with NEW type points.
	 
	 \todo Make this function callable even outside the interslice method.
	 
	 \param W_midpt Skeleton witness set coming from the precedin(containing) interslice call.
	 \param V the vertex set to which this decomposition refers.
	 \param pi_in The linear projection being used to decompose.
	 \param solve_options The current solver configuration.
	 */
	void merge(witness_set & W_midpt, vertex_set & V,
			   vec_mp * pi_in,
			   solver_configuration & solve_options);
	
	
	/**
	 \brief The main call for decomposing an indepentent curve.
	 
	 Includes deflation, etc.
	 
	 \param V The vertex set which runs through everything.
	 \param W The input witness set for the curve, computed probably by Bertini's tracktype:1.
	 \param pi A set of pointers to vec_mp's containing the projection being used.
	 \param program_options The current state of Bertini_real.
	 \param solve_options The current state of the solver routines.
	 */
	void main(vertex_set & V,
			  witness_set & W,
			  vec_mp *pi,
			  BR_configuration & program_options,
			  solver_configuration & solve_options);
	
	
	
	
	
	/**
	\brief the main function for computing the self-intersection of a non-self-conjugate dimension 1 component.
	 
	 \param W		the witness set
	 \param pi	the projection for this decomposition
	 \param V		the vertex set being passed around.  C indexes into here.
	 \param num_vars		the total number of variables in the problem.
	 \param program_options		main structure holding configuration
	 \param	solve_options			structure holding options to pass to a solver.
	 
	 */
	void 	computeCurveNotSelfConj(const witness_set		& W,
									vec_mp					pi,
									vertex_set				&V,
									int						num_vars,
									BR_configuration		& program_options,
									solver_configuration	& solve_options);
	
	
	
	/**
	 \brief the main function for computing cell decom for a curve.  only for use on a self-conjugate component.
	 
	 \param W	witness_set containing linears, patches, points, etc.  much info and calculation performed from this little guy.
	 \param pi the set of projections to use.  in this case, there should be only 1.
	 \param V		vertex set structure into which to place collected data.
	 \param options program configuration.
	 \param solve_options solver configuration.
	 */
	void computeCurveSelfConj(const witness_set & W,
							  vec_mp *pi,
							  vertex_set &V,
							  BR_configuration & options,
							  solver_configuration & solve_options);
	
	
	
	
	/**
	 \brief The main method for slicing above critical points, halfway between them, and connecting the dots.
	 
	 This method peforms three functions.  It slices the curve halfway between the critical points, and tracks the midpoints-upstairs to the bounding critical values, added new points where necessary.  It also calls merge if desired.
	 
	 \return SUCCESSFUL
	 \param W_curve The main witness set for the curve.
	 \param W_crit_real Witness set containing the real critical points for the curve, previously computed.
	 \param pi The projection being used to decompose.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 \param V The vertex set being used to contain the computed points.
	 */
	int interslice(const witness_set & W_curve,
				   const witness_set & W_crit_real,
				   vec_mp *pi,
				   BR_configuration & program_options,
				   solver_configuration & solve_options,
				   vertex_set & V);
	
	
	
	/**
	 \brief Compute actual critical points of the curve, by regeneration and the left-nullspace method.
	 
	 \return SUCCESSFUL
	 \param W_curve Witness set for the curve, from which we will regenerate
	 \param pi the linear projection being used.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 \param W_crit_real The computed value, containing the real critical points of the curve.
	 */
	int compute_critical_points(const witness_set & W_curve,
								vec_mp *pi,
								BR_configuration & program_options,
								solver_configuration & solve_options,
								witness_set & W_crit_real);
	
	
	
	
	
	/**
	 \brief Compute the intersection of the curve with a sphere containing the supplied critical points.
	 
	 Calling the sphere intersection solver, this method finds the points of intersection between the sphere and curve.
	 
	 \return SUCCESSFUL flag.
	 \param W_crit_real The previously computed real critical points.
	 \param W supplied witness set for the curve.
	 \param program_options The current state of Bertini_real
	 \param solve_options The current state of the solver.
	 */
	int get_additional_critpts(witness_set *W_crit_real,
							   const witness_set & W,
							   BR_configuration & program_options,
							   solver_configuration & solve_options);
	
	
	/**
	 \brief send a curve to a single MPI target
	 
	 \param target Who to send to.
	 \param mpi_config The state of MPI.
	 */
	void send(int target, parallelism_config & mpi_config);
	
	
	/**
	 \brief receive a curve from a single MPI source
	 
	 \param source Who to receive from.
	 \param mpi_config The state of MPI.
	 */
	void receive(int source, parallelism_config & mpi_config);
	
	
	
	
	
	
	/**
	 \brief Initialize a curve for sampling by the adaptive method.
	 
	 \see curve_decomposition::adaptive_sampler
	 
	 \ingroup samplermethods
	 
	 \return The number 0.
	 */
	int adaptive_set_initial_sample_data();
	
	/**
	\brief Sample a curve using adaptive method based on distance between computed samples, and a maximum number of refinement passes.
	 
	 \todo Add a better summary of the adaptive cufve method.
	
	 \ingroup samplermethods
	 
	 \param V the vertex set containing the points of the decomposition.
	 \param sampler_options The current state of the sampler program.
	 \param solve_options The current state of the solver and tracker configuration.
	*/
	void adaptive_sampler(vertex_set &V,
						  sampler_configuration & sampler_options,
						  solver_configuration & solve_options);
	
	
	
	/**
	 \brief sets up refinement flags to YES for every interval, for first pass of adaptive refinement.
	 
	 \param num_refinements The number of intervals to refine.
	 \param refine_flags The mutable vector of bools indicating whether to refine a particular interval.
	 \param current_indices Indices of points between which to refine (or not).
	 \param V The vertex set storing the points.
	 \param current_edge The index of which edge is being refined.
	 \param sampler_options The current state of the program.
	 */
	void adaptive_set_initial_refinement_flags(int & num_refinements,
											   std::vector<bool> & refine_flags,
											   std::vector<int> & current_indices,
											   vertex_set &V,
											   int current_edge, sampler_configuration & sampler_options);
	
	
	/**
	 \brief Initialize for the fixed-number curve sampler method.
	 
	 \param target_num_samples The number of samples per edge, including boundary points.
	 \return The number 0.
	 */
	int fixed_set_initial_sample_data(int target_num_samples);
	
	
	/**
	 \brief Sample a curve so it has an equal number of points per edge, including boundary points.
	 
	 \todo Add a description ofthe method with picture right here.
	 
	 \param V the vertex set containing the points computed.
	 \param sampler_options The current state of the sampler program.
	 \param solve_options The current state of the solver.
	 \param target_num_samples The number of points to get on each edge, including boundary points.
	 */
	void fixed_sampler(vertex_set &V,
					   sampler_configuration & sampler_options,
					   solver_configuration & solve_options,
					   int target_num_samples);
	
	
	/**
	 \brief Dump the curve's sampler data to a file.
	 
	 \param samplingName The name of the file to write.
	 */
	void  output_sampling_data(boost::filesystem::path samplingName);
	
	
	
	
	
	
	
	/**
	 \brief Reset the curve to an empty set.
	 */
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








/**
 \brief Check whether a projection is valid or not, by testing the rank of the jacobian at a random point.
 
 \return a boolean integer indicating whether the projection is valid.
 \param W witness set containing required information
 \param projection The linear projection being checked.
 \param solve_options The current state of the solver.
 */
int verify_projection_ok(const witness_set & W,
						 vec_mp * projection,
						 solver_configuration & solve_options);


/**
 \brief Check whether a projection is valid or not, by testing the rank of the jacobian at a random point.
 
 \return a boolean integer indicating whether the projection is valid.
 \param W witness set containing required information
 \param randomizer The way the system is randomized. 
 \param projection The linear projection being checked.
 \param solve_options The current state of the solver.
 */
int verify_projection_ok(const witness_set & W,
						 system_randomizer * randomizer,
						 vec_mp * projection,
						 solver_configuration & solve_options);





#endif
