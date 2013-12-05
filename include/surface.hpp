
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


#include <boost/filesystem.hpp>


#ifndef SURFACE_H
#define SURFACE_H

#include "missing_bertini_headers.hpp"


#include "fileops.hpp"


#include "nullspace_left.hpp"

#include "solver_midpoint_tracker.hpp"

#include "checkSelfConjugate.hpp"
#include "curve.hpp"
#include "data_type.hpp"
#include "isosingular.hpp"
#include "programConfiguration.hpp"

#include "solver_sphere.hpp"



/**
 the face data type..
 */
class face : public cell
{
public:
	
  std::vector<int>	left;  ///< index into vertices
  std::vector<int>	right; ///< index into vertices
	
	int top;				///<  index into edges.  indexes into crit_curve.edges
	int bottom;			///<  index into edges
	
	int num_left;		///<  counters
	int num_right;	///<
	
	
	comp_mp left_crit_val; ///<
	comp_mp right_crit_val; ///<
	
	int system_type_top;
	int system_type_bottom;
	
	int crit_slice_index; ///< which midpoint this face came from.
	
	face() : cell()
	{
		init();
	} ///<  constructor
	
	~face(){
		clear_mp(left_crit_val);
		clear_mp(right_crit_val);
	}  ///<  destructor
	
	face(const face & other){ ///<  copy
		init();
		copy(other);
	}
	
	face& operator=(const face & other){ ///<  assignment

		init();
		
		copy(other);
		return *this;
	}
	
	
	void send(int target, parallelism_config & mpi_config);
	
	void receive(int source, parallelism_config & mpi_config);
	
	
	
	
private:
	
	
	
	void init()
	{
		system_type_top = UNSET;
		system_type_bottom = UNSET;
		
		
		init_mp2(left_crit_val,1024);
		init_mp2(right_crit_val,1024);
		
		num_left = num_right = 0;
		top = bottom = -1;
		crit_slice_index = -1;
	}
	
	
	void copy(const face & other)
	{
		cell::copy(other);
		
		this->system_type_bottom = other.system_type_bottom;
		this->system_type_top = other.system_type_top;
		
		
		this->crit_slice_index = other.crit_slice_index;
		
		this->left = other.left;
		this->right = other.right;
		
		this->top = other.top;
		this->bottom = other.bottom;
		
		this->num_left = other.num_left;
		this->num_right = other.num_right;
		
		set_mp(this->left_crit_val, other.left_crit_val);
		set_mp(this->right_crit_val, other.right_crit_val);
	}
	
};







/**
 surface decomposition.
 
 includes methods to add vertices, look up vertices, etc
 */
class surface_decomposition : public decomposition
{
	
	std::vector<edge> edges;
	std::vector<face> faces;
	//these counters keep track of the number of things
	
	int      num_edges;
	int      num_faces;
	
	
	
public:
	
	std::vector< curve_decomposition > mid_slices;
	std::vector< curve_decomposition > crit_slices;
	curve_decomposition crit_curve;
	curve_decomposition sphere_curve;
	

	
	

	
	
	void print(boost::filesystem::path base);
	
	void print_faces(boost::filesystem::path outputfile);
	
	

	

	
	void add_face(const face & F)
	{
		faces.push_back(F);
		num_faces++;
	}
	
	void add_edge(const edge & E)
	{
		edges.push_back(E);
		num_edges++;
	}
	
	

	
	
	void main(vertex_set & V,
						const witness_set & W,
						vec_mp *pi,
						BR_configuration & program_options,
						solver_configuration & solve_options);
	
	
	void beginning_stuff(const witness_set & W_surf,
																							BR_configuration & program_options,
																							solver_configuration & solve_options);
	
	
	void compute_slices(const witness_set W_surf,
											vertex_set & V,
											vec_mp projection_values_downstairs,
											std::vector< curve_decomposition > & slices,
											BR_configuration & program_options,
											solver_configuration & solve_options,
											bool rerun_empty,
											std::string kindofslice);
	
	
    
    void compute_critcurve_witness_set(witness_set & W_critcurve,
                                                              const witness_set & W_surf,
                                                              BR_configuration & program_options,
                                       solver_configuration & solve_options);
    
    void compute_critcurve_critpts(witness_set & W_critcurve_crit, // the computed value
                                   const witness_set & W_surf, // input witness set
                                   const witness_set & W_critcurve,
                                   BR_configuration & program_options,
                                   solver_configuration & solve_options);
    
    
	
	void compute_critical_curve(const witness_set & W_critcurve,
                                const witness_set & W_critpts,
                                vertex_set & V,
                                BR_configuration & program_options,
                                solver_configuration & solve_options);
	
	
	

	
	
	
	void compute_sphere_witness_set(const witness_set & W_surf,
																	witness_set & W_intersection_sphere,//
																	BR_configuration & program_options,//
																	solver_configuration & solve_options);
	
	void compute_sphere_crit(const witness_set & W_intersection_sphere,
													 witness_set & W_sphere_crit,
													 BR_configuration & program_options,
													 solver_configuration & solve_options);
	

	
	void compute_bounding_sphere(const witness_set & W_intersection_sphere,
															 const witness_set & W_crit,
															 vertex_set & V,
															 BR_configuration & program_options,
															solver_configuration & solve_options);
	
	
	void connect_the_dots(vertex_set & V,
												 vec_mp *pi,
												 BR_configuration & program_options,
												 solver_configuration & solve_options);
	
	void serial_connect(vertex_set & V, midpoint_config & md_config, solver_configuration & solve_options, BR_configuration & program_options);
	
	void master_connect(vertex_set & V, midpoint_config & md_config, solver_configuration & solve_options, BR_configuration & program_options);
	
	void worker_connect(solver_configuration & solve_options, BR_configuration & program_options);
	
	void master_face_requester(int ii, int jj, int next_worker, parallelism_config & mpi_config);
	
	void worker_face_requester(int & ii, int & jj, parallelism_config & mpi_config);
	
	face make_face(int ii, int jj, vertex_set & V,
				   midpoint_config & md_config,
				   solver_configuration & solve_options, BR_configuration & program_options);
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	surface_decomposition() : decomposition()
	{
		this->init();
	}
	
	surface_decomposition & operator=(const surface_decomposition& other){
		this->init();
		this->copy(other);
		return *this;
	}
	
	surface_decomposition(const surface_decomposition & other){
		this->init();
		this->copy(other);
	}
	
	
	~surface_decomposition()
	{
		this->clear();
	}
	
	
	
	
	
	
	void send(int target, parallelism_config & mpi_config);
	
	
	void receive(int source, parallelism_config & mpi_config);
	
protected:
	
	
	void copy(const surface_decomposition & other)
	{
		decomposition::copy(other);
		
		this->faces = other.faces;
		this->edges = other.edges;
		
		this->num_edges = other.num_edges;
		this->num_faces = other.num_faces;
		
		this->mid_slices = other.mid_slices;
		this->crit_slices = other.crit_slices;
		this->crit_curve = other.crit_curve;
	}
	
	
	void init()
	{

		num_edges = 0;
		num_faces = 0;
		dimension = 2;
		
		
	}
	
	void clear()
	{
		
	}
};








void create_sliced_system(boost::filesystem::path input_file, boost::filesystem::path output_file,
													vec_mp * linears, int num_to_add,
													const witness_set & W);


void create_sphere_system(boost::filesystem::path input_file, boost::filesystem::path output_file,
													comp_mp sphere_diameter,
													vec_mp sphere_center,
													const witness_set & W);





#endif


