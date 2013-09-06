
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
	
	int index; ///< which midpoint this face came from.
	
	face() : cell()
	{
		init();
	} ///<  constructor
	
	~face(){clear_mp(left_crit_val); clear_mp(right_crit_val);}  ///<  destructor
	
	face(const face & other){ ///<  copy
		init();
		copy(other);
	}
	
	face& operator=(const face & other){ ///<  assignment

		init();
		
		copy(other);
		return *this;
	}
	
	void init()
	{
		init_mp(left_crit_val);
		init_mp(right_crit_val);
		
		num_left = num_right = 0;
		top = bottom = -1;
		index = -1;
	}
	
	
	void copy(const face & other)
	{
		cell::copy(other);
		
		this->index = other.index;
		
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
	
	
	surface_decomposition() : decomposition()
	{
		init();
	}
	
	surface_decomposition & operator=(const surface_decomposition& other){
		init();
		copy(other);
		return *this;
	}
	
	surface_decomposition(const surface_decomposition & other){
		init();
		copy(other);
	}
	
	void print(boost::filesystem::path base);
	
	void print_faces(boost::filesystem::path outputfile);
	
	
	void init()
	{
		decomposition::init();
		num_edges = 0;
		num_faces = 0;
		dimension = 2;
	}
	
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
	
	
	
	void compute_critical_curve(witness_set & W_curve_crit,
															const witness_set & W_surf,
															vertex_set & V,
															BR_configuration & program_options,
															solver_configuration & solve_options);
	
	
	
	void connect_the_dots(vertex_set & V,
												vec_mp *pi,
												BR_configuration & program_options,
												solver_configuration & solve_options);
	
	
};








void create_sliced_system(boost::filesystem::path input_file, boost::filesystem::path output_file,
													vec_mp * linears, int num_to_add,
													const witness_set & W);






#endif


