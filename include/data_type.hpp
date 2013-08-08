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
#include <signal.h>


#include <iostream>
#include <ios>
#include <string>
#include <fstream>

#include <vector>
#include <sstream>
#include <map>

#ifndef _DATA_TYPE_H
#define _DATA_TYPE_H

#include "missing_bertini_headers.hpp"


#include "boost/filesystem.hpp"
#include "fileops.hpp"



/*** low-level data types. ***/


//The following lets us use words instead of numbers to indicate vertex type.
enum {UNSET=-10, CRITICAL=0, NEW=1, MIDPOINT=2, ISOLATED=-1, SAMPLE_POINT=3};

enum {SUCCESSFUL=0, CRITICAL_FAILURE=-10, TOLERABLE_FAILURE=-1};




void * br_malloc(size_t size);

void br_exit(int errorCode);

void deliberate_segfault();






class function
{
	std::string func;  //symbolic representation of function (straight from input file).
										 // this class is woefully incomplete.
};







// CURVE CELL DECOMP DATA TYPES
class vertex
{
public:
	
  point_mp pt_mp;
  comp_mp  projVal_mp;
	
  int type;  //See enum above.
	
	
	vertex()
	{
		init_mp(this->projVal_mp);
		init_point_mp(this->pt_mp,1);
		this->pt_mp->size = 1;
		this->type = UNSET;
	};
	
	~vertex()
	{
		clear_vec_mp(this->pt_mp);
		clear_mp(this->projVal_mp);
	};
	
	vertex & operator=(const vertex & other)
	{
		
		init_mp(this->projVal_mp);
		init_point_mp(this->pt_mp,1); this->pt_mp->size = 1;
		
		vec_cp_mp(this->pt_mp, other.pt_mp);
		set_mp(this->projVal_mp, other.projVal_mp);
		this->type = other.type;
		return *this;
	};
	
	vertex(const vertex& other)
	{
		init_mp(this->projVal_mp);
		init_point_mp(this->pt_mp,1);
		this->pt_mp->size = 1;
		
		vec_cp_mp(this->pt_mp, other.pt_mp);
		set_mp(this->projVal_mp, other.projVal_mp);
		this->type = other.type;
	}
};



/**
 the main structure for storing vertices.  
 there are methods in place to add vertices, and perform lookups.
 */
class vertex_set
{
public: 
	std::vector<vertex> vertices;  //Isolated real points.
	int num_vertices;
	
	int num_natural_variables;  ///< the number of natural variables appearing in the problem to solve.
	
	void print_to_screen(); ///< operator for displaying information to screen
	
	int add_vertex(const vertex new_vertex);
	int setup_vertices(boost::filesystem::path INfile);
	
	mpf_t abs;
	mpf_t zerothresh;
	comp_mp diff;
	vec_mp checker_1;
	vec_mp checker_2;
	
	vertex_set(){
		init();
	}
	
	vertex_set(int num_vars){
		init();
		
		this->num_natural_variables = num_vars;
		
		change_size_vec_mp(checker_1, num_vars);
		change_size_vec_mp(checker_2, num_vars);
		checker_1->size = checker_2->size = num_vars;
	}
	

	

	
	
	vertex_set & operator=( const vertex_set& other) {
		copy(other);
		return *this;
	}
	
	vertex_set(const vertex_set &other)
	{
		init();
		copy(other);
	}
	
	~vertex_set()
	{
		clear();
	}
	
	
	void print(boost::filesystem::path outputfile);
	
	
private:
	
	void init()
	{
		this->num_vertices = 0;
		this->num_natural_variables = 0;
		
		init_vec_mp(checker_1,0);
		init_vec_mp(checker_2,0);
		
		
		
		init_mp(this->diff);

		mpf_init(abs);
		mpf_init(zerothresh);
		mpf_set_d(zerothresh, 1e-8);
	}
	
	
	void copy(const vertex_set &other)
	{
		this->num_vertices = other.num_vertices;
		this->num_natural_variables = other.num_natural_variables;
		
		vec_cp_mp(this->checker_1,other.checker_1);
		vec_cp_mp(this->checker_2,other.checker_2);
		
	}
	
	void clear()
	{
		clear_vec_mp(checker_1);
		clear_vec_mp(checker_2);
	}

};




class cell
{
	
private:
//	int n;
//	function homotopy;
	
public:
	
};



/**
 the edge data type.  has three indices: left, right, midpt.
 */
class edge : public cell
{
public:
  int left;  ///< index into vertices
  int right; ///< index into vertices
	int midpt; ///<  index into vertices
	
	edge() {
		left = right = midpt = -1;
	}
	
	edge(int left_, int midpt_, int right_){
		this->left = left_;
		this->right = right_;
		this->midpt = midpt_;
	}
	
	// other defaults are correct for this type.
	
};




/**
 the face data type..
 */
class face : public cell
{
public:
	
  std::vector<int>	left;  ///< index into vertices
  std::vector<int>	right; ///< index into vertices
	int top; ///<  index into edges
	int bottom; ///<  index into edges
	
	int num_left;  ///<  counters
	int num_right; ///< 
	
	comp_mp left_crit_val; ///< 
	comp_mp right_crit_val; ///< 
	
	int interior_pt; ///< index into vertex set
	
	face(){init_mp(left_crit_val); init_mp(right_crit_val);} ///<  constructor
	~face(){clear_mp(left_crit_val); clear_mp(right_crit_val);}  ///<  destructor
	
	face(const face & other){ ///<  copy
		init_mp(left_crit_val); init_mp(right_crit_val);
		
		set_mp(this->left_crit_val, other.left_crit_val);
		set_mp(this->right_crit_val, other.right_crit_val);
	}
	
	face& operator=(const face & other){ ///<  assignment
		init_mp(left_crit_val); init_mp(right_crit_val);
		
		set_mp(this->left_crit_val, other.left_crit_val);
		set_mp(this->right_crit_val, other.right_crit_val);
		return *this;
	}
	
	
};







class decomposition
{

public:
	std::map< int , int > counters;
	std::map< int , std::vector< int > > indices;
	
	int num_variables;
	int dimension;
	
	int num_curr_projections;
	vec_mp	*pi_mp; // the projections
	
	
	mat_mp randomizer_matrix;
	
	int num_patches;
	vec_mp *patch;

	
	
	decomposition(){
		pi_mp = NULL;
		patch = NULL;
		init_mat_mp(randomizer_matrix, 0, 0);
		randomizer_matrix->rows = randomizer_matrix->cols = 0;
		
		num_curr_projections = num_patches = 0;
		num_variables = 0;
		dimension = -1;
	}
	
	~decomposition()
	{
		if (num_curr_projections>0){
		for (int ii=0; ii<num_curr_projections; ii++) 
			clear_vec_mp(pi_mp[ii]);
		free(pi_mp);
		}
		
		if (num_patches>0){
			for (int ii=0; ii<num_patches; ii++)
				clear_vec_mp(patch[ii]);
			free(patch);
		}
		
		counters.clear();
		indices.clear();
		
		clear_mat_mp(randomizer_matrix);
	}
	
	decomposition & operator=(const decomposition& other){
		
		this->counters = other.counters;
		this->indices = other.indices;
		this->num_variables = other.num_variables;
		this->dimension = other.dimension;
		this->num_curr_projections = other.num_curr_projections;
		init_mat_mp(this->randomizer_matrix, 0, 0); randomizer_matrix->rows = randomizer_matrix->cols = 0;
		mat_cp_mp(this->randomizer_matrix, other.randomizer_matrix);
		
		this->pi_mp = (vec_mp *) br_malloc(other.dimension * sizeof(vec_mp));
		for (int ii = 0; ii<other.dimension; ii++) {
			init_vec_mp(this->pi_mp[ii],other.pi_mp[ii]->size);
			this->pi_mp[ii]->size = other.pi_mp[ii]->size;
			vec_cp_mp(this->pi_mp[ii], other.pi_mp[ii])
		}
		
		this->num_patches = other.num_patches;
		this->patch = (vec_mp *) br_malloc(other.num_patches * sizeof(vec_mp));
		for (int ii = 0; ii<other.num_patches; ii++) {
			init_vec_mp(this->patch[ii],other.patch[ii]->size);
			this->patch[ii]->size = other.patch[ii]->size;
			vec_cp_mp(this->patch[ii], other.patch[ii])
		}
		return *this;
	}
	
	decomposition(const decomposition & other){
		this->counters = other.counters;
		this->indices = other.indices;
		this->num_variables = other.num_variables;
		this->dimension = other.dimension;
		this->num_curr_projections = other.num_curr_projections;
		init_mat_mp(this->randomizer_matrix, 0, 0); randomizer_matrix->rows = randomizer_matrix->cols = 0;
		mat_cp_mp(this->randomizer_matrix, other.randomizer_matrix);
		
		this->pi_mp = (vec_mp *) br_malloc(other.dimension * sizeof(vec_mp));
		for (int ii = 0; ii<other.dimension; ii++) {
			init_vec_mp(this->pi_mp[ii],other.pi_mp[ii]->size);
			this->pi_mp[ii]->size = other.pi_mp[ii]->size;
			vec_cp_mp(this->pi_mp[ii], other.pi_mp[ii])
		}
		
		this->num_patches = other.num_patches;
		this->patch = (vec_mp *) br_malloc(other.num_patches * sizeof(vec_mp));
		for (int ii = 0; ii<other.num_patches; ii++) {
			init_vec_mp(this->patch[ii],other.patch[ii]->size);
			this->patch[ii]->size = other.patch[ii]->size;
			vec_cp_mp(this->patch[ii], other.patch[ii])
		}
	}
	
	void add_projection(vec_mp proj){
		if (this->num_curr_projections==0) {
			pi_mp = (vec_mp *) br_malloc(sizeof(vec_mp));
		}
		else{
			this->pi_mp = (vec_mp *)brealloc(this->pi_mp, (this->num_curr_projections+1) * sizeof(vec_mp));
		}
		
		init_vec_mp(this->pi_mp[num_curr_projections],proj->size);
		this->pi_mp[num_curr_projections]->size = proj->size;
		
		vec_cp_mp(pi_mp[num_curr_projections], proj);
		num_curr_projections++;
	}
	
	
	
	void add_patch(vec_mp new_patch){
		if (this->num_patches==0) {
			this->patch = (vec_mp *) br_malloc(sizeof(vec_mp));
		}
		else{
			this->patch = (vec_mp *)brealloc(this->patch, (this->num_patches+1) * sizeof(vec_mp));
		}
		
		init_vec_mp(this->patch[num_patches],new_patch->size);
		this->patch[this->num_patches]->size = new_patch->size;
		
		vec_cp_mp(this->patch[this->num_patches], new_patch);
		this->num_patches++;
		std::cout << "adding " << this->num_patches << "th patch to decomp\n";
	}
	
	
	int add_vertex(vertex_set &V, vertex source_vertex);
	
	int index_in_vertices(vertex_set &V,
												vec_mp testpoint, comp_mp projection_value,
												tracker_config_t T);
	
	int index_in_vertices_with_add(vertex_set &V,
																 vec_mp testpoint, comp_mp projection_value,
																 tracker_config_t T);
	
	int setup(boost::filesystem::path INfile,
						boost::filesystem::path & inputName,
						boost::filesystem::path directoryName);
	
	void print(boost::filesystem::path input_deflated_Name, boost::filesystem::path outputfile);
}; // end decomposition




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
	
	curve_decomposition()
	{
		num_edges = 0;
		dimension = 1;
	}
}; // end curve_decomposition






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
	surface_decomposition()
	{
		num_edges = 0;
		num_faces = 0;
	}
	
	curve_decomposition crit_curve;
};






typedef struct
{
	int num_variables;
	int num_edges;
	
	int *num_samples_each_edge;
	
	int **sample_indices;
}
sample_data;





void clear_sample(sample_data *S, int MPType);










//function prototypes for bertini_real data clearing etc.

void norm_of_difference(mpf_t result, vec_mp left, vec_mp right);

void dehomogenize(point_d *result, point_d dehom_me);
void dehomogenize(point_d *result, point_d dehom_me, int num_variables);

void dehomogenize(point_mp *result, point_mp dehom_me);
void dehomogenize(point_mp *result, point_mp dehom_me, int num_variables);


void dot_product_d(comp_d result, vec_d one, vec_d two);
void dot_product_mp(comp_mp result, vec_mp one, vec_mp two);

void projection_value_homogeneous_input(comp_d result, vec_d input, vec_d projection);
void projection_value_homogeneous_input(comp_mp result, vec_mp input, vec_mp projection);


int isSamePoint_inhomogeneous_input(point_d left, point_d right);
int isSamePoint_inhomogeneous_input(point_mp left, point_mp right);



int isSamePoint_homogeneous_input(point_d left, point_d right);

int isSamePoint_homogeneous_input(point_mp left, point_mp right);



void print_point_to_screen_matlab(vec_d M, std::string name);
void print_point_to_screen_matlab(vec_mp M, std::string name);
void print_matrix_to_screen_matlab(mat_d M, std::string name);
void print_matrix_to_screen_matlab(mat_mp M, std::string name);

void print_comp_matlab(comp_mp M,std::string name);
void print_comp_matlab(comp_d M,std::string name);

void print_path_retVal_message(int retVal);

/**
retrieves the number of variables from the PPD by taking the sum of the sizes, plus the sum of the types.
*/
int get_num_vars_PPD(preproc_data PPD);


void cp_patch_mp(patch_eval_data_mp *PED, patch_eval_data_mp PED_input);
void cp_patch_d(patch_eval_data_d *PED, patch_eval_data_d PED_input);
void cp_preproc_data(preproc_data *PPD, preproc_data PPD_input);

void sort_increasing_by_real(vec_mp *projections_sorted, int **index_tracker, vec_mp projections_input);

void make_randomization_matrix_based_on_degrees(mat_mp randomization_matrix, int ** randomized_degrees,
																								int num_variables, int num_funcs);
int compare_integers_decreasing(const void * left_in, const void * right_in);

void send_patch_mp   (patch_eval_data_mp * patch);
void receive_patch_mp(patch_eval_data_mp * patch);


void send_patch_d   (patch_eval_data_d * patch);
void receive_patch_d(patch_eval_data_d * patch);


void send_preproc_data(preproc_data *PPD);
void receive_preproc_data(preproc_data *PPD);
#endif

