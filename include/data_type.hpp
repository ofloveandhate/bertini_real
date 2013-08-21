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

enum {SUCCESSFUL=0, CRITICAL_FAILURE=-10, TOLERABLE_FAILURE=-1};



//The following lets us use words instead of numbers to indicate vertex type.
enum {UNSET= 100, CRITICAL, NEW, MIDPOINT, ISOLATED, SAMPLE_POINT};



// enum for worker mode choice
enum {NULLSPACE = 3000, LINPRODTODETJAC, DETJACTODETJAC, LINTOLIN, MULTILIN};

enum {TERMINATE = 2000, INITIAL_STATE};

enum {PARSING = 1000, TYPE_CONFIRMATION, DATA_TRANSMISSION, NUMPACKETS};

enum {INACTIVE = 500, ACTIVE};

enum {VEC_MP = 4000, VEC_D, MAT_MP, MAT_D, COMP_MP, COMP_D, INDICES};

std::string enum_lookup(int flag);

void * br_malloc(size_t size);

void *br_realloc(void *ptr, size_t size);

void br_exit(int errorCode);

void deliberate_segfault();


void print_point_to_screen_matlab(vec_d M, std::string name);
void print_point_to_screen_matlab(vec_mp M, std::string name);
void print_matrix_to_screen_matlab(mat_d M, std::string name);
void print_matrix_to_screen_matlab(mat_mp M, std::string name);

void print_comp_matlab(comp_mp M,std::string name);
void print_comp_matlab(comp_d M,std::string name);









class function
{
	std::string func;  //symbolic representation of function (straight from input file).
										 // this class is woefully incomplete.
};



class witness_set
{
	
public:
	
	//begin data members
	
	vec_mp *L_mp;
	vec_mp *patch_mp;
  point_mp *pts_mp;
	
	int dim;
  int comp_num;
	int incidence_number;
	
	int num_variables;
	int num_synth_vars;
	
	int num_pts;
	int num_linears;
	int num_patches;
	
	std::vector< std::string > variable_names;
	
	boost::filesystem::path input_filename;
	function input_file;
	// end data members
	
	
	
	// overloaded operators
	
	
	// default constructor
	
	witness_set(){
		this->input_filename = "unset_filename";
		
		this->num_variables = 0;
		this->num_synth_vars = 0;
		
		
		this->num_patches = this->num_linears = 0;
		this->num_pts = 0;
		
		this->patch_mp = NULL;
		this->L_mp = NULL;
		this->pts_mp = NULL;
		
		this->incidence_number = -1;
		this->comp_num = this->dim = -1;
	};
	
	
	~witness_set(){ // the destructor
		
		if (this->num_linears>0) {
			for (int ii =0; ii<this->num_linears; ii++) {
				clear_vec_mp(this->L_mp[ii]);
			}
			free(this->L_mp);
		}
		
		if (this->num_patches>0) {
			for (int ii =0; ii<this->num_patches; ii++) {
				clear_vec_mp(this->patch_mp[ii]);
			}
			free(this->patch_mp);
		}
		
		
		if (this->num_pts>0) {
			for (int ii =0; ii<this->num_pts; ii++) {
				clear_vec_mp(this->pts_mp[ii]);
			}
			free(this->pts_mp);
		}
		
		
		this->num_variables = 0;
		this->num_synth_vars = 0;
		
		
		this->num_patches = this->num_linears = 0;
		this->num_pts = 0;
		
		this->patch_mp = this->L_mp = this->pts_mp = NULL;
		
		this->incidence_number = -1;
		this->comp_num = this->dim = -1;
		
	};
	
	
	// assignment
	witness_set& operator=( const witness_set& other) {
		
		// i am pretty sure that this code has leaks /  errors, in that when you assign a non-empty witness set, this attempts to br_malloc, not br_realloc.
		
		this->dim = other.dim;
		this->comp_num = other.comp_num;
		this->incidence_number = other.incidence_number;
		
		
		this->num_variables = other.num_variables;
		this->num_synth_vars = other.num_synth_vars;
		
		this->num_pts = other.num_pts;
		this->num_linears = other.num_linears;
		this->num_patches = other.num_patches;
		
		
		this->variable_names = other.variable_names;
		
		if (this->num_linears>0) {
			this->L_mp = (vec_mp *)br_malloc(other.num_linears*sizeof(vec_mp));
			for (int ii=0; ii<other.num_linears; ii++) {
				init_vec_mp2(this->L_mp[ii],1,1024); this->L_mp[ii]->size = 1;
				vec_cp_mp(this->L_mp[ii], other.L_mp[ii]);
			}
		}
		
		if (this->num_patches>0) {
			this->patch_mp = (vec_mp *)br_malloc(other.num_patches*sizeof(vec_mp));
			for (int ii=0; ii<other.num_patches; ii++) {
				init_vec_mp2(this->patch_mp[ii],1,1024); this->patch_mp[ii]->size = 1;
				vec_cp_mp(this->patch_mp[ii], other.patch_mp[ii]);
			}
		}
		
		if (this->num_pts>0) {
			this->pts_mp = (vec_mp *)br_malloc(other.num_pts*sizeof(vec_mp));
			for (int ii=0; ii<other.num_pts; ii++) {
				init_vec_mp2(this->pts_mp[ii],1,1024); this->pts_mp[ii]->size = 1;
				vec_cp_mp(this->pts_mp[ii], other.pts_mp[ii]);
			}
		}
		
		
    return *this;
  };
	
	
	//copy operator.  must be explicitly declared because the underlying c structures use pointers.
	witness_set(const witness_set & other){
		
		this->patch_mp = NULL;
		this->L_mp = NULL;
		this->pts_mp = NULL;
		
		this->dim = other.dim;
		this->comp_num = other.comp_num;
		this->incidence_number = other.incidence_number;
		
		
		this->num_variables = other.num_variables;
		this->num_synth_vars = other.num_synth_vars;
		
		this->num_pts = other.num_pts;
		this->num_linears = other.num_linears;
		this->num_patches = other.num_patches;
		
		
		this->variable_names = other.variable_names;
		
		if (this->num_linears>0) {
			this->L_mp = (vec_mp *)br_malloc(other.num_linears*sizeof(vec_mp));
			for (int ii=0; ii<other.num_linears; ii++) {
				init_vec_mp2(this->L_mp[ii],1,1024); this->L_mp[ii]->size = 1;
				vec_cp_mp(this->L_mp[ii], other.L_mp[ii]);
			}
		}
		
		if (this->num_patches>0) {
			this->patch_mp = (vec_mp *)br_malloc(other.num_patches*sizeof(vec_mp));
			for (int ii=0; ii<other.num_patches; ii++) {
				init_vec_mp2(this->patch_mp[ii],1,1024); this->patch_mp[ii]->size = 1;
				vec_cp_mp(this->patch_mp[ii], other.patch_mp[ii]);
			}
		}
		
		if (this->num_pts>0) {
			this->pts_mp = (vec_mp *)br_malloc(other.num_pts*sizeof(vec_mp));
			for (int ii=0; ii<other.num_pts; ii++) {
				init_vec_mp2(this->pts_mp[ii],1,1024); this->pts_mp[ii]->size = 1;
				vec_cp_mp(this->pts_mp[ii], other.pts_mp[ii]);
			}
		}
		
	};
	
	
	void only_natural_vars();
	void only_first_vars(int num_vars);
	void sort_for_real(tracker_config_t T);
	void sort_for_unique(tracker_config_t T);
	
	
	int  witnessSetParse(const boost::filesystem::path witness_set_file, const int num_vars);
	
	void reset()
	{
		reset_names();
		
		reset_points();
		
		reset_linears();
		
		reset_patches();
		
		this->input_filename = "unset_filename";
		this->num_variables = 0;
		this->num_synth_vars = 0;
		
		this->incidence_number = -1;
		this->comp_num = this->dim = -1;
	};
	
	
	void reset_names()
	{
		variable_names.clear();
	}
	
	
	void reset_points()
	{
		for (int ii =0; ii<this->num_pts; ii++)
			clear_vec_mp(this->pts_mp[ii]);
		
		if (this->num_pts>0)
			free(this->pts_mp);
		
		
		this->num_pts = 0;
		this->pts_mp = NULL;
	}
	
	void reset_linears()
	{
		for (int ii =0; ii<this->num_linears; ii++)
			clear_vec_mp(this->L_mp[ii]);
		
		if (this->num_linears>0)
			free(this->L_mp);
		
		
		this->num_linears = 0;
		this->L_mp = NULL;
	}
	
	
	void reset_patches()
	{
		for (int ii =0; ii<this->num_patches; ii++)
			clear_vec_mp(this->patch_mp[ii]);
		
		if (this->num_patches>0)
			free(this->patch_mp);
		
		
		this->num_patches = 0;
		this->patch_mp = NULL;
	};
	
	
	
	void add_patch(vec_mp new_patch);
	void add_point(vec_mp new_point);
	void add_linear(vec_mp new_linear);
	
	
	void merge(const witness_set & W_in);///< merges W_in into this
	
	
	void get_variable_names(); ///< reads variable names from names.out
	
	
	void print_to_screen(); ///< prints some information about the witness set to the screen
	
	void print_to_file();
	
	
	/**
	 writes the linears in point form to file filename
	 
	 \param filename the name of the file to be written.
	 */
	void write_linears(boost::filesystem::path filename);
	
	/**
	 writes the patches in point form to file filename
	 
	 \param filename the name of the file to be written.
	 */
	void print_patches(boost::filesystem::path filename);
	void read_patches_from_file(boost::filesystem::path filename);
	
	
	void write_homogeneous_coordinates(boost::filesystem::path filename);
	void write_dehomogenized_coordinates(boost::filesystem::path filename);
	
	void compute_downstairs_crit_midpts(vec_mp crit_downstairs,
																			vec_mp midpoints_downstairs,
																			std::vector< int > & index_tracker,
																			vec_mp pi);
};
// end the double types





void cp_names(witness_set *W_out, witness_set & W_in);
void cp_linears(witness_set *W_out, witness_set & W_in);
void cp_patches(witness_set *W_out, witness_set & W_in);













// CURVE CELL DECOMP DATA TYPES
class vertex
{
public:
	
  point_mp pt_mp;
  comp_mp  projVal_mp;
	
  int type;  //See enum above.
	
	
	vertex()
	{
		init();
	};
	
	~vertex()
	{
		clear_vec_mp(this->pt_mp);
		clear_mp(this->projVal_mp);
	};
	
	vertex & operator=(const vertex & other)
	{
		init();
		
		vec_cp_mp(this->pt_mp, other.pt_mp);
		set_mp(this->projVal_mp, other.projVal_mp);
		this->type = other.type;
		return *this;
	};
	
	vertex(const vertex& other)
	{
		init();
		
		vec_cp_mp(this->pt_mp, other.pt_mp);
		set_mp(this->projVal_mp, other.projVal_mp);
		this->type = other.type;
	}
	
	void print()
	{
		print_point_to_screen_matlab(pt_mp,"point");
		print_comp_matlab(projVal_mp,"projVal");
		std::cout << "type: " << type << std::endl;
	}
	
private:
	
	void init()
	{
		init_mp2(this->projVal_mp,1024);
		init_point_mp2(this->pt_mp,1,1024);
		this->pt_mp->size = 1;
		this->type = UNSET;
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
		
		this->vertices = other.vertices;
		
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
	
	edge()
	{
		left = right = midpt = -1;
	}
	
	edge(int left_, int midpt_, int right_)
	{
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
	int component_num;
	
	int num_curr_projections;
	vec_mp	*pi; // the projections
	
	
	mat_mp randomizer_matrix;
	
	int num_patches;
	vec_mp *patch;

	boost::filesystem::path input_filename;
//	function input_file;
	
	decomposition(){
		init();
	}
	
	~decomposition()
	{
		if (num_curr_projections>0){
		for (int ii=0; ii<num_curr_projections; ii++) 
			clear_vec_mp(pi[ii]);
		free(pi);
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
		
		init();
		
		copy(other);
		
		return *this;
	}
	
	decomposition(const decomposition & other){
		
		init();
		
		copy(other);
	}
	
	

	void add_projection(vec_mp proj){
		if (this->num_curr_projections==0) {
			pi = (vec_mp *) br_malloc(sizeof(vec_mp));
		}
		else{
			this->pi = (vec_mp *)br_realloc(this->pi, (this->num_curr_projections+1) * sizeof(vec_mp));
		}
		
		init_vec_mp(this->pi[num_curr_projections],proj->size);
		this->pi[num_curr_projections]->size = proj->size;
		
		vec_cp_mp(pi[num_curr_projections], proj);
		num_curr_projections++;
	}
	
	void add_patch(vec_mp new_patch){
		if (this->num_patches==0) {
			this->patch = (vec_mp *) br_malloc(sizeof(vec_mp));
		}
		else{
			this->patch = (vec_mp *)br_realloc(this->patch, (this->num_patches+1) * sizeof(vec_mp));
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
	
	virtual void print(boost::filesystem::path outputfile);
	
	void init(){
		pi = NULL;
		patch = NULL;
		init_mat_mp(randomizer_matrix, 0, 0);
		randomizer_matrix->rows = randomizer_matrix->cols = 0;
		
		num_curr_projections = num_patches = 0;
		num_variables = 0;
		dimension = -1;
		component_num = -1;
	}
	
	void copy(const decomposition & other)
	{
		
		
		
		this->input_filename = other.input_filename;
//		this->input_file = other.input_file;
		
		this->counters = other.counters;
		this->indices = other.indices;
		this->num_variables = other.num_variables;
		this->dimension = other.dimension;
		this->component_num = other.component_num;
		
		this->num_curr_projections = other.num_curr_projections;
		this->pi = (vec_mp *) br_malloc(other.num_curr_projections * sizeof(vec_mp));
		for (int ii = 0; ii<other.num_curr_projections; ii++) {
			init_vec_mp(this->pi[ii],other.pi[ii]->size);
			this->pi[ii]->size = other.pi[ii]->size;
			vec_cp_mp(this->pi[ii], other.pi[ii])
		}
		
		
		init_mat_mp2(this->randomizer_matrix, 0, 0,1024); randomizer_matrix->rows = randomizer_matrix->cols = 0;
		mat_cp_mp(this->randomizer_matrix, other.randomizer_matrix);
		
		
		
		this->num_patches = other.num_patches;
		this->patch = (vec_mp *) br_malloc(other.num_patches * sizeof(vec_mp));
		for (int ii = 0; ii<other.num_patches; ii++) {
			init_vec_mp(this->patch[ii],other.patch[ii]->size);
			this->patch[ii]->size = other.patch[ii]->size;
			vec_cp_mp(this->patch[ii], other.patch[ii])
		}
	}
	

	
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
	
	std::vector< curve_decomposition > midpoint_slices;
	std::vector< curve_decomposition > critpoint_slices;
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
	
	void init()
	{
		decomposition::init();
		num_edges = 0;
		num_faces = 0;
		dimension = 2;
	}
	
	void copy(const surface_decomposition & other)
	{
		this->faces = other.faces;
		this->edges = other.edges;
		
		this->num_edges = other.num_edges;
		this->num_faces = other.num_faces;
		
		this->midpoint_slices = other.midpoint_slices;
		this->critpoint_slices = other.critpoint_slices;
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





void print_path_retVal_message(int retVal);

/**
retrieves the number of variables from the PPD by taking the sum of the sizes, plus the sum of the types.
*/
int get_num_vars_PPD(preproc_data PPD);


void cp_patch_mp(patch_eval_data_mp *PED, patch_eval_data_mp PED_input);
void cp_patch_d(patch_eval_data_d *PED, patch_eval_data_d PED_input);
void cp_preproc_data(preproc_data *PPD, preproc_data PPD_input);

void sort_increasing_by_real(vec_mp *projections_sorted, std::vector< int > & index_tracker, vec_mp projections_input);

void make_randomization_matrix_based_on_degrees(mat_mp randomization_matrix, int ** randomized_degrees,
																								int num_variables, int num_funcs);
int compare_integers_decreasing(const void * left_in, const void * right_in);

void send_patch_mp   (patch_eval_data_mp * patch);
void receive_patch_mp(patch_eval_data_mp * patch);


void send_patch_d   (patch_eval_data_d * patch);
void receive_patch_d(patch_eval_data_d * patch);


void send_preproc_data(preproc_data *PPD);
void receive_preproc_data(preproc_data *PPD);


void send_vec_mp(vec_mp b, int target);
void receive_vec_mp(vec_mp b, int source);

void send_vec_d(vec_d b, int target);
void receive_vec_d(vec_d b, int source);
#endif

