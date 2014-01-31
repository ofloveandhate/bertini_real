#ifndef _DATA_TYPE_H
#define _DATA_TYPE_H

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

#include <set>
#include <iostream>
#include <ios>
#include <string>
#include <fstream>

#include <vector>
#include <sstream>
#include <map>

#define DEFAULT_MAX_PREC 1024

#include "missing_bertini_headers.hpp"

#include <boost/timer/timer.hpp>
#include "boost/filesystem.hpp"

#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/bind.hpp>

#include "fileops.hpp"

#define SAMEPOINTTOL 1e-6

class parallelism_config; // a forward declaration
class BR_configuration;   // a forward declaration

/*** low-level data types. ***/

enum {SUCCESSFUL=0, CRITICAL_FAILURE=-10, TOLERABLE_FAILURE=-1};


enum {SYSTEM_CRIT = -1600, SYSTEM_SPHERE};

//The following lets us use words instead of numbers to indicate vertex type.
enum {UNSET= 100, CRITICAL, SEMICRITICAL, MIDPOINT, ISOLATED, NEW, SAMPLE_POINT, REMOVED, PROBLEMATIC};



// enum for worker mode choice
enum {NULLSPACE = 3000, LINPRODTODETJAC, DETJACTODETJAC, LINTOLIN, MULTILIN, MIDPOINT_SOLVER, SPHERE_SOLVER};

enum {TERMINATE = 2000, INITIAL_STATE};

enum {PARSING = 1000, TYPE_CONFIRMATION, DATA_TRANSMISSION, NUMPACKETS};

enum {INACTIVE = 500, ACTIVE};

enum {VEC_MP = 4000, VEC_D, MAT_MP, MAT_D, COMP_MP, COMP_D, VEC_RAT, MAT_RAT, COMP_RAT, INDICES,
	DECOMPOSITION, CURVE, SURFACE, EDGE, CELL, FACE, UNUSED};

std::string enum_lookup(int flag);

void *br_malloc(size_t size);

void *br_realloc(void *ptr, size_t size);

void br_exit(int errorCode);

void deliberate_segfault();




template <typename key_type, typename value_type>
value_type map_lookup_with_default(const  std::map <key_type,value_type> & mc_mapperson, const key_type & lookup_key, const value_type& default_value ) {
	typename std::map<key_type,value_type>::const_iterator it = mc_mapperson.find( lookup_key );
	if ( it == mc_mapperson.end() ) {
		return default_value;
	}
	else {
		return it->second;
	}
}















bool is_identity(mat_d M);
bool is_identity(mat_mp M);

//function prototypes for bertini_real data clearing etc.

void norm_of_difference(mpf_t result, vec_mp left, vec_mp right);

void dehomogenize(point_d *result, point_d dehom_me);
void dehomogenize(point_d *result, point_d dehom_me, int num_variables);

void dehomogenize(point_mp *result, point_mp dehom_me);
void dehomogenize(point_mp *result, point_mp dehom_me, int num_variables);


void dot_product_d(comp_d result, vec_d one, vec_d two);
void dot_product_mp(comp_mp result, vec_mp one, vec_mp two);

void dot_product_mindim(comp_d result, vec_d left, vec_d right);
void dot_product_mindim(comp_mp result, vec_mp left, vec_mp right);

void projection_value_homogeneous_input(comp_d result, vec_d input, vec_d projection);
void projection_value_homogeneous_input(comp_mp result, vec_mp input, vec_mp projection);


int isSamePoint_inhomogeneous_input(point_d left, point_d right);
int isSamePoint_inhomogeneous_input(point_mp left, point_mp right);



int isSamePoint_homogeneous_input(point_d left, point_d right);
int isSamePoint_homogeneous_input(point_mp left, point_mp right);


void real_threshold(comp_mp blabla, double threshold);
void real_threshold(vec_mp blabla, double threshold);
void real_threshold(mat_mp blabla, double threshold);


void print_path_retVal_message(int retVal);

/**
 retrieves the number of variables from the PPD by taking the sum of the sizes, plus the sum of the types.
 */
int get_num_vars_PPD(preproc_data PPD);


void cp_patch_mp(patch_eval_data_mp *PED, patch_eval_data_mp PED_input);
void cp_patch_d(patch_eval_data_d *PED, patch_eval_data_d PED_input);
void cp_preproc_data(preproc_data *PPD, const preproc_data & PPD_input);

void clear_post_process_t(post_process_t * endPoint, int num_vars);


void print_tracker(const tracker_config_t * T);


int sort_increasing_by_real(vec_mp projections_sorted, std::vector< int > & index_tracker, vec_mp projections_input);

void make_randomization_matrix_based_on_degrees(mat_mp randomization_matrix, std::vector< int > & randomized_degrees,
																								int num_variables, int num_funcs);
int compare_integers_decreasing(const void * left_in, const void * right_in);
int compare_integers_increasing(const void * left_in, const void * right_in);

void send_patch_mp   (patch_eval_data_mp * patch);
void receive_patch_mp(patch_eval_data_mp * patch);


void send_patch_d   (patch_eval_data_d * patch);
void receive_patch_d(patch_eval_data_d * patch);


void send_preproc_data(preproc_data *PPD);
void receive_preproc_data(preproc_data *PPD);

//
//void send_vec_mp(vec_mp b, int target);
//void receive_vec_mp(vec_mp b, int source);
//
//void send_vec_d(vec_d b, int target);
//void receive_vec_d(vec_d b, int source);



void send_mat_d(mat_d A, int target);
void receive_mat_d(mat_d A, int source);



void send_mat_mp(mat_mp A, int target);
void receive_mat_mp(mat_mp A, int source);

void send_mat_rat(mat_d A_d, mat_mp A_mp, mpq_t ***A_rat, int target);
void receive_mat_rat(mat_d A_d, mat_mp A_mp, mpq_t ***A_rat, int source);

void send_vec_d(vec_d b, int target);
void receive_vec_d(vec_d b, int source);


void send_vec_mp(vec_mp b, int target);
void receive_vec_mp(vec_mp b, int source);

void send_vec_rat(mpq_t ***b, int size, int target);
void receive_vec_rat(mpq_t ***b, int size, int source);


void send_comp_d(comp_d c, int target);
void receive_comp_d(comp_d c, int source);


void send_comp_num_d(comp_d *c, int num, int target);
void receive_comp_num_d(comp_d *c, int num, int source);


void send_comp_mp(comp_mp c, int target);
void receive_comp_mp(comp_mp c, int source);


void send_comp_num_mp(comp_mp *c, int num, int target);
void receive_comp_num_mp(comp_mp *c, int num, int source);

void send_comp_rat(mpq_t c[2], int target);
void receive_comp_rat(mpq_t c[2], int source);

void send_comp_num_rat(mpq_t c[][2], int num, int target);
void receive_comp_num_rat(mpq_t c[][2], int num, int source);






void print_point_to_screen_matlab(const vec_d M, std::string name);
void print_point_to_screen_matlab(const vec_mp M, std::string name);
void print_matrix_to_screen_matlab(const mat_d M, std::string name);
void print_matrix_to_screen_matlab(const mat_mp M, std::string name);

void print_comp_matlab(const comp_mp M,std::string name);
void print_comp_matlab(const comp_d M,std::string name);





class mpi_data_class
{
	
};


namespace bertini {
	class complex : public mpi_data_class
	{
		
	};
};

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
		init();
	};
	
	
	~witness_set(){ // the destructor
		
		clear();
		
	};
	
	
	// assignment
	witness_set& operator=( const witness_set& other) {
		
		init();
		
		copy(other);
    return *this;
  };
	
	
	//copy operator.  must be explicitly declared because the underlying c structures use pointers.
	witness_set(const witness_set & other){
		
		init();
		
		copy(other);
	};
	
	void init()
	{
		this->input_filename = "unset_filename";
		
		this->num_variables = 0;
		this->num_synth_vars = 0;
		
		
		this->num_patches = 0;
		this->num_linears = 0;
		this->num_pts = 0;
		
		this->patch_mp = NULL;
		this->L_mp = NULL;
		this->pts_mp = NULL;
		
		this->incidence_number = -1;
		this->comp_num = this->dim = -1;
	}
	
	void copy(const witness_set & other)
	{
		this->input_filename = other.input_filename;

		this->dim = other.dim;
		this->comp_num = other.comp_num;
		this->incidence_number = other.incidence_number;
		
		this->num_variables = other.num_variables;
		this->num_synth_vars = other.num_synth_vars;
		
		this->variable_names = other.variable_names;
		
		copy_points(other);
        copy_linears(other);
		copy_patches(other);
	}
	
    
    void copy_linears(const witness_set & other) {
        for (int ii=0; ii<other.num_linears; ii++)
			add_linear(other.L_mp[ii]);
    }
    
    void copy_points(const witness_set & other) {
        for (int ii=0; ii<other.num_pts; ii++)
			add_point(other.pts_mp[ii]);
    }
    
    void copy_patches(const witness_set & other) {
        for (int ii=0; ii<other.num_patches; ii++)
			add_patch(other.patch_mp[ii]);
    }
	
    int num_natural_vars() const
    {
        return num_variables - num_synth_vars;
    }
    
    
	void only_natural_vars();
	void only_first_vars(int num_vars);
	void sort_for_real(tracker_config_t T);
	void sort_for_unique(tracker_config_t T);
	void sort_for_inside_sphere(comp_mp radius, vec_mp center);
	
	int  witnessSetParse(const boost::filesystem::path witness_set_file, const int num_vars);
	
	void reset()
	{
		clear();
	};
	
	
	
	void clear()
	{
		
		reset_names();
		
		reset_points();
		
		reset_linears();
		
		reset_patches();
		
		init();
	}
	
	
	void reset_names()
	{
		variable_names.resize(0);
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
		
		if (this->num_patches>0) {
			free(this->patch_mp);
		}
		
		this->num_patches = 0;
		this->patch_mp = NULL;
	};
	
	
	
	void add_patch(vec_mp new_patch);
	void add_point(vec_mp new_point);
	void add_linear(vec_mp new_linear);
	
	void cp_names(const witness_set & W_in);
	void cp_linears(const witness_set & W_in);
	void cp_patches(const witness_set & W_in);
	
	
	void merge(const witness_set & W_in);///< merges W_in into this
	
	
	void get_variable_names(); ///< reads variable names from names.out
	
	
	void print_to_screen() const; ///< prints some information about the witness set to the screen
	
	void print_to_file(boost::filesystem::path filename) const;
	
	/**
	 writes the linears in point form to file filename
	 
	 \param filename the name of the file to be written.
	 */
	void write_linears(boost::filesystem::path filename) const;
	
	/**
	 writes the patches in point form to file filename
	 
	 \param filename the name of the file to be written.
	 */
	void print_patches(boost::filesystem::path filename) const;
	void read_patches_from_file(boost::filesystem::path filename);
	
	
	void write_homogeneous_coordinates(boost::filesystem::path filename) const;
	void write_dehomogenized_coordinates(boost::filesystem::path filename) const;
	
    void send(parallelism_config & mpi_config, int target);
    void receive(parallelism_config & mpi_config);
};


















// CURVE CELL DECOMP DATA TYPES
class vertex
{
public:
	
	point_mp pt_mp;
	
	
	vec_mp  projection_values;
	
	int type;  //See enum above.
	int removed;
	int input_filename_index;
	
	vertex()
	{
		init();
	};
	
	~vertex()
	{
		clear();
		
	};
	
	vertex & operator=(const vertex & other)
	{
		init();
		copy(other);
		return *this;
	}
	
	vertex(const vertex& other)
	{
		init();
		copy(other);
	}
	
	void print() const
	{
		print_point_to_screen_matlab(pt_mp,"point");
		print_point_to_screen_matlab(projection_values,"projection_values");
		std::cout << "type: " << type << std::endl;
	}
	
	
	void send(int target, parallelism_config & mpi_config);
	
	void receive(int source, parallelism_config & mpi_config);
	
private:
	
	void clear()
	{
		clear_vec_mp(this->pt_mp);
		clear_vec_mp(this->projection_values);
	}
	
	void copy(const vertex & other)
	{
		vec_cp_mp(this->pt_mp, other.pt_mp);
		vec_cp_mp(this->projection_values, other.projection_values);
		this->type = other.type;
		
		this->input_filename_index = other.input_filename_index;
		this->removed = other.removed;
	}
	
	
	void init()
	{
		init_point_mp2(this->projection_values,0,1024);
		init_point_mp2(this->pt_mp,0,1024);
		this->pt_mp->size = this->projection_values->size = 0;
		this->type = UNSET;
		
		this->removed = 0;
		this->input_filename_index = -1;
	}
};



/**
 the main structure for storing vertices.  
 there are methods in place to add vertices, and perform lookups.
 */
class vertex_set
{
public:
	
	vec_mp *projections;
	int num_projections;
	int curr_projection;
	
	int curr_input_index;
	std::vector< boost::filesystem::path > filenames;
	
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
		
		set_num_vars(num_vars);
	}
	

	

	
	
	vertex_set & operator=( const vertex_set& other) {
		init();
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
	
	
	
	void set_num_vars(int num_vars)
	{
		this->num_natural_variables = num_vars;
		
		change_size_vec_mp(checker_1, num_vars);
		change_size_vec_mp(checker_2, num_vars);
		checker_1->size = checker_2->size = num_vars;
	}
	
	
	
	
	void print(boost::filesystem::path outputfile) const;
	
	
	
	
	int set_curr_input(boost::filesystem::path el_nom){
		
		int nom_index = -1;
        
		if ( el_nom.string().compare("unset_filename")==0 )
		{
			std::cout << "trying to set curr_input from unset_filename" << std::endl;
			
			br_exit(-59251);
		}
		
		
		int counter = 0;
		for (std::vector<boost::filesystem::path>::iterator ii = filenames.begin(); ii!= filenames.end(); ++ii) {
			
			if (*ii == el_nom) {
				nom_index = counter;
				break;
			}
			counter++;
		}
		
		if (nom_index==-1) {
			filenames.push_back(el_nom);
			nom_index = counter;
		}
		
		curr_input_index = nom_index;
		return nom_index;
	}
	
	
    
    /**
     \param W witness set containing points of which we wish to retrieve projections values.
     \return  crit_downstairs the projection values of the input witness set, sorted for uniqueness and increasingness.
     \return midpoints_downstairs the bisection of each interval in crit_downstairs.
     \return the indices of the points in W.
     \param pi the projection we are retrieving projection values with respect to.
     
     */
    
    int search_for_point(vec_mp testpoint);
    int search_for_active_point(vec_mp testpoint);
    int search_for_removed_point(vec_mp testpoint);
    
    int compute_downstairs_crit_midpts(const witness_set & W,
                                       vec_mp crit_downstairs,
                                       vec_mp midpoints_downstairs,
                                       std::vector< int > & index_tracker,
                                       vec_mp pi);
    
    
    void assert_projection_value(const std::set< int > & relevant_indices, comp_mp new_value);
    void assert_projection_value(const std::set< int > & relevant_indices, comp_mp new_value, int proj_index);
	
	int set_curr_projection(vec_mp new_proj){
        
        int proj_index = get_proj_index(new_proj);
        
        if (proj_index==-1) {
            int init_size = new_proj->size;
            new_proj->size = this->num_natural_variables;
            
			proj_index = add_projection(new_proj);
            
            new_proj->size = init_size;
		}
		
        curr_projection = proj_index;
        
//        std::cout << "curr_projection is now: " << curr_projection << std::endl;
		return proj_index;
	}
	
    
    int get_proj_index(vec_mp proj) const
    {
        int init_size = proj->size;
        proj->size = this->num_natural_variables;
        
        
        int proj_index = -1;
        
        for (int ii=0; ii<num_projections; ii++) {
			if (isSamePoint_inhomogeneous_input(projections[ii],proj)) {
				proj_index = ii;
				break;
			}
		}
        
        proj->size = init_size;
        
        return proj_index;
    }
    
    
	
	int add_projection(vec_mp proj){
		
		if (this->num_projections==0) {
			projections = (vec_mp *) br_malloc(sizeof(vec_mp));
		}
		else{
			this->projections = (vec_mp *)br_realloc(this->projections, (this->num_projections+1) * sizeof(vec_mp));
		}
		
		
		init_vec_mp2(this->projections[num_projections],num_natural_variables,DEFAULT_MAX_PREC);
		this->projections[num_projections]->size = num_natural_variables;
		
		
		if (proj->size != num_natural_variables) {
			vec_mp tempvec;
			init_vec_mp2(tempvec,num_natural_variables,DEFAULT_MAX_PREC); tempvec->size = num_natural_variables;
			for (int kk=0; kk<num_natural_variables; kk++) {
				set_mp(&tempvec->coord[kk], &proj->coord[kk]);
			}
			vec_cp_mp(projections[num_projections], tempvec);
			clear_vec_mp(tempvec);
		}
		else
		{
			vec_cp_mp(projections[num_projections], proj);
		}
		

		num_projections++;
		
		return num_projections;
	}
	
	
	
	
	void send(int target, parallelism_config & mpi_config);
	
	
	
	void receive(int source, parallelism_config & mpi_config);
	
	
private:
	
	void init()
	{
		
		num_projections = 0;
		projections = NULL;
		curr_projection = -1;
		
		curr_input_index = -2;
		
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
		this->curr_projection = other.curr_projection;
		if (this->num_projections==0 && other.num_projections>0) {
			projections = (vec_mp *) br_malloc(other.num_projections*sizeof(vec_mp));
		}
		else if(this->num_projections>0 && other.num_projections>0) {
			projections = (vec_mp *) br_realloc(projections,other.num_projections*sizeof(vec_mp));
		}
		else if (this->num_projections>0 && other.num_projections==0){
			for (int ii=0; ii<this->num_projections; ii++) {
				clear_vec_mp(projections[ii]);
			}
			free(projections);
		}
		
		for (int ii=0; ii<other.num_projections; ii++) {
			if (ii>=this->num_projections){
				init_vec_mp2(projections[ii],1,1024); projections[ii]->size = 1;
			}
			vec_cp_mp(projections[ii],other.projections[ii]);
		}
		
		this->num_projections = other.num_projections;
		
		
		curr_input_index = other.curr_input_index;
		filenames = other.filenames;
		
		
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
	int midpt; ///< index into vertex set
	
	void copy(const cell & other){
		this->midpt = other.midpt;
	}
	
	void send(int target, parallelism_config & mpi_config)
	{
		int buffer = midpt;
		MPI_Send(&buffer, 1, MPI_INT, target, CELL, MPI_COMM_WORLD);
	}
	
	void receive(int source, parallelism_config & mpi_config)
	{
		MPI_Status statty_mc_gatty;
		int buffer;
		MPI_Recv(&buffer, 1, MPI_INT, source, CELL, MPI_COMM_WORLD, &statty_mc_gatty);
		midpt = buffer;
	}
	
	
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
	
	std::vector< int > removed_points;
	
	
	
	edge() : cell()
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
	
	inline bool is_degenerate()
	{
		if ((left == right) || (left==midpt) || (right==midpt))
			return true;
		else
			return false;
	}
	
	void send(int target, parallelism_config & mpi_config)
	{
		int * buffer = new int[4];
		
		buffer[0] = left;
		buffer[1] = midpt;
		buffer[2] = right;
		buffer[3] = removed_points.size();
		
		MPI_Send(buffer, 4, MPI_INT, target, EDGE, MPI_COMM_WORLD);
		
		delete [] buffer;
		
		buffer = new int[removed_points.size()];
		for (int ii=0; ii<removed_points.size(); ii++) {
			buffer[ii] = removed_points[ii];
		}
		MPI_Send(buffer, removed_points.size(), MPI_INT, target, EDGE, MPI_COMM_WORLD);
		delete [] buffer;
		

	}
	
	void receive(int source, parallelism_config & mpi_config)
	{
		MPI_Status statty_mc_gatty;
		int * buffer = new int[4];
		MPI_Recv(buffer, 4, MPI_INT, source, EDGE, MPI_COMM_WORLD, &statty_mc_gatty);
		
		left  = buffer[0];
		midpt = buffer[1];
		right = buffer[2];
		
		int temp_num_removed = buffer[3];
		
		
		delete [] buffer;
		
		buffer = new int[temp_num_removed];
		MPI_Recv(buffer, temp_num_removed, MPI_INT, source, EDGE, MPI_COMM_WORLD, &statty_mc_gatty);
		for (int ii=0; ii<temp_num_removed; ii++) {
			removed_points.push_back(buffer[ii]);
		}
		
		delete [] buffer;
		
	//		int left;  ///< index into vertices
	//		int right; ///< index into vertices
	//		int midpt; ///<  index into vertices
	//		
	//		std::vector< int > removed_points;
	}
};












class decomposition
{

public:
	std::map< int , int > counters;
	std::map< int , std::vector< int > > indices;
	
	witness_set W;
	
	int num_variables;
	int dimension;
	int component_num;
	
	int num_curr_projections;
	vec_mp	*pi; // the projections
	
	std::vector< int > randomized_degrees;
	mat_mp randomizer_matrix;
	
	int num_patches;
	vec_mp *patch;

	vec_mp sphere_center;
	comp_mp sphere_radius;
	bool have_sphere_radius;
	
	boost::filesystem::path input_filename;
//	function input_file;
	
	
	
    

	void add_projection(vec_mp proj){
		if (this->num_curr_projections==0) {
			pi = (vec_mp *) br_malloc(sizeof(vec_mp));
		}
		else{
			this->pi = (vec_mp *)br_realloc(this->pi, (this->num_curr_projections+1) * sizeof(vec_mp));
		}
		
		init_vec_mp2(this->pi[num_curr_projections],proj->size,DEFAULT_MAX_PREC);
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
		
		init_vec_mp2(this->patch[num_patches],new_patch->size,DEFAULT_MAX_PREC);
		this->patch[this->num_patches]->size = new_patch->size;
		
		vec_cp_mp(this->patch[this->num_patches], new_patch);
		this->num_patches++;
	}
	
	
    int add_witness_set(const witness_set & W, int add_type, vertex_set & V);
    
	int index_in_vertices(vertex_set &V,
                          vec_mp testpoint);
	
	int index_in_vertices_with_add(vertex_set &V,
                                   vertex vert);
	
	int setup(boost::filesystem::path INfile);
	
	virtual void print(boost::filesystem::path outputfile);
	
	
	
	
	int read_sphere(const boost::filesystem::path & bounding_sphere_filename);
	
	void compute_sphere_bounds(const witness_set & W_crit);
	
	void copy_sphere_bounds(const decomposition & other)
	{
		set_mp(this->sphere_radius, other.sphere_radius);
		vec_cp_mp(this->sphere_center, other.sphere_center);
		this->have_sphere_radius = true;
	}
	
	
	
	void output_main(const boost::filesystem::path base);
	
	
	
	
	
	
	
	void reset()
	{
		clear();
		init();
	}
	
	
	
	
	decomposition(){
		init();
	}
	
	virtual ~decomposition()
	{
		this->clear();
	}
	
	decomposition & operator=(const decomposition& other){
		
		this->init();
		
		this->copy(other);
		
		return *this;
	}
	
	decomposition(const decomposition & other){
		
		this->init();
		
		this->copy(other);
	}
	
	
	
	void send(int target, parallelism_config & mpi_config);
	void receive(int source, parallelism_config & mpi_config);
	
protected:
	
	int add_vertex(vertex_set &V, vertex source_vertex);
	
	
	
	void init(){
		

		input_filename = "unset";
		pi = NULL;
		patch = NULL;
		init_mat_mp2(randomizer_matrix, 0, 0,1024);
		randomizer_matrix->rows = randomizer_matrix->cols = 0;
		
		randomized_degrees.clear();
		
		num_curr_projections = num_patches = 0;
		num_variables = 0;
		dimension = -1;
		component_num = -1;
		
		init_mp2(sphere_radius,DEFAULT_MAX_PREC);
		init_vec_mp2(sphere_center,0,1024);
		sphere_center->size = 0;
		have_sphere_radius = false;
		
		set_one_mp(sphere_radius);
		neg_mp(sphere_radius,sphere_radius);
		

		


		
		
		
	}
	
	
	void copy(const decomposition & other)
	{
		

		
		
		
		this->randomized_degrees = other.randomized_degrees;
		
		
		this->input_filename = other.input_filename;
		
		this->counters = other.counters;
		this->indices = other.indices;
		this->num_variables = other.num_variables;
		this->dimension = other.dimension;
		this->component_num = other.component_num;
		
		this->num_curr_projections = other.num_curr_projections;
		this->pi = (vec_mp *) br_malloc(other.num_curr_projections * sizeof(vec_mp));
		for (int ii = 0; ii<other.num_curr_projections; ii++) {
			init_vec_mp2(this->pi[ii],other.pi[ii]->size,DEFAULT_MAX_PREC);
			this->pi[ii]->size = other.pi[ii]->size;
			vec_cp_mp(this->pi[ii], other.pi[ii])
		}
		
		
		init_mat_mp2(this->randomizer_matrix, 0, 0,1024); randomizer_matrix->rows = randomizer_matrix->cols = 0;
		mat_cp_mp(this->randomizer_matrix, other.randomizer_matrix);
		
		
		
		this->num_patches = other.num_patches;
		this->patch = (vec_mp *) br_malloc(other.num_patches * sizeof(vec_mp));
		for (int ii = 0; ii<other.num_patches; ii++) {
			init_vec_mp2(this->patch[ii],other.patch[ii]->size,DEFAULT_MAX_PREC);
			this->patch[ii]->size = other.patch[ii]->size;
			vec_cp_mp(this->patch[ii], other.patch[ii])
		}
		
		copy_sphere_bounds(other);
		
		return;
	}
	
	void clear()
	{
		

		randomized_degrees.clear();
		
		if (num_curr_projections>0){
			for (int ii=0; ii<num_curr_projections; ii++)
				clear_vec_mp(pi[ii]);
			free(pi);
		}
		num_curr_projections = 0;
		
		if (num_patches>0){
			for (int ii=0; ii<num_patches; ii++)
				clear_vec_mp(patch[ii]);
			free(patch);
		}
		num_patches = 0;
		
		counters.clear();
		indices.clear();
		randomized_degrees.clear();
		
		clear_mat_mp(randomizer_matrix);
		
		clear_mp(sphere_radius);
		clear_vec_mp(sphere_center);
		

		
	}
	
	
	
	

	
}; // end decomposition














namespace color {
	
	int color_to_int(const std::string c);
	
	
	std::string bold(std::string new_color);
	
	std::string dark(std::string new_color);
	
	
	std::string underline(std::string new_color);
	
	
	std::string background(std::string new_color);
	
	
	std::string strike(std::string new_color);
			 
	std::string console_default();
	
	std::string black();
	
	std::string red();
	
	std::string green();

	
	std::string brown();
	
	std::string blue();
	
	std::string magenta();
	
	std::string cyan();
	
	std::string gray();
	
	
	//black - 30
	//red - 31
	//green - 32
	//brown - 33
	//blue - 34
	//magenta - 35
	//cyan - 36
	//lightgray - 37
	
}


#endif

