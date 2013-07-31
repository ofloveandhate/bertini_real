	
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>


#include <iostream>
#include <ios>
#include <string>
#include <fstream>

#include <vector>
#include <sstream>



#ifndef WITNESS_SET_H_
#define WITNESS_SET_H_


extern "C" {
#include "polysolve.h"
}


#include "data_type.hpp"

#include "fileops.hpp"
#include "missing_bertini_headers.hpp"




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
	
	int MPType; //  will indicate the type of the solve.  both fields will have data, but one will be used.
	
	std::vector< std::string > variable_names;
	// end data members
	
	
	
	// overloaded operators
	
	
	// default constructor
	
	witness_set(){
		
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
		
		this->MPType = other.MPType;
		
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
		
		this->MPType = other.MPType;
		
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
	

	
	void only_first_vars(int num_vars);
	void sort_for_real(tracker_config_t T);
	void sort_for_unique(tracker_config_t T);
	
	
	int  witnessSetParse(const boost::filesystem::path witness_set_file, const int num_vars);

	void reset()
	{
		for (int ii =0; ii<this->num_linears; ii++) {
			clear_vec_mp(this->L_mp[ii]);
		}
		if (this->num_linears>0) {
			free(this->L_mp);
		}
		
		
		for (int ii =0; ii<this->num_patches; ii++) {
			clear_vec_mp(this->patch_mp[ii]);
		}
		if (this->num_patches>0) {
			free(this->patch_mp);
		}
		
		
		for (int ii =0; ii<this->num_pts; ii++) {
			clear_vec_mp(this->pts_mp[ii]);
		}
		if (this->num_pts>0) {
			free(this->pts_mp);
		}
		
		this->num_variables = 0;
		this->num_synth_vars = 0;
		
		this->num_pts = this->num_patches = this->num_linears = 0;
		
		this->patch_mp = NULL;
		this->L_mp = NULL;
		this->pts_mp = NULL;
		
		this->incidence_number = -1;
		this->comp_num = this->dim = -1;
	};
	
	void reset_patches()
	{
				
		
		for (int ii =0; ii<this->num_patches; ii++) {
			clear_vec_mp(this->patch_mp[ii]);
		}
		if (this->num_patches>0) {
			free(this->patch_mp);
		}
		
		this->num_patches = 0;
		
		this->patch_mp = NULL;
	};
	
	
	
	void add_patch(vec_mp new_patch);
	void add_point(vec_mp new_point);
	void add_linear(vec_mp new_linear);
	
	
	void merge(const witness_set & W_in);///< merges W_in into this
	
	/**
	 parses variable names from names.out
	 */
	void get_variable_names();
	
	/**
	 prints some information about the witness set to the screen
	 */
	void print_to_screen();
	
	void print_to_file();
	
	
	/**
	 writes the linears in point form to file filename
	 
	 \param filename the name of the file to be written.
	 */
	void write_linears(boost::filesystem::path filename);
	void print_patches(boost::filesystem::path filename);
	void read_patches_from_file(boost::filesystem::path filename);
	
	
	void write_homogeneous_coordinates(boost::filesystem::path filename);
	void write_dehomogenized_coordinates(boost::filesystem::path filename);
	
	
};  
// end the double types





void cp_names(witness_set *W_out, witness_set & W_in);
void cp_linears(witness_set *W_out, witness_set & W_in);
void cp_patches(witness_set *W_out, witness_set & W_in);






//
//void sort_for_real(witness_set *W_out,
//										witness_set & W_in,
//										tracker_config_t T);
//
//void sort_for_unique(witness_set *W_out,
//										 witness_set & W_in,
//										 tracker_config_t T);


#endif
