#ifndef BERTINI_EXTENSIONS_H
#define BERTINI_EXTENSIONS_H

#include "fileops.hpp"

#include "bertini_headers.hpp"


#define TOL_DOUBLE_PRECISION 1e-13
#define LARGECHANGE_DOUBLEPRECISION 1e14


#define TOL_MP 1e-40
#define LARGECHANGE_MP 1e50

enum {VEC_MP = 4000, VEC_D, MAT_MP, MAT_D, COMP_MP, COMP_D, VEC_RAT, MAT_RAT, COMP_RAT, INDICES, DECOMPOSITION, CURVE, SURFACE, EDGE, CELL, FACE, UNUSED, VERTEX_SET, WITNESS_SET, VERTEX};

void *br_malloc(size_t size);

void *br_realloc(void *ptr, size_t size);



bool is_identity(mat_d M);
bool is_identity(mat_mp M);

//function prototypes for bertini_real data clearing etc.

void norm_of_difference(mpf_t result, vec_mp left, vec_mp right);

void dehomogenize(point_d *result, point_d dehom_me);
void dehomogenize(point_d *result, point_d dehom_me, int num_variables);

void dehomogenize(point_mp *result, point_mp dehom_me);
void dehomogenize(point_mp *result, point_mp dehom_me, int num_variables);

void nonconj_transpose(mat_d Res, mat_d M);
void nonconj_transpose(mat_mp Res, mat_mp M);

void dot_product_d(comp_d result, vec_d one, vec_d two);
void dot_product_mp(comp_mp result, vec_mp one, vec_mp two);

void dot_product_mindim(comp_d result, vec_d left, vec_d right);
void dot_product_mindim(comp_mp result, vec_mp left, vec_mp right);


int take_determinant_d(comp_d determinant, mat_d source_matrix);
int take_determinant_mp(comp_mp determinant, mat_mp source_matrix);


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





#endif