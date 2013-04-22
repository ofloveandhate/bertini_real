

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



#ifndef _SAMPLER_H
#define _SAMPLER_H

#include "fileops.h"
#include "polysolve.h"
#include "data_type.h"
#include "programStartup.h"
#include "witnessSetParse.h"
#include "lintolinSolver.h"

int  setup_curveDecomp(int argC, char *args[], char **inputName, char **witnessSetName, char **RandMatName,  char **samplingNamenew, curveDecomp_d *C,int *num_vars);
int setup_edges(edge_d **E,char *INfile,int *num_vars, char **inputName, char *directoryName, int MPType);
int setup_vertices(vertex_d **V,char *INfile,int num_vars, int MPType);
int  Load_sampling_data(sample_d *S, curveDecomp_d C,int num_vars,int MPType);
void  output_sampling_data(sample_d S,char *samplingName,int num_vars, int MPType);
void  set_witness_set_d(witness_set_d *W, vec_d L,vec_d pts,int num_vars);
void set_witness_set_mp(witness_set_d *W, vec_mp L,vec_mp pts,int num_vars);
void generate_new_sampling_pts(sample_d *S_new,mat_mp n_minusone_randomizer_matrix,sample_d S_old, curveDecomp_d C,witness_set_d W, int num_vars,unsigned int currentSeed,int  MPType);
void generate_new_sampling_pts_d(sample_d *S_new,mat_mp n_minusone_randomizer_matrix,sample_d S_old, curveDecomp_d C,witness_set_d W, int num_vars,unsigned int currentSeed,int  MPType);
void generate_new_sampling_pts_mp(sample_d *S_new,mat_mp n_minusone_randomizer_matrix,sample_d S_old, curveDecomp_d C,witness_set_d W, int num_vars,unsigned int currentSeed,int  MPType);

void read_rand_matrix(char *INfile, mat_mp n_minusone_randomizer_matrix);

#endif



