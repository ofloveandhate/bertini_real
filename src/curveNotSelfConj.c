#include "curveNotSelfConj.h"




void computeCurveNotSelfConj(witness_set W_in,
                             vec_mp         pi_mp,
                             curveDecomp_d	*C,
														 vertex_set			*V,
                             int           num_vars,
                             char          *input_file,
														 program_configuration * program_options,
														 solver_configuration * solve_options)
/***************************************************************\
 * USAGE: compute the isolated points for Non-Self-Conjugate case            *
 * ARGUMENTS: witness set, projection, # of variables, name of input file*
 * RETURN VALUES: Curve decomposition*
 * NOTES:                                                        *
 \***************************************************************/
{
	
	// num_vars includes the homogeneous variable
	
	
  int ii,jj,strLength,num_sols,*declarations = NULL;
  char *strSys,*bertini_command="bertini";
  FILE *IN = NULL;
  vec_mp cur_sol,cur_sol_bar;
  init_vec_mp(cur_sol,num_vars); cur_sol->size = num_vars;
	set_one_mp(&cur_sol->coord[0]);
	
  init_vec_mp(cur_sol_bar,num_vars); cur_sol_bar->size = num_vars;
	set_one_mp(&cur_sol_bar->coord[0]);
	
  IN = fopen(input_file, "r");
  partitionParse(&declarations, IN, "func_input_nsc", "config_nsc",1);
	

  //generate input file
  diag_homotopy_input_file("input_NSC", "func_input_nsc","func_inputbar","config_nsc",W_in.L_d[0],num_vars-1);
  //generate start file
	diag_homotopy_start_file("start",  W_in);

  strLength = 1 + snprintf(NULL, 0, "%s input_NSC", bertini_command);
  strSys = (char *)bmalloc(strLength * sizeof(char));
  sprintf(strSys, "%s input_NSC", bertini_command);
  //run bertini
	
	copyfile("witness_data","witness_data_0");
	
	printf("%s\n",strSys);
  system(strSys);
	
  //read the real solutions
  IN = safe_fopen_read("real_solutions");
  
  fscanf(IN, "%d\n\n", &num_sols);

	
	
	vertex temp_vertex; init_vertex(&temp_vertex);
	change_size_vec_mp(temp_vertex.pt_mp,num_vars); temp_vertex.pt_mp->size = num_vars;

	temp_vertex.type = ISOLATED;
	
	
	comp_mp projection_value;  init_mp(projection_value);
	
	
	for(ii=0;ii<num_sols;ii++) {
		for(jj=0;jj<num_vars-1;jj++){
			mpf_inp_str(cur_sol->coord[jj+1].r, IN, 10);
			mpf_inp_str(cur_sol->coord[jj+1].i, IN, 10);

			mpf_inp_str(cur_sol_bar->coord[jj+1].r, IN, 10);
			mpf_inp_str(cur_sol_bar->coord[jj+1].i, IN, 10);
		}
    
    //check if x=x_bar
  
    if (isSamePoint_homogeneous_input_mp(cur_sol,cur_sol_bar)) { // x=x_bar
			
			vec_cp_mp(temp_vertex.pt_mp,cur_sol);
			
			dot_product_mp(projection_value, temp_vertex.pt_mp, pi_mp);
			
			if (curve_index_in_vertices(C,V,
																	temp_vertex.pt_mp,
																	projection_value,
																	solve_options->T)
					==-1)
				curve_add_vertex(C,V,temp_vertex);      
		}
		
	}
	
	clear_mp(projection_value);
  //set the number of vertices
  fclose(IN);
  //clear
  clear_vec_d(cur_sol);
  clear_vec_d(cur_sol_bar);
  free(declarations);
  
	// delete temporary files
  remove("func_input_nsc");
  remove("config_nsc");
  remove("func_inputbar");
  remove("var_names");
	
	printf("renaming\n");
	rename("witness_data_0","witness_data");
}




// MISC. FUNCTIONS


void diag_homotopy_input_file(char  *outputFile, 
                              char  *funcInputx,
                              char  *funcInputy, 
                              char  *configInput,
                              vec_d L,
                              int   num_vars)
/***************************************************************\
* USAGE: setup input file to do diagonal homotopy             *
* ARGUMENTS: name of output file, function & configuration input*
* RETURN VALUES: none                                           *
* NOTES:                                                        *
\***************************************************************/
{
  char ch,**str,*fmt = NULL;
  int ii,jj,size;
  mat_d A;
  FILE *OUT = fopen(outputFile, "w"), *IN = NULL;
  if (OUT == NULL)
  {
    printf("\n\nERROR: '%s' is an inproper name of a file!!\n\n\n", outputFile);
    bexit(ERROR_FILE_NOT_EXIST);
  }
  str=(char **)bmalloc(num_vars*sizeof(char *));
  for(ii=0;ii<num_vars;ii++)
    str[ii]=(char*)bmalloc(sizeof(char)*256);
  
  // setup configurations in OUT
  fprintf(OUT, "CONFIG\n");
  IN = fopen(configInput, "r");
  if (IN == NULL)
  {
    fclose(OUT);
    printf("\n\nERROR: '%s' does not exist!!!\n\n\n", configInput);
    bexit(ERROR_FILE_NOT_EXIST);
  }
  fclose(IN);
  fprintf(OUT, "USERHOMOTOPY: 1;\nDeleteTempFiles: 1;\nEND;\nINPUT\n");

  // setup variables in OUT
  IN = fopen(funcInputx, "r");
  if (IN == NULL)
  {
    fclose(OUT);
    printf("\n\nERROR: '%s' does not exist!!!\n\n\n", funcInputx);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  while ((ch = fgetc(IN)) != EOF )
    fprintf(OUT, "%c", ch);
  fclose(IN);

  IN = fopen(funcInputy, "r");
  if (IN == NULL)
  {
    fclose(OUT);
    printf("\n\nERROR: '%s' does not exist!!!\n\n\n", funcInputy);
    bexit(ERROR_FILE_NOT_EXIST);
  }
  //setup the function name in OUT
  while ((ch = fgetc(IN)) != EOF )
    fprintf(OUT, "%c", ch);
  fclose(IN);
  IN = fopen("var_names", "r");
  if (IN == NULL)
  {
    fclose(OUT);
    printf("\n\nERROR: var_names does not exist!!!\n\n\n");
    bexit(ERROR_FILE_NOT_EXIST);
  }
  ii=0;jj=0;
  while ((ch = fgetc(IN)) != EOF)
  {
    if(ch!='\n')
      str[ii][jj++]=ch;
    else
    {
      str[ii++][jj]='\0';
      jj=0;
    }
  }
  //setup the linear equations 
  // find the size needed
  size = 1 + snprintf(NULL, 0, "%%.%dlf+%%.%dlf*I", 15, 15);
  // allocate size
  fmt = (char *)bmalloc(size * sizeof(char));
  // setup fmt
  sprintf(fmt, "%%.%dlf+%%.%dlf*I", 15, 15);
  // output the linear function L and L_bar
  for (ii = 0; ii < L->size; ii++)
  {
    fprintf(OUT, "L%d = ",ii);
    // print output
    fprintf(OUT, fmt, L->coord[ii].r, L->coord[ii].i); 
    fprintf(OUT, ";\n");

    fprintf(OUT, "Lbar%d = ",ii);
    // print output
    fprintf(OUT, fmt, L->coord[ii].r, -L->coord[ii].i); 
    fprintf(OUT, ";\n");
  } 
  fprintf(OUT, "\n");
  //Generate a random matrix A and output to input file.
  init_mat_d(A, 2, num_vars);
	make_matrix_random_d(A, 2, num_vars);
	print_matrix_to_screen_matlab(A,"A");
//  get_random_mat_d(A,2, num_vars);
  for (ii = 0; ii < 2; ii++)
    for(jj=0;jj<num_vars;jj++)
    {
      fprintf(OUT, "A%d%d = ",ii,jj);
      // print output
      fprintf(OUT, fmt, A->entry[ii][jj].r, A->entry[ii][jj].i); 
      fprintf(OUT, ";\n");
    } 
  //setup the diagonal homotopy functions
  fprintf(OUT, "\nL=t*(");
  //(Lx-1)*t+(1-t)*A[0]*(x-x_bar)
  for(ii=0;ii<num_vars;ii++)
  {
    fprintf(OUT, "L%d*", ii);
    fprintf(OUT, "%s", str[ii]);
    fprintf(OUT, "+");
  }
  fprintf(OUT, "-1)+(1-t)*(");
  for(ii=0;ii<num_vars;ii++)
  {
    fprintf(OUT, "A0%d*", ii);
    fprintf(OUT, "(%s-%s", str[ii],str[ii]);
    fprintf(OUT, "bar)+");
  }
  //(L_bar x_bar-1)*t+(1-t)*A[1]*(x-x_bar)
  fprintf(OUT, "0);\nLbar=t*(");
  for(ii=0;ii<num_vars;ii++)
  {
    fprintf(OUT, "Lbar%d*", ii);
    fprintf(OUT, "%sbar", str[ii]);
    fprintf(OUT, "+");
  }
  fprintf(OUT, "-1)+(1-t)*(");
  for(ii=0;ii<num_vars;ii++)
  {
    fprintf(OUT, "A1%d*", ii);
    fprintf(OUT, "(%s-%s", str[ii],str[ii]);
    fprintf(OUT, "bar)+");
  }
  fprintf(OUT, "0);\nEND;");
  fclose(OUT);
  //free
  for(ii=0;ii<num_vars;ii++)
    free(str[ii]);
  free(str);
  free(fmt);
  clear_mat_d(A);
  return;
}



void diag_homotopy_start_file(char                 *startFile,
                              witness_set  W)
/***************************************************************\
 * USAGE: setup start file to do diagonal homotopy             *
 * ARGUMENTS: name of output file, start points & number of variables*
 * RETURN VALUES: none                                           *
 * NOTES:                                                        *
 \***************************************************************/
{
  FILE *OUT = fopen(startFile, "w");
  char *fmt = NULL;
  int ii,jj,kk;
	int size,digits=15;
  if (OUT == NULL)
  {
    printf("\n\nERROR: '%s' is an improper name of a file!!\n\n\n", startFile);
    bexit(ERROR_FILE_NOT_EXIST);
  }
  size = 1 + snprintf(NULL, 0, "%%.%de %%.%de\n", digits, digits);
  // allocate size
  fmt = (char *)bmalloc(size * sizeof(char));
  // setup fmt & fmtb
  sprintf(fmt, "%%.%de %%.%de\n", digits, digits);
	// output the number of start points
  fprintf(OUT,"%d\n\n",W.num_pts*W.num_pts);
	
	
	comp_mp temp; init_mp(temp);
	
	
	for (kk=0; kk<W.num_pts; kk++){
  for (ii=0; ii<W.num_pts; ii++) { // output {w \bar{w}}'
		
    vec_mp result; init_vec_mp(result,0);
		vec_mp result2; init_vec_mp(result2,0);
		
    dehomogenize_mp(&result,W.pts_mp[ii]);
		dehomogenize_mp(&result2,W.pts_mp[kk]);
		
		
    for(jj=0; jj<W.num_variables-1;jj++) {
			print_mp(OUT, 15, &result->coord[jj]); fprintf(OUT, "\n");
			
			conjugate_mp(temp, &result2->coord[jj] )
			print_mp(OUT, 15, temp); fprintf(OUT, "\n");
//      fprintf(OUT, fmt, result->coord[jj].r,  result->coord[jj].i);
    }
		
	
		
    fprintf(OUT,"\n");
  }// re: ii
	}// re: kk
	
  free(fmt);
  fclose(OUT);
}






