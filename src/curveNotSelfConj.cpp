#include "curveNotSelfConj.hpp"




void computeCurveNotSelfConj(witness_set W_in,
                             vec_mp         pi_mp,
                             curveDecomp_d	*C,
														 vertex_set			*V,
                             int           num_vars,
                             boost::filesystem::path input_file,
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
	
	
  int ii,jj,num_sols,*declarations = NULL;
  
	std::string bertini_system_command = "bertini input_NSC";
	
  FILE *IN = NULL;
  vec_mp cur_sol,cur_sol_bar;
  init_vec_mp(cur_sol,num_vars); cur_sol->size = num_vars;
	set_one_mp(&cur_sol->coord[0]);
	
  init_vec_mp(cur_sol_bar,num_vars); cur_sol_bar->size = num_vars;
	set_one_mp(&cur_sol_bar->coord[0]);
	
  IN = safe_fopen_read(input_file.c_str());
  partitionParse(&declarations, IN, const_cast<char *>("func_input_nsc"), const_cast<char *>("config_nsc"),1);
	

  //generate input file
  diag_homotopy_input_file("input_NSC", "func_input_nsc","func_inputbar","config_nsc",W_in.L_d[0],num_vars-1);
  //generate start file
	diag_homotopy_start_file("start",  W_in);


  //run bertini
		
	copyfile("witness_data","witness_data_0");
	
  system(bertini_system_command.c_str());

	rename("witness_data_0","witness_data");
	
	
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
			
			dot_product_mp(temp_vertex.projVal_mp, temp_vertex.pt_mp, pi_mp);
			
			if (curve_index_in_vertices(C,V,
																	temp_vertex.pt_mp,
																	temp_vertex.projVal_mp,
																	solve_options->T) ==-1
					)
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
	
	
}




// MISC. FUNCTIONS


void diag_homotopy_input_file(boost::filesystem::path outputFile,
                              boost::filesystem::path funcInputx,
                              boost::filesystem::path funcInputy, 
                              boost::filesystem::path configInput,
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
	
	FILE *IN = NULL;
	

  str=(char **)bmalloc(num_vars*sizeof(char *));
  for(ii=0;ii<num_vars;ii++)
    str[ii]=(char*)bmalloc(sizeof(char)*256);
  
	FILE *OUT = safe_fopen_write(outputFile.c_str());
	
  // setup configurations in OUT
  fprintf(OUT, "CONFIG\n");
  IN = safe_fopen_read(configInput.c_str());
  fclose(IN);
	
  fprintf(OUT, "USERHOMOTOPY: 1;\nDeleteTempFiles: 1;\nEND;\nINPUT\n");

  // setup variables in OUT
  IN = safe_fopen_read(funcInputx.c_str());
  while ((ch = fgetc(IN)) != EOF )
    fprintf(OUT, "%c", ch);
  fclose(IN);

	//setup the function name in OUT
	IN = safe_fopen_read(funcInputy.c_str());
  while ((ch = fgetc(IN)) != EOF )
    fprintf(OUT, "%c", ch);
  fclose(IN);
	
	
  IN = safe_fopen_read("var_names");
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



void diag_homotopy_start_file(boost::filesystem::path startFile,
                              witness_set  W)
/***************************************************************\
 * USAGE: setup start file to do diagonal homotopy             *
 * ARGUMENTS: name of output file, start points & number of variables*
 * RETURN VALUES: none                                           *
 * NOTES:                                                        *
 \***************************************************************/
{
  
  char *fmt = NULL;
  int ii,jj,kk;
	int size,digits=15;
	
	
	FILE *OUT = safe_fopen_write(startFile.c_str());
  
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






