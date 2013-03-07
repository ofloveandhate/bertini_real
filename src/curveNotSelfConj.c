#include "Curve_NotSelfConj.h"





void get_random_mat_d(mat_d A, 
                      int   n,
                      int   m)
/***************************************************************\
* USAGE: obtain a complex random matrix in double precision     *
* ARGUMENTS:                                                    *
* RETURN VALUES: a complex random matrix A. The size is n by m          *
* NOTES:                                                        *
\***************************************************************/
{ 
  if (n > 0 && m>0)
  { 
    int i,j;
    double sum=0.0;
    // set sum to 0
    for (i = 0; i < n; i++)
		for(j = 0; j< m;j++)
	    { // get a random number
    	  get_comp_rand_d(&A->entry[i][j]);
      	  // add on norm squared
      	  sum=A->entry[i][j].r ;
  	      sum += (A->entry[i][j].r) * (A->entry[i][j].r) 
		  	   + (A->entry[i][j].i) * (A->entry[i][j].i);
    } 
    // compute 1/norm of vector
    sum=1/sqrt(sum);
    // normalize
    for (i = 0; i < n; i++)
      for(j = 0; j < m; j++)
      {
        (A->entry[i][j]).r *= sum;
        (A->entry[i][j]).i *= sum;
      }
  }
  return;
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
  int i,j,size;
  mat_d A;
  FILE *OUT = fopen(outputFile, "w"), *IN = NULL;
  if (OUT == NULL)
  {
    printf("\n\nERROR: '%s' is an inproper name of a file!!\n\n\n", outputFile);
    bexit(ERROR_FILE_NOT_EXIST);
  }
  str=(char **)bmalloc(num_vars*sizeof(char *));
  for(i=0;i<num_vars;i++)
    str[i]=(char*)bmalloc(sizeof(char)*256);
  
  // setup configurations in OUT
  fprintf(OUT, "CONFIG\n");
  IN = fopen(configInput, "r");
  if (IN == NULL)
  {
    fclose(OUT);
    printf("\n\nERROR: '%s' does not exist!!!\n\n\n", configInput);
    bexit(ERROR_FILE_NOT_EXIST);
  }
  while ((ch = fgetc(IN)) != EOF)
    fprintf(OUT, "%c", ch);
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
  fprintf(OUT, "variable ");
  while ((ch = fgetc(IN)) != ' ');
  i=0;j=0;
  //store variables in str strucutre. str[i] for i-th variable.
  while ((ch = fgetc(IN)) != '\n')
  {
    if(ch!=',' && ch!=';')
    { 
      str[i][j++]=ch;
      fprintf(OUT, "%c", ch);
    }
    else
    {
      str[i++][j]='\0';
      j=0;
      fprintf(OUT, ","); 
    }
  }
  for(i=0;i<num_vars;i++)
  {
    fprintf(OUT, "%sb", str[i]);
    if(i!=num_vars-1)
      fprintf(OUT,",");
    else
      fprintf(OUT,";\n");
  }

  // setup system in OUT
  //setup the function name in OUT
  while ((ch = fgetc(IN)) != ';')
    fprintf(OUT, "%c", ch);
  fprintf(OUT, ",");
  rewind(IN);
  while ((ch = fgetc(IN)) != '\n');
  while ((ch = fgetc(IN)) != ' ');
  while ((ch = fgetc(IN)) != '\n')
  {
    if(ch==',' || ch==';')
      fprintf(OUT, "b,");
    else
      fprintf(OUT, "%c", ch);
  }
  fprintf(OUT, "L,Lb;\n");	  

  fprintf(OUT,"pathvariable s;\n");
  fprintf(OUT,"parameter t;\n");
  fprintf(OUT,"t = s;\n\n");
  //output f(x) in OUT
  
  while ( (ch = fgetc(IN)) != EOF)
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
  while ((ch = fgetc(IN)) != '\n');
  while ((ch = fgetc(IN)) != '\n');
  while ((ch = fgetc(IN)) != EOF )
    fprintf(OUT, "%c", ch);
  fclose(IN);
  //setup the linear equations 
  // find the size needed
  size = 1 + snprintf(NULL, 0, "%%.%de+%%.%de*I", 15, 15);
  // allocate size
  fmt = (char *)bmalloc(size * sizeof(char));
  // setup fmt
  sprintf(fmt, "%%.%de+%%.%de*I", 15, 15);
  // output the linear function L and L_bar
  for (i = 0; i < L->size; i++)
  {
    fprintf(OUT, "L%d = ",i);
    // print output
    fprintf(OUT, fmt, L->coord[i].r, L->coord[i].i); 
    fprintf(OUT, ";\n");

    fprintf(OUT, "Lb%d = ",i);
    // print output
    fprintf(OUT, fmt, L->coord[i].r, -L->coord[i].i); 
    fprintf(OUT, ";\n");
  } 
  fprintf(OUT, "\n");
  //Generate a random matrix A and output to input file.
  init_mat_d(A, 2, num_vars);
  get_random_mat_d(A,2, num_vars);
  for (i = 0; i < 2; i++)
    for(j=0;j<num_vars;j++)
    {
      fprintf(OUT, "A%d%d = ",i,j);
      // print output
      fprintf(OUT, fmt, A->entry[i][j].r, A->entry[i][j].i); 
      fprintf(OUT, ";\n");
    } 
  //setup the diagonal homotopy functions
  fprintf(OUT, "\nL=t*(");
  //(Lx-1)*t+(1-t)*A[0]*(x-x_bar)
  for(i=0;i<num_vars;i++)
  {
    fprintf(OUT, "L%d*", i);
    fprintf(OUT, "%s", str[i]);
    fprintf(OUT, "+");
  }
  fprintf(OUT, "-1)+(1-t)*(");
  for(i=0;i<num_vars;i++)
  {
    fprintf(OUT, "A0%d*", i);
    fprintf(OUT, "(%s-%s", str[i],str[i]);
    fprintf(OUT, "b)+");
  }
  //(L_bar x_bar-1)*t+(1-t)*A[1]*(x-x_bar)
  fprintf(OUT, "0);\nLb=t*(");
  for(i=0;i<num_vars;i++)
  {
    fprintf(OUT, "Lb%d*", i);
    fprintf(OUT, "%sb", str[i]);
    fprintf(OUT, "+");
  }
  fprintf(OUT, "-1)+(1-t)*(");
  for(i=0;i<num_vars;i++)
  {
    fprintf(OUT, "A1%d*", i);
    fprintf(OUT, "(%s-%s", str[i],str[i]);
    fprintf(OUT, "b)+");
  }
  fprintf(OUT, "0);\nEND;");
  fclose(OUT);
  //free
  for(i=0;i<num_vars;i++)
    free(str[i]);
  free(str);
  free(fmt);
  clear_mat_d(A);
  return;
}
void diag_homotopy_start_file(char                 *startFile, 
                              witness_point_set_d  W, 
                              int                  num_vars)
/***************************************************************\
* USAGE: setup start file to do diagonal homotopy             *
* ARGUMENTS: name of output file, start points & number of variables*
* RETURN VALUES: none                                           *
* NOTES:                                                        *
\***************************************************************/

{
  FILE *OUT = fopen(startFile, "w");
  char *fmt = NULL;
  int i,j,size,digits=15;
  if (OUT == NULL)
  {
    printf("\n\nERROR: '%s' is an inproper name of a file!!\n\n\n", startFile);
    bexit(ERROR_FILE_NOT_EXIST);
  }
  size = 1 + snprintf(NULL, 0, "%%.%de %%.%de\n", digits, digits);
  // allocate size
  fmt = (char *)bmalloc(size * sizeof(char));
  // setup fmt & fmtb
  sprintf(fmt, "%%.%de %%.%de\n", digits, digits);
   // output the number of start points	
  fprintf(OUT,"%d\n\n",W.num_pts);
  for (i = 0; i < W.num_pts; i++)
  { // output {w \bar{w}}'
    for(j=0; j<num_vars;j++)
      fprintf(OUT, fmt, W.pts[i]->coord[j].r, W.pts[i]->coord[j].i);
    for(j=0; j<num_vars;j++)
      fprintf(OUT, fmt, W.pts[i]->coord[j].r, -W.pts[i]->coord[j].i);
    fprintf(OUT,"\n");
  }
  free(fmt);
  fclose(OUT);

}
void computeCurveNotSelfConj(witness_set_d Wnew, 
                             vec_d         pi, 
                             curveDecomp_d *C,
                             int           num_vars,
                             char          *input_file)
/***************************************************************\
* USAGE: compute the isolated points for Non-Self-Conjugate case            *
* ARGUMENTS: witness set, projection, # of variables, name of input file*
* RETURN VALUES: Curve decomposition*
* NOTES:                                                        *
\***************************************************************/
{
  int i,j,k,m,strLength,num_sols,*declarations = NULL,isadd;
  char *strSys,*bertini_command="bertini";
  FILE *IN = NULL;
  vec_d cur_sol,cur_sol_bar;
  double d_norm=0.0,TOL=1e-15,tempD;
  comp_d temp;
  init_vec_d(cur_sol,num_vars);
  init_vec_d(cur_sol_bar,num_vars);
   
  IN = fopen(input_file, "r");
  partitionParse(&declarations, IN, "func_input_real", "config_real");

  printf("\nCompute Non Self Conjugate Curve\n");
  //generate input file
  diag_homotopy_input_file("input_NSC", "func_input_real","func_inputb","config_real",Wnew.L,num_vars);
  //generate start file
  diag_homotopy_start_file("start", Wnew.W, num_vars);
  strLength = 1 + snprintf(NULL, 0, "%s input_NSC", bertini_command);
  strSys = (char *)bmalloc(strLength * sizeof(char));
  sprintf(strSys, "%s input_NSC", bertini_command);
  //run bertini
  system(strSys);
  //read the real solutions
  IN = fopen("real_solutions", "r");
  if (IN == NULL)
  {
    printf("\n\nERROR: Bertini was unable to compute the diagonal homotopy.\n\n\n");
    bexit(ERROR_CONFIGURATION);
  }
  fscanf(IN, "%d", &num_sols);
  //insert the real solutions into V0
  C->V0=(vertex_d *)bmalloc(num_sols*sizeof(vertex_d));
  for(i=0,k=0;i<num_sols;i++)
  {
    for(j=0;j<num_vars;j++)
      fscanf(IN, "%lf %lf", &(cur_sol->coord[j].r), &(cur_sol->coord[j].i));
    for(j=0;j<num_vars;j++)
      fscanf(IN, "%lf %lf", &(cur_sol_bar->coord[j].r), &(cur_sol_bar->coord[j].i));
    //check if x=x_bar
    for(j=0;j<num_vars;j++)
    {
      sub_d(&(cur_sol_bar->coord[j]), &(cur_sol_bar->coord[j]), &(cur_sol->coord[j]));
      if(d_norm<d_abs_d(&(cur_sol_bar->coord[j])))
        d_norm=d_abs_d(&(cur_sol_bar->coord[j]));
    }
    if(d_norm<TOL)//x=x_bar
    {
      // check if V0 has current vertex already.
      isadd=1;
      for(j=0;j<k;j++) 
      {
        d_norm=0.0;
        for(m=0;m<num_vars;m++)
        {
          sub_d(&(cur_sol_bar->coord[m]), &(cur_sol->coord[m]), &(C->V0[j].pt->coord[m]));
          if(d_norm<d_abs_d(&(cur_sol_bar->coord[m])))
            d_norm=d_abs_d(&(cur_sol_bar->coord[m]));
        }
        if(d_norm<TOL)
        {
          isadd=0;break;//V0 has this point already.
        }
      }
      if(isadd)
      {
         // add points to V0
        init_point_d(C->V0[k].pt,num_vars);
        point_cp_d(C->V0[k].pt,cur_sol);
        //compute the projection
        set_zero_d(C->V0[k].projVal);
        for(j=0;j<num_vars;j++)
        {
          mul_d2(temp,&(C->V0[k].pt->coord[j]),&(pi->coord[j]),tempD);
          add_d(C->V0[k].projVal,temp,C->V0[k].projVal);
        }
        k++;
      }
    }
  }
  //set the number of vertices
  C->num_V0=k;
  fclose(IN);
  //clear
  clear_vec_d(cur_sol);
  clear_vec_d(cur_sol_bar);
  free(declarations);
  
   // delete temporary files
  remove("func_input_real");
  remove("config_real");

}



