#include "checkSelfConjugate.h"







int checkSelfConjugate(witness_set_d W,
                       int           num_vars,
                       char          *input_file)
/***************************************************************\
 * USAGE: check if component is self conjugate                  *
 * ARGUMENTS: witness set, # of variables and name of input file*
 * RETURN VALUES: 1-self conjugate; 0-non self conjugate   *
 * NOTES:                                                        *
 \***************************************************************/
{
  int strLength = 0, digits = 15, *declarations = NULL;
  char *SysStr = NULL,*fmt = NULL, *bertini_command="bertini";
  FILE *IN = NULL, *OUT=NULL;
	
	
	//make the command string to run
	strLength = 1 + snprintf(NULL, 0, "%s input_membership_test", bertini_command);
  SysStr = (char *)bmalloc(strLength * sizeof(char));
  sprintf(SysStr, "%s input_membership_test", bertini_command);
	
	
	
	// make the format to write the member_points file
  strLength = 1 + snprintf(NULL, 0, "%%.%dle %%.%dle\n", digits, digits);
  // allocate size
  fmt = (char *)bmalloc(strLength * sizeof(char));
  // setup fmt
  sprintf(fmt, "%%.%dle %%.%dle\n", digits, digits);
	
	
  // setup input file
	IN = safe_fopen_read(input_file);
  partitionParse(&declarations, IN, "func_input_real", "config_real",0); // the 0 means not self conjugate
	fclose(IN);

	
	
	//check existence of the required witness_data file.
	IN = safe_fopen_read("witness_data");
	fclose(IN);
	
	//perhaps this could be used to get the witness_data file, but there would be problems
//  membership_test_input_file("input_membership_test", "func_input_real", "config_real",1);

	
	
	

	//only need to do this once.  we put the point and its conjugate into the same member points file and run the membership test simultaneously with one bertini call.
  membership_test_input_file("input_membership_test", "func_input_real", "config_real",3);
	
	
	//setup  member_points file, including both the first witness point, and its complex-conjugate
	write_member_points_sc(W.W.pts[0],fmt);
	
	// Do membership test
	printf("*\n%s\n*\n",SysStr);
  system(SysStr);
	
	int *component_numbers;
	component_numbers = (int *)bmalloc(2*sizeof(int));
	
	read_incidence_matrix(component_numbers);

  free(declarations);
	free(component_numbers);
	
	
  // delete temporary files
  remove("func_input_real");
  remove("config_real");
	remove("incidence_matrix");
	remove("member_points");
	if (component_numbers[0]==component_numbers[1]) {
		return 1;
	}
	else
	{
		return 0;
	}

	
}





int get_component_number(witness_set_d W,
												 int           num_vars,
												 char          *input_file)
/***************************************************************\
 * USAGE: check if component is self conjugate                  *
 * ARGUMENTS: witness set, # of variables and name of input file*
 * RETURN VALUES: 1-self conjugate; 0-non self conjugate   *
 * NOTES:                                                        *
 \***************************************************************/
{
  int strLength = 0, digits = 15, *declarations = NULL;
  char *SysStr = NULL,*fmt = NULL, *bertini_command="bertini";
  FILE *IN = NULL, *OUT=NULL;
	
	
	//make the command string to run
	strLength = 1 + snprintf(NULL, 0, "%s input_membership_test", bertini_command);
  SysStr = (char *)bmalloc(strLength * sizeof(char));
  sprintf(SysStr, "%s input_membership_test", bertini_command);
	
	
	
	// make the format to write the member_points file
  strLength = 1 + snprintf(NULL, 0, "%%.%dle %%.%dle\n", digits, digits);
  // allocate size
  fmt = (char *)bmalloc(strLength * sizeof(char));
  // setup fmt
  sprintf(fmt, "%%.%dle %%.%dle\n", digits, digits);
	
	
  // setup input file
	IN = safe_fopen_read(input_file);
  partitionParse(&declarations, IN, "func_input_real", "config_real",0); // the 0 means not self conjugate
	fclose(IN);
	
	
	
	//check existence of the required witness_data file.
	IN = safe_fopen_read("witness_data");
	fclose(IN);
	
	//perhaps this could be used to get the witness_data file, but there would be problems
	//  membership_test_input_file("input_membership_test", "func_input_real", "config_real",1);
	
	
	
	
	
	//only need to do this once.  we put the point and its conjugate into the same member points file and run the membership test simultaneously with one bertini call.
  membership_test_input_file("input_membership_test", "func_input_real", "config_real",3);
	
	
	//setup  member_points file, including both the first witness point, and its complex-conjugate
	write_member_points_singlept(W.W.pts[0],fmt);
	
	// Do membership test
	printf("*\n%s\n*\n",SysStr);
  system(SysStr);
	
	int *component_numbers;
	component_numbers = (int *)bmalloc(2*sizeof(int));
	
	read_incidence_matrix(component_numbers);
	
  free(declarations);
	free(component_numbers);
	
	
  // delete temporary files
  remove("func_input_real");
  remove("config_real");
	remove("incidence_matrix");
	remove("member_points");
	return component_numbers[0];
	
	
	
}



int write_member_points_singlept(point_d point_to_write, char * fmt){
	FILE *OUT = NULL;
	int ii;
	
	
	remove("member_points");
	OUT = safe_fopen_write("member_points");

	
	vec_d result;
	init_vec_d(result,0);
	dehomogenize(&result,point_to_write);
	
	fprintf(OUT,"1\n\n");
	for(ii=0;ii<result->size;ii++)
	{
		fprintf(OUT, fmt, result->coord[ii].r, result->coord[ii].i);
	}
	
	fclose(OUT);
	
	
	clear_vec_d(result);
	
	return 0;
}



int write_member_points_sc(point_d point_to_write, char * fmt){
	FILE *OUT = NULL;
	int ii;
	
	
	remove("member_points");
	OUT = safe_fopen_write("member_points");

	
	vec_d result;
	init_vec_d(result,0);
	dehomogenize(&result,point_to_write);
	
	fprintf(OUT,"2\n\n");
	for(ii=0;ii<result->size;ii++)
	{
		fprintf(OUT, fmt, result->coord[ii].r, result->coord[ii].i);
	}
	
	fprintf(OUT,"\n");
	for(ii=0;ii<result->size;ii++)
	{
		fprintf(OUT, fmt, result->coord[ii].r, -result->coord[ii].i);
	}
	fclose(OUT);
	
	clear_vec_d(result);
	
	return 0;
}




void read_incidence_matrix(int *component_numbers){
	FILE *IN = NULL;
	int ii,jj;
	int num_nonempty_codims, num_components, num_pts, codim;
	int component_indicator, component_number,total_num_components=0;

  //read incidence_matrix and see if it is self-conjugate
  IN = fopen("incidence_matrix", "r");
  if (IN == NULL)
  {
    printf("\n\nERROR: indicence_matrix file not produced.\n");
    bexit(ERROR_FILE_NOT_EXIST);
  }
	
	
  fscanf(IN, "%d",&num_nonempty_codims);      // number of nonempty codimensions
	
	for (ii = 0; ii<num_nonempty_codims; ii++) {
		fscanf(IN, "%d",&codim);      // codimension  (iterated for each codimension?)
		fscanf(IN, "%d",&num_components);  // number of components (is whatever)
		total_num_components = total_num_components+num_components;
	}

	
  fscanf(IN, "%d",&num_pts);    // number of points (should be 1)
	
	
	//and then a binary matrix indicating membership on which component
	//from the appendices of the book:
	//	% Binary matrix with one row per point, columns corresponding to components.
	//	%                0 if given point is not on given component;
	//	%                1 else .
	for (jj=0; jj<num_pts; jj++) {
		component_number=-10;
		for(ii=0;ii<total_num_components;ii++)  // i don't this is correct if there is more than one nonempty codimension.
			//todo: check this iterator limit  !!!
		{
			fscanf(IN, "%d", &component_indicator);
			if (component_indicator==1) {  //then is on this component
				component_number = ii;
			}
		}
		if (component_number==-10) {
			printf("it appears the membership test FAILED.\n");
//			exit(-1);
		}
		component_numbers[jj]=component_number;
	}
	
	fclose(IN);


	
	return;
}



void read_incidence_matrix_wrt_number(int *component_numbers, int given_incidence_number){
	FILE *IN = NULL;
	int ii,jj;
	int num_nonempty_codims, num_components, num_pts, codim;
	int component_indicator, component_number,total_num_components=0;
	
  //read incidence_matrix and see if it is self-conjugate
  IN = fopen("incidence_matrix", "r");
  if (IN == NULL)
  {
    printf("\n\nERROR: indicence_matrix file not produced.\n");
    bexit(ERROR_FILE_NOT_EXIST);
  }
	
	
  fscanf(IN, "%d",&num_nonempty_codims);      // number of nonempty codimensions
	
	for (ii = 0; ii<num_nonempty_codims; ii++) {
		fscanf(IN, "%d",&codim);      // codimension  (iterated for each codimension?)
		fscanf(IN, "%d",&num_components);  // number of components (is whatever)
		total_num_components = total_num_components+num_components;
	}
	
	
  fscanf(IN, "%d",&num_pts);    // number of points (should be 1)
	
	
	//and then a binary matrix indicating membership on which component
	//from the appendices of the book:
	//	% Binary matrix with one row per point, columns corresponding to components.
	//	%                0 if given point is not on given component;
	//	%                1 else .
	for (jj=0; jj<num_pts; jj++) {
		component_number=-10;
		for(ii=0;ii<total_num_components;ii++)  // i don't know if this is correct if there is more than one nonempty codimension.
		{
			fscanf(IN, "%d", &component_indicator);
			if (component_indicator==1 && given_incidence_number==ii) {  //then is on this component
				component_number = 1;
			}
		}
		if (component_number==-10) {
			printf("it appears the membership test FAILED.\n");
			//			exit(-1);
		}
		component_numbers[jj]=component_number;
	}
	
	fclose(IN);

	
	return;
}





void membership_test_input_file(char *outputFile,
                                char *funcInput,
                                char *configInput,
                                int  tracktype)
/***************************************************************\
 * USAGE: setup input file to membership test                    *
 * ARGUMENTS: name of output file, function & configuration input*
 * RETURN VALUES: none                                           *
 * NOTES:                                                        *
 \***************************************************************/
{
  char ch;
  FILE *OUT = fopen(outputFile, "w"), *IN = NULL;
  if (OUT == NULL)
  {
    printf("\n\nERROR: could not open '%s' to write!\n\n\n", outputFile);
    bexit(ERROR_FILE_NOT_EXIST);
  }
	
  // setup configurations in OUT
  fprintf(OUT, "CONFIG\n");
  IN = fopen(configInput, "r");
  if (IN == NULL)
  {
    fclose(OUT);
    printf("\n\nERROR: '%s' does not exist!!!\n\n\n", funcInput);
    bexit(ERROR_FILE_NOT_EXIST);
  }
	
	while ((ch = fgetc(IN)) != EOF)
    fprintf(OUT, "%c", ch);
  fclose(IN);
  if(tracktype==3)
    fprintf(OUT, "TrackType: %d;\nDeleteTempFiles: 1;\nEND;\nINPUT\n",tracktype);
  else
    fprintf(OUT, "TrackType: %d;\nDeleteTempFiles: 1;\nEND;\nINPUT\n",tracktype);
	
  // setup system in OUT
  IN = fopen(funcInput, "r");
  if (IN == NULL)
  {
    fclose(OUT);
    printf("\n\nERROR: '%s' does not exist!!!\n\n\n", funcInput);
    bexit(ERROR_FILE_NOT_EXIST);
  }
  while ((ch = fgetc(IN)) != EOF)
    fprintf(OUT, "%c", ch);
  fclose(IN);
  fprintf(OUT, "END;\n");
  fclose(OUT);
	
  return;
}






