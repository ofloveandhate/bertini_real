#include "checkSelfConjugate.h"







int checkSelfConjugate(witness_set_d W,
                       int           num_vars,
                       char          *input_file)
/***************************************************************\
 * USAGE: check if it is self conjugate                  *
 * ARGUMENTS: witness set, # of variables and name of input file*
 * RETURN VALUES: 1-self conjugate; 0-non self conjugate   *
 * NOTES:                                                        *
 \***************************************************************/
{
	int component_number_original, component_number_conjugated;
  int IsSelfConj=1, strLength = 0, digits = 15, *declarations = NULL;
  char *SysStr = NULL,*fmt = NULL, *bertini_command="bertini";
  FILE *IN = NULL;
	
	
	//make the command string to run
	strLength = 1 + snprintf(NULL, 0, "%s input_membership_test", bertini_command);
  SysStr = (char *)bmalloc(strLength * sizeof(char));
  sprintf(SysStr, "%s input_membership_test", bertini_command);
	
	
	
	// make the format to write the member_points file
  strLength = 1 + snprintf(NULL, 0, "%%.%de -%%.%de\n", digits, digits);
  // allocate size
  fmt = (char *)bmalloc(strLength * sizeof(char));
  // setup fmt
  sprintf(fmt, "%%.%dlf %%.%dlf\n", digits, digits);
	
	
  // setup input file
  IN = fopen(input_file, "r");
  partitionParse(&declarations, IN, "func_input", "config");
	fclose(IN);
	printf("done partitionParse\n");

	
	//check existence of the required witness_data file.
	IN = fopen("witness_data","r");
	if (IN == NULL) {
		printf("missing witness_data file\n");
		bexit(ERROR_FILE_NOT_EXIST);
	}
	else{
		printf("correctly detected witness_data file\n");
	}
	fclose(IN);
	
	//perhaps this could be used to get the witness_data file, but there would be problems
//  membership_test_input_file("input_membership_test", "func_input_real", "config_real",1);
//  system(SysStr);
	
	
	

	
	

	

	//only need to do this once, and then we feed it two different member_points files
  membership_test_input_file("input_membership_test", "func_input", "config",3);
	
	
	
	//setup  member_points file
	write_member_points(W.W.pts[0],num_vars,fmt);
	
	// Do membership test
	printf("*\n%s\n*\n",SysStr);
  system(SysStr);
	component_number_original = read_incidence_matrix();
	
	
	
	
	//now write the conjugated point, and check if its component number matches.
	
	//setup  member_points file
	write_member_points_conjugated(W.W.pts[0],num_vars,fmt);
	
	// Do membership test
	printf("*\n%s\n*\n",SysStr);
  system(SysStr);
	component_number_conjugated = read_incidence_matrix();
	
	
	
	
  free(declarations);
	
  // delete temporary files
  remove("func_input_real");
  remove("config_real");
	

	if (component_number_original == component_number_conjugated) {
		return 1;
	}
	else
	{
		return 0;
	}

  return IsSelfConj;
}


int write_member_points(point_d point_to_write, int num_vars, char * fmt){
	FILE *OUT = NULL;
	int ii;
	
	
	remove("member_points");
	OUT = fopen("member_points", "w");
	if (OUT == NULL)
	{
		printf("\n\nERROR: Can not open 'member_points' to write points!\n");
		bexit(ERROR_FILE_NOT_EXIST);
	}
	
	fprintf(OUT,"1\n\n");
	for(ii=0;ii<num_vars;ii++)
	{
		fprintf(OUT, fmt, point_to_write->coord[ii].r, point_to_write->coord[ii].i);
	}
	fclose(OUT);
	
	//old way of printing the points:
//    fprintf(OUT, fmt, WS.W.pts[0]->coord[i].r, -WS.W.pts[0]->coord[i].i);
	return 0;
}


int write_member_points_conjugated(point_d point_to_write, int num_vars, char * fmt){
	FILE *OUT = NULL;
	int ii;
	
	
	remove("member_points");
	OUT = fopen("member_points", "w");
	if (OUT == NULL)
	{
		printf("\n\nERROR: Can not open 'member_points' to write points!\n");
		bexit(ERROR_FILE_NOT_EXIST);
	}
	
	fprintf(OUT,"1\n\n");
	for(ii=0;ii<num_vars;ii++)
	{ //                                             difference right here, put the minus sign
		fprintf(OUT, fmt, point_to_write->coord[ii].r, -point_to_write->coord[ii].i);
	}
	fclose(OUT);
	
	//old way of printing the points:
	//    fprintf(OUT, fmt, WS.W.pts[0]->coord[i].r, -WS.W.pts[0]->coord[i].i);
	return 0;
}


int read_incidence_matrix(){
	FILE *IN = NULL;
	int ii;
	int num_nonempty_codims, num_components, num_pts, codim;
	int component_indicator, component_number;

  //read incidence_matrix and see if it is self-conjugate
  IN = fopen("incidence_matrix", "r");
  if (IN == NULL)
  {
    printf("\n\nERROR: indicence_matrix file not produced.\n");
    bexit(ERROR_FILE_NOT_EXIST);
  }
  fscanf(IN, "%d",&num_nonempty_codims);      // number of nonempty codimensions
	
	// [ this is repeated for each nonempty codimension
  fscanf(IN, "%d",&codim);      // codimension  (iterated for each codimension?)
  fscanf(IN, "%d",&num_components);  // number of components (is whatever)
	// ] end repeat
	
  fscanf(IN, "%d",&num_pts);    // number of points (should be 1)
	
	//and then a binary matrix indicating membership on which component
	//from the appendices of the book:
	//	% Binary matrix with one row per point, columns corresponding to components.
	//	%                0 if given point is not on given component;
	//	%                1 else .
	
	for(ii=0;ii<num_components;ii++)  // i don't this is correct if there is more than one nonempty codimension.
		//todo: check this iterator limit  !!!
	{
		fscanf(IN, "%d", &component_indicator);
		if (component_indicator==1) {  //then is on this component
			component_number = ii;
		}
		
	}
	fclose(IN);

	return component_number;
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
	
	
  //fprintf(OUT, "CONFIG\n");
  while ((ch = fgetc(IN)) != EOF)
    fprintf(OUT, "%c", ch);
  fclose(IN);
  if(tracktype==3)
    fprintf(OUT, "TrackType: %d;\nDeleteTempFiles: 1;\nEND;\nINPUT\n",tracktype);
  else
    fprintf(OUT, "TrackType: %d;\nDeleteTempFiles: 0;\nEND;\nINPUT\n",tracktype);
	
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






