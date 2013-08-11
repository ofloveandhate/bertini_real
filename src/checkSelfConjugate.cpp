#include "checkSelfConjugate.hpp"







int checkSelfConjugate(witness_set & W,
                       int           num_vars,
											 BR_configuration & program_options,
											 boost::filesystem::path input_file)
/***************************************************************\
 * USAGE: check if component is self conjugate                  *
 * ARGUMENTS: witness set, # of variables and name of input file*
 * RETURN VALUES: 1-self conjugate; 0-non self conjugate   *
 * NOTES:                                                        *
 \***************************************************************/
{
  int *declarations = NULL;
	std::string bertini_command=program_options.bertini_command;
	bertini_command.append(" input_membership_test ");
	bertini_command.append(program_options.stifle_text);
  FILE *IN = NULL;
	

	
	
  // setup input file
  partition_parse(&declarations, input_file, "func_input_real", "config_real",0); // the 0 means not self conjugate

	
	
	//check existence of the required witness_data file.
	IN = safe_fopen_read("witness_data");
	fclose(IN);
	
	
	
	
	
	

	//only need to do this once.  we put the point and its conjugate into the same member points file and run the membership test simultaneously with one bertini call.
  membership_test_input_file("input_membership_test", "func_input_real", "config_real",3);
	
	
	//setup  member_points file, including both the first witness point, and its complex-conjugate
	write_member_points_sc(W.pts_mp[0]);//,fmt
	
	// Do membership test
	std::cout << "*\n" << program_options.bertini_command << "\n*" << std::endl;
  system(bertini_command.c_str());
	
	int *component_numbers;
	component_numbers = (int *)br_malloc(2*sizeof(int));
	
	read_incidence_matrix(component_numbers);

  free(declarations);
	
	
	
  // delete temporary files

	
	if (component_numbers[0]==component_numbers[1]) {
		printf("component IS self conjugate\n");
		free(component_numbers);
		return 1;
	}
	else
	{
		printf("component is NOT self conjugate\n");
		free(component_numbers);
		return 0;
	}

	
}





int get_incidence_number(const witness_set & W,
												 int           num_vars,
												 BR_configuration & program_options,
												 boost::filesystem::path input_file)
/***************************************************************\
 * USAGE: check if component is self conjugate                  *
 * ARGUMENTS: witness set, # of variables and name of input file*
 * RETURN VALUES: 1-self conjugate; 0-non self conjugate   *
 * NOTES:                                                        *
 \***************************************************************/
{
  
	std::string system_command=program_options.bertini_command;
	system_command.append(" input_membership_test ");
	system_command.append(program_options.stifle_text);
  
	
  // setup input file
	int *declarations = NULL;
  partition_parse(&declarations, input_file, "func_input_real", "config_real",0); // the 0 means not self conjugate
	free(declarations);
	
	
	//check existence of the required witness_data file.
	FILE *IN = NULL;
	IN = safe_fopen_read("witness_data"); fclose(IN);
	

	
	//only need to do this once.  we put the point and its conjugate into the same member points file and run the membership test simultaneously with one bertini call.
  membership_test_input_file("input_membership_test", "func_input_real", "config_real",3);
	
	//setup  member_points file, including both the first witness point, and its complex-conjugate
	write_member_points_singlept(W.pts_mp[0]);//,fmt
	// Do membership test
	std::cout << "*\n" << system_command << "\n*" << std::endl;
  system(system_command.c_str());
	
	int component_number;
	
	read_incidence_matrix(&component_number);
	
	return component_number;
}



int write_member_points_singlept(point_mp point_to_write){
	FILE *OUT = NULL;
	int ii;
	
	
	remove("member_points");
	OUT = safe_fopen_write("member_points");

	
	vec_mp result;
	init_vec_mp(result,0);
	dehomogenize(&result,point_to_write);
	
	fprintf(OUT,"1\n\n");
	for(ii=0;ii<result->size;ii++){
		print_mp(OUT,0,&result->coord[ii]);
		fprintf(OUT,"\n");
	}
//		fprintf(OUT, fmt, result->coord[ii].r, result->coord[ii].i);
	
	
	fclose(OUT);
	
	
	clear_vec_mp(result);
	
	return 0;
}



int write_member_points_sc(point_mp point_to_write){
	FILE *OUT = NULL;
	int ii;
	
	
	remove("member_points");
	OUT = safe_fopen_write("member_points");

	
	vec_mp result;
	init_vec_mp(result,0);
	dehomogenize(&result,point_to_write);
	
	fprintf(OUT,"2\n\n");
	for(ii=0;ii<result->size;ii++)
	{
//		fprintf(OUT, fmt, result->coord[ii].r, result->coord[ii].i);
		print_mp(OUT,0,&result->coord[ii]);  fprintf(OUT,"\n");
	}
	
	
	comp_mp temp; init_mp(temp);
	fprintf(OUT,"\n");
	for(ii=0;ii<result->size;ii++)
	{
		conjugate_mp(temp, &result->coord[ii]);
		print_mp(OUT,0,temp);  fprintf(OUT,"\n");
//		fprintf(OUT, fmt, result->coord[ii].r, -result->coord[ii].i);
	}
	fclose(OUT);
	
	clear_vec_mp(result);  clear_mp(temp);
	
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
//	printf("reading incidence for %d pts\n",num_pts);

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
			printf("it appears the candidate point lies on NO COMPONENT.\n");
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
		component_number=0;
		for(ii=0;ii<total_num_components;ii++)  // i don't know if this is correct if there is more than one nonempty codimension.
		{
			fscanf(IN, "%d", &component_indicator);
			if (component_indicator==1 && given_incidence_number==ii) {  //then is on this component
				component_number++;
			}
		}
		if (component_number==0) {
//			printf("test point did not lie on any component at all.\n");
			//			exit(-1);
		}
		component_numbers[jj]=component_number;
	}
	
	fclose(IN);

	
	return;
}





void membership_test_input_file(boost::filesystem::path outputFile,
                                boost::filesystem::path funcInput,
                                boost::filesystem::path configInput,
                                int  tracktype)
/***************************************************************\
 * USAGE: setup input file to membership test                    *
 * ARGUMENTS: name of output file, function & configuration input*
 * RETURN VALUES: none                                           *
 * NOTES:                                                        *
 \***************************************************************/
{
  char ch;
  FILE *OUT = safe_fopen_write(outputFile), *IN = NULL;

	
  // setup configurations in OUT
  fprintf(OUT, "CONFIG\n");
  IN = safe_fopen_read(configInput);

	
	while ((ch = fgetc(IN)) != EOF)
    fprintf(OUT, "%c", ch);
  fclose(IN);
  if(tracktype==3)
    fprintf(OUT, "TrackType: %d;\nDeleteTempFiles: 1;\nEND;\nINPUT\n",tracktype);
  else
    fprintf(OUT, "TrackType: %d;\nDeleteTempFiles: 1;\nEND;\nINPUT\n",tracktype);
	
  // setup system in OUT
  IN = safe_fopen_read(funcInput);

  while ((ch = fgetc(IN)) != EOF)
    fprintf(OUT, "%c", ch);
  fclose(IN);
  fprintf(OUT, "END;\n");
  fclose(OUT);
	
  return;
}




