#include "dataToSet.hpp"


#define PREDETERMINEDLENGTH 4096

int main(int argC, char *args[]){
	
	
	std::string inputName, directoryName, outputName;
	
	
	std::ifstream IN;
	std::ofstream OUT;
	
	int num_vars,
	num_nonempty_codims,
	max_prec_used,
	multiplicity,
	component_number,
	deflations_needed,
	typeflag,
	current_codimension,
	num_points_this_dim,
	corank,
	numbertype,  // indicates double, multipleprecision, or rational.
	garbage,  //for reading in stuff which will be discarded.
	ii,jj,kk,mm,   // counters for loops
	num_rows_randomization, num_cols_randomization,
	num_linears;
	
	
	double condition_number,
	smallest_nonsing_value,
	largest_nonsing_value,
	real,imag,
	hom_variable_constant_r, hom_variable_constant_i;
	
	std::string tmpstr;
	std::stringstream converter;
	
	dataToSetStartup(argC, args, &inputName, &directoryName);  //sets inputName, directoryName
	
	//safely opens the input file
	open_input_file(inputName, IN);
	
	//remove the previous files with directoryName as basis.  will not remove files which start with a period.
	purge_previous_directory(const_cast<char *>( directoryName.c_str())); //this still works fine, after porting to c++
	
	//make directory to hold the data.
	mkdir(directoryName.c_str(),0777);
	
	
	//get the number of variables, and number of nonempty codimensions.
	
	getline(IN,tmpstr);
	converter << tmpstr;
	converter >> num_vars;
	converter.clear();
	converter.str("");
	
	
	getline(IN,tmpstr);
	converter << tmpstr;
	converter >> num_nonempty_codims;
	converter.clear();
	converter.str("");
	


	vec_d solution, prev_approx, dehomogenized_solution;
	init_vec_d(solution,num_vars);
	init_vec_d(prev_approx,num_vars);
	init_vec_d(dehomogenized_solution,num_vars-1);
	
	
	int component_counter[num_nonempty_codims]; // initialize these both to 0.
	memset( component_counter, 0, num_nonempty_codims*sizeof(int) );
	int codim_indicator[num_nonempty_codims];
	memset( codim_indicator, 0, num_nonempty_codims*sizeof(int) );

	

	std::map< std::pair < int, int >, witness_set> gathered_data;
	
	
	//enter the loop for the points.
	for (ii=0; ii<num_nonempty_codims; ii++) {
		
		getline(IN,tmpstr);
		converter << tmpstr;
		converter >> current_codimension;
		converter.clear();
		converter.str("");
		
		getline(IN,tmpstr);
		converter << tmpstr;
		converter >> num_points_this_dim;
		converter.clear();
		converter.str("");
		

		codim_indicator[ii] = current_codimension;
		
		for (jj=0; jj<num_points_this_dim; jj++) {
			
			getline(IN,tmpstr);
			converter << tmpstr;
			converter >> max_prec_used;
			converter.clear();
			converter.str("");
			
					
			std::vector < std::pair < double, double > > new_solution;
			new_solution.resize(num_vars);
			
			for (kk=0; kk<num_vars; kk++) {
				//here we get the actual solutions we are after.

				getline(IN,tmpstr);
				converter << tmpstr;
				converter >> real;
				converter >> imag;
				converter.clear();
				converter.str("");
				
				new_solution[kk] = std::pair < double, double > (real,imag);
			}
			

			//now have the solution in memory.
			
			
			//get the previous approximation
			getline(IN,tmpstr);
			converter << tmpstr;
			converter >> max_prec_used;
			converter.clear();
			converter.str("");
			
			for (kk=0; kk<num_vars; kk++) {
				// the previous approximations to the solutions
				getline(IN,tmpstr);
				converter << tmpstr;
				converter >> real;
				converter >> imag;
				converter.clear();
				converter.str("");
				prev_approx->coord[kk].r = real;
				prev_approx->coord[kk].i = imag;
			}
			
			getline(IN,tmpstr);
			converter << tmpstr;
			converter >> condition_number;
			converter.clear();
			converter.str("");
			
			getline(IN,tmpstr);
			converter << tmpstr;
			converter >> corank;
			converter.clear();
			converter.str("");
			
			getline(IN,tmpstr);
			converter << tmpstr;
			converter >> smallest_nonsing_value;
			converter.clear();
			converter.str("");
			
			getline(IN,tmpstr);
			converter << tmpstr;
			converter >> largest_nonsing_value;
			converter.clear();
			converter.str("");
			
			getline(IN,tmpstr);
			converter << tmpstr;
			converter >> typeflag;
			converter.clear();
			converter.str("");
			
			getline(IN,tmpstr);
			converter << tmpstr;
			converter >> multiplicity;
			converter.clear();
			converter.str("");
			
			getline(IN,tmpstr);
			converter << tmpstr;
			converter >> component_number;
			converter.clear();
			converter.str("");
			
			getline(IN,tmpstr);
			converter << tmpstr;
			converter >> deflations_needed;
			converter.clear();
			converter.str("");
			
			gathered_data[std::pair< int, int> (current_codimension, component_number)].add_raw_solution(current_codimension, component_number, new_solution);
			
			if (component_number+1 > component_counter[ii]) {
				component_counter[ii]  = component_number+1;  // keeps track of the number of components per dimension
			}
			
		}
		
		
	}//  re: first loop, to get points.
	
	


	getline(IN,tmpstr);
	converter << tmpstr;
	converter >> garbage;
	converter.clear();
	converter.str("");
	
	getline(IN,tmpstr);//waste a line

	getline(IN,tmpstr);
	converter << tmpstr;
	converter >> numbertype;
	converter.clear();
	converter.str("");
	

	
	

	
	for (mm=0; mm<num_nonempty_codims; mm++) {
		//BEGIN REPEATING BLOCK, one per nonempty codimension
		
		
		
		
		
		//initialize scope data members.
		std::vector < std::string > randomization_matrix,
				homogenization_matrix,
				homogenization_patch_eqn,
				linear_slice,
				patch_coefficents;
		std::string tmpstr;
		
		

		

		// the cursor must be at the correct position here, or things get messed up.
		getline(IN,tmpstr);
		converter << tmpstr;
		converter >> num_rows_randomization;
		converter >> num_cols_randomization;
		converter.clear();
		converter.str("");
		
		for (ii = 0; ii < num_rows_randomization*num_cols_randomization; ii++) {
			getline(IN,tmpstr);
			randomization_matrix.push_back(tmpstr);
		}
		

		
		
		
		//MATRIX W FOR HOMOGENIZATION
		//  same length as A, the randomization matrix
		int A[num_rows_randomization][num_cols_randomization];

		
		for (ii=0; ii<num_rows_randomization; ii++) {
			getline(IN,tmpstr);
			converter << tmpstr;
				
			for (jj=0; jj<num_cols_randomization; jj++) {

				converter >> A[ii][jj];
			}
			converter.clear();
			converter.str("");
		}

		
		getline(IN,tmpstr);//waste the blank line
		
		getline(IN,tmpstr);
		converter << tmpstr;
		converter >> garbage;
		converter.clear();
		converter.str("");
		
//		//VECTOR H FOR HOMOGENIZATION
		for (ii = 0; ii<garbage; ii++) {
			getline(IN,tmpstr);
			homogenization_patch_eqn.push_back(tmpstr);
		}


		
		
		
//		//   HOMVARCONST
		getline(IN,tmpstr);

		getline(IN,tmpstr);
		converter << tmpstr;
		converter >> hom_variable_constant_r;
		converter >> hom_variable_constant_i;
		converter.clear();
		converter.str("");

//		//MATRIX B FOR LINEAR SLICE COEFFICIENTS
		getline(IN,tmpstr);
		converter << tmpstr;
		converter >> num_linears;
		converter >> garbage;
		converter.clear();
		converter.str("");

//		// the cursor must be at the correct position here, or things get messed up.
		
		for (ii = 0; ii < num_linears; ii++) {
			std::vector< std::string >  tmpslice;
			for (jj=0; jj< num_vars; jj++) {
				getline(IN,tmpstr);
				tmpslice.push_back(tmpstr);
			}
			
			gathered_data[std::pair< int, int> (current_codimension, component_number)].add_linear_slice(current_codimension, component_number, tmpslice);
		}
		

		
		getline(IN,tmpstr);//waste blank line

		
//		// PATCH COEFFICIENTS
		for (ii = 0; ii<num_vars; ii++) {
			getline(IN,tmpstr);
			patch_coefficents.push_back(tmpstr);
		}
	
		//END REPEATED BLOCK (ONE FOR EACH NONEMPTY CODIM).

	}
	
	
	std::map < std::pair < int, int >, witness_set >::iterator iter;
	for (iter=gathered_data.begin(); iter!=gathered_data.end();  iter++) {
		iter->second.dehomogenize_and_write(directoryName);
	}
	
	
	std::cout << std::endl;
	
	return 0;
}  //re: main









void make_specific_output_name(std::string *outputName, std::string directoryName, int current_codimension,int component_number){
	
	std::stringstream converter;
	converter << directoryName << "/" << "codim_" << current_codimension << "_comp_" << component_number;
	*outputName = converter.str();
	return;
}



int dataToSetStartup(int argC, char *args[], std::string *inputName, std::string *directoryName)
/***************************************************************\
 * USAGE:    prepares the variables inputname and directoryName
 *      for use later in the program
 * ARGUMENTS:                                                    *
 * RETURN VALUES:           an integer                           *
 * NOTES:                                                        *
 \***************************************************************/
{
	
	
	std::cout << "\n\n\nthis is data2set, a utility for turning witness_data into " << std::endl <<
		"witness_set files, to be used with BertiniReal\n\n\n" << std::endl;
	
	
	// setup inputName & directoryName
	if (argC >= 2 && args[1] != NULL)
	{ // inputName is args[1]
		
		*inputName = args[1];
		
		// setup directoryName
		if (argC >= 3 && args[2] != NULL)
		{ // directoryName is args[2]
			*directoryName = args[2];
		}
		else
		{ // default to 'start'
			*directoryName = "witness_set";
		}
	}
	else
	{ // default to 'input' & 'start'
		*inputName = "witness_data";
		*directoryName = "witness_set";
	}
	
	return 0;
}




//opens an input file.  exits the program if the file DNE
int open_input_file(std::string inputName, std::ifstream & IN){
	
	IN.clear();
	IN.open(inputName.c_str());
	
	if (!IN.is_open()) {
		std::cerr << "unable to open input file " << inputName << std::endl;
		bexit(ERROR_FILE_NOT_EXIST);
	}
	else{
//		std::cout << "successfully opened input " << inputName << std::endl;
	}
	
	return 0;
}


//opens a file to write.  overwrites previous file if exists.
int open_output_file(std::string outputName, std::ofstream & OUT){
	
	OUT.clear();
	OUT.open(outputName.c_str());
	
	if (!OUT.is_open()) {
		std::cerr << "unable to open file " << outputName << " to write!\n";		bexit(ERROR_FILE_NOT_EXIST);
	}
	else{
	}
	
	return 0;
}


//opens a file to write, in append mode.  starts new file if necessary.
int open_output_file_append(std::string outputName, std::ofstream & OUT){
	
	OUT.clear();
	OUT.open(outputName.c_str(),std::ios::app);
	
	
	if (!OUT.is_open()) {
		std::cerr << "unable to open file " << outputName << " to append!\n";
		bexit(ERROR_FILE_NOT_EXIST);
	}
	else{
	}
	
	return 0;
}



//will remove all files which do not start with a . from directory.
void purge_previous_directory(char *directoryName)
{
  DIR *dp;
  struct dirent *ep;
	char tempname[PREDETERMINEDLENGTH];
	
	
  dp = opendir (directoryName);
  if (dp != NULL)
	{
		while ( (ep = readdir (dp)) )
			if (ep->d_name[0] != '.') {
				sprintf(tempname,"%s/%s",directoryName, ep->d_name);
				remove(tempname);
			}
		
		
		(void) closedir (dp);
//		std::cout << "\ndone deleting previous files\n\n\n";
		
	}
  else
	{
//    printf ("no previous directory to delete\n");
	}
	
	
	return;
}




