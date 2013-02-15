
#include <stdio.h>
#include <stdlib.h>
//#include <stdarg.h>
//#include <string.h>
//#include <math.h>
//#include <gmp.h>
//#include <time.h>
//#include <float.h>
//#include <limits.h>
//#include <sys/stat.h>
#include <dirent.h>
#include <mpfr.h>
#include <mpf2mpfr.h>

#include <errno.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <list>

#include <iostream>
#include <ios>
#include <string>
#include <fstream>
#include <unistd.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <sstream>
#include <cmath>

#ifndef _DATA_TO_SET_H
#define _DATA_TO_SET_H

extern "C" {
	#include "polysolve.h"
}
//extern "C" {
//	#include "memory.h"
//}
//
extern "C" {
	#include "data_type.h"
}

//void dehomogenize(vec_d dehomogenized_solution, vec_d solution, int num_vars);

void purge_previous_directory(char *directoryName);

void make_specific_output_name(std::string *specificOutputName,
															 std::string directoryName,
															 int component_number,
															 int current_codimension);


int dataToSetStartup(int argC, char *args[], std::string *inputName, std::string *outputName);

int open_input_file(std::string inputName, std::ifstream & IN);  //perhaps these ought to be moved into a more general file?
int open_output_file(std::string outputName, std::ofstream & OUT);
int open_output_file_append(std::string outputName, std::ofstream & OUT);


class witness_set {
	
public:
	
	witness_set(){
		this->component_number = -10;
		this->codimension = -10;
	}
	

	
	void add_raw_solution(int _dimmy, int _pony, std::vector < std::pair < double, double > > new_solution){
		this->solutions_raw.push_back(new_solution);
		this->component_number = _pony;
		this->codimension = _dimmy;
	};
	
	void add_linear_slice( int _dimmy, int _pony, std::vector < std::string > new_slice){
		this->linear_slice.push_back(new_slice);
		this->component_number = _pony;
		this->codimension = _dimmy;
	};
	
	
	void dehomogenize_and_write(std::string directoryName);
	

	

private:
	int component_number;
	int codimension;
	std::vector < std::vector < std::string > > linear_slice;
	std::vector < std::vector < std::pair < double, double > > > solutions_raw;
};

void witness_set::dehomogenize_and_write(std::string directoryName){
	
	std::string specificOutputName;
	std::ofstream OUT;
	
	make_specific_output_name(&specificOutputName, directoryName, this->codimension, this->component_number);
	
	open_output_file(specificOutputName, OUT);
	
	std::cout << "writing data to " << specificOutputName <<  std::endl;
	
	OUT << this->solutions_raw.size() << std::endl << std::endl;
	
	OUT.precision(16);
	
	for (int ii=0; ii< int(this->solutions_raw.size()) ; ii++) {

		_comp_d a,b,res;
		for (int jj = 0; jj<int(solutions_raw[ii].size())-1; jj++) {
			a.r = solutions_raw[ii][jj+1].first;
			a.i = solutions_raw[ii][jj+1].second;
			b.r = solutions_raw[ii][0].first;
			b.i = solutions_raw[ii][0].second;
			div_d(&res, &a , &b);
			
			OUT << res.r << " " << res.i << std::endl;
		}
		OUT << std::endl;

	}
	
	OUT << this->linear_slice.size() << std::endl << std::endl;
	
	for (int ii=0; ii< int(this->linear_slice.size()) ; ii++) {
		for (int jj = 0; jj<int(linear_slice[ii].size()); jj++) {

			OUT << linear_slice[ii][jj] << std::endl;
		}
		OUT << std::endl;
		
	}
	
	
	OUT.close();
}; //re: dehomogenize_and_write
//
//class datatracker {
//
//public:
//	
//
//	
//	datatracker(int _num_vars){
//		this->num_vars = _num_vars;
//	};
//	
//
//	
//	void add_raw_solution(int _codimension, int _component_number, std::vector < std::pair < double, double > > new_solution){
//		
//		
//		
//		
//		this->mega_data[std::pair< int, int> (_codimension,_component_number) ].num_points++;
//		this->mega_data[_codimension,_component_number].solutions_raw.push_back(new_solution);
//	};
//	
//	vo
//	
//	
//private:
//	int num_points;
//	int codimension;
//	int component_number;
//	int num_vars;
//	
//	std::map < std::pair < int, int > , witness_set> mega_data;
//	
//
//};


#endif



