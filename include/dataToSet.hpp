
#include <dirent.h>

#include <sys/stat.h>

#include <iostream>
#include <ios>
#include <string>
#include <fstream>

#include <vector>
#include <map>
#include <sstream>

#include "mpreal.h"  // <---- this software graciously free to use for free projects, which this is.


#ifndef _DATA_TO_SET_H
#define _DATA_TO_SET_H

#define DTSVERSION 1e-16
extern "C" {
#include "polysolve.h"
}

extern "C" {
#include "data_type.h"
}


void purge_previous_directory(char *directoryName);

void make_specific_output_name(std::string *specificOutputName,
															 std::string directoryName,
															 int component_number,
															 int current_codimension);


int dataToSetStartup(int argC, char *args[], std::string *inputName, std::string *outputName);

int open_input_file(std::string inputName, std::ifstream & IN);  //perhaps these ought to be moved into a more general file?
int open_output_file(std::string outputName, std::ofstream & OUT);
int open_output_file_append(std::string outputName, std::ofstream & OUT);


/**
 *  @class "witness_set"
 *
 * contains a parsed witness_data file from Bertini's tracktype: 1 solver.
 *
 * \brief parsed witness_data file from Bertini's tracktype: 1 solver
 **/
class witness_set {
	
public:
	
	witness_set(){
		this->component_number = -10;
		this->codimension = -10;
	}
	
	
	
	void add_raw_solution(int _dimmy, int _pony, std::vector < std::pair < std::string, std::string > > new_solution){
		this->solutions_raw.push_back(new_solution);
		this->component_number = _pony;
		this->codimension = _dimmy;
	};
	
	void add_linear_slice( int _dimmy, int _pony, std::vector < std::string > new_slice, int _numtype){
		this->linear_slice.push_back(new_slice);
		this->component_number = _pony;
		this->codimension = _dimmy;
		this->numtype = _numtype;
	};
	
	void add_patch_equation( int _dimmy, int _pony, std::vector < std::string > new_patch, int _numtype){
		this->patch_equation.push_back(new_patch);
		this->component_number = _pony;
		this->codimension = _dimmy;
		this->numtype = _numtype;
	};
	
	void write_to_file(std::string directoryName, bool dehomogenize);
	
	int degree(){
		return this->solutions_raw.size();
	};
	
	
private:
	int component_number;
	int codimension;
	int numtype;
	std::vector < std::vector < std::string > > linear_slice;
	std::vector < std::vector < std::string > > patch_equation;
	std::vector < std::vector < std::pair < std::string, std::string > > > solutions_raw;
};


typedef std::pair< mpfr::mpreal, mpfr::mpreal > complex_number;

//divides z by w.
complex_number complex_divide(complex_number z, complex_number w){
	complex_number result;
	
	mpfr::mpreal denomenator = w.first*w.first + w.second*w.second;
	
	result.first = z.first*w.first + z.second*w.second;
	result.second = z.second*w.first - z.first*w.second;
	
	result.first = result.first/denomenator;
	result.second = result.second/denomenator;
	
	return result;
}


void witness_set::write_to_file(std::string directoryName, bool dehomogenize){
	
	using mpfr::mpreal;

	int digits = 50;
	std::stringstream converter;
	
	mpreal::set_default_prec(mpfr::digits2bits(digits));
	std::cout.setf (std::ios_base::scientific);
	

	
	std::string specificOutputName;
	std::ofstream OUT;
	
	make_specific_output_name(&specificOutputName, directoryName, this->codimension, this->component_number);
	
	open_output_file(specificOutputName, OUT);
	OUT.setf(std::ios_base::scientific);
	

	std::cout << "writing data to " << specificOutputName <<  std::endl;

	
	
	OUT << this->solutions_raw.size() << " " << this->codimension << " " << this->component_number;
	OUT << std::endl << std::endl;  // print the number of solutions, and the number of variables.
	
	
	for (int ii=0; ii< int(this->solutions_raw.size()) ; ii++)
	{
		complex_number coordinate;
		if (dehomogenize==true){
			complex_number dehomogenizing_variable;
			
			converter << solutions_raw[ii][0].first;
			converter >> dehomogenizing_variable.first;
			converter.clear();
			converter.str("");
			
			converter << solutions_raw[ii][0].second;
			converter >> dehomogenizing_variable.second;
			converter.clear();
			converter.str("");
			

			for (int jj = 1; jj<int(solutions_raw[ii].size()); jj++)
			{
				complex_number result;
				
				converter << solutions_raw[ii][jj].first;
				converter >> coordinate.first;
				converter.clear();
				converter.str("");
				
				converter << solutions_raw[ii][jj].second;
				converter >> coordinate.second;
				converter.clear();
				converter.str("");
								
				result = complex_divide(coordinate,dehomogenizing_variable);
				
				digits = std::max(solutions_raw[ii][jj].first.size(),solutions_raw[ii][jj].second.size());
				OUT.precision(digits);
				OUT << result.first << " " << result.second << std::endl;
			}
		}
		else // dehomogenize == false
		{
			for (int jj = 0; jj<int(solutions_raw[ii].size()); jj++)
			{
				//write the solutions straight as the strings obtained in the read portion of the program.
				OUT << solutions_raw[ii][jj].first << " " << solutions_raw[ii][jj].second << std::endl;
			}
		}
		OUT << std::endl;
	}
		


	
	OUT << this->linear_slice.size();
	if (this->linear_slice.size()>0)
	{
		OUT << " " << this->linear_slice[0].size();
	}
	else
	{
		OUT << " 0";
	}
	OUT << std::endl << std::endl;
	
	std::string tmpstr_real, tmpstr_imag;
	
	for (int ii=0; ii< int(this->linear_slice.size()) ; ii++) {
		for (int jj = 0; jj<int(linear_slice[ii].size()); jj++) {
		
			
			if (this->numtype==2)
			{

				std::string convertme_numerator, convertme_denomenator;
				std::string searchforme = "/";
				
				
				size_t found;
				std::string crap;
				
				
				mpreal numerator, denomenator;
				
				converter << linear_slice[ii][jj];
				converter >> tmpstr_real;
				converter >> tmpstr_imag;
				converter.clear();
				converter.str("");
				
				found=tmpstr_real.find(searchforme);
				
				
				convertme_numerator.resize(found);
				convertme_denomenator.resize(tmpstr_real.size()-found-1);
				for (int mm = 0; mm<int(found)+1; mm++) {
					convertme_numerator[mm] = tmpstr_real[mm];
				}
				
				for (int mm = found+1; mm<int(tmpstr_real.size()); mm++) {
					convertme_denomenator[mm-found-1] = tmpstr_real[mm];
				}
				
				
				digits = std::max(convertme_denomenator.size(),convertme_numerator.size());

				
				
				
				converter << convertme_numerator;
				converter >> numerator;
				converter.clear();
				converter.str("");
				
				converter << convertme_denomenator;
				converter >> denomenator;
				converter.clear();
				converter.str("");
				
				OUT.precision(digits);
				OUT << numerator/denomenator << " ";
				

				
				
				found=tmpstr_imag.find(searchforme);
				convertme_numerator.resize(found);
				convertme_denomenator.resize(tmpstr_imag.size()-found-1);
				for (int mm = 0; mm<int(found)+1; mm++) {
					convertme_numerator[mm] = tmpstr_imag[mm];
				}
				
				for (int mm = found+1; mm<int(tmpstr_imag.size()); mm++) {
					convertme_denomenator[mm-found-1] = tmpstr_imag[mm];
				}
				

				digits = std::max(convertme_denomenator.size(),convertme_numerator.size());

				converter << convertme_numerator;
				converter >> numerator;
				converter.clear();
				converter.str("");
				
				converter << convertme_denomenator;
				converter >> denomenator;
				converter.clear();
				converter.str("");
				

				OUT.precision(digits);
				OUT << numerator/denomenator << std::endl;
				
				
				
			}
			else
			{ // simply write the strings, no need to convert data types or set precision.
				converter << linear_slice[ii][jj];
				converter >> tmpstr_real;
				converter >> tmpstr_imag;
				converter.clear();
				converter.str("");
				OUT << tmpstr_real << " " << tmpstr_imag << std::endl;
			}
			
			
		} // re: for jj
		OUT << std::endl; //  write the line which separates linears
	} // re: for ii
	
	
	/////////
	//   print the patch equations
	/////////
	
	
	OUT << this->patch_equation.size();
	if (this->patch_equation.size()>0)
	{
		OUT << " " << this->patch_equation[0].size();
	}
	else
	{
		OUT << " 0";
	}
	OUT << std::endl << std::endl;
	
	for (int ii=0; ii< int(this->patch_equation.size()) ; ii++) {
		for (int jj = 0; jj<int(patch_equation[ii].size()); jj++) {
			

			std::string convertme_numerator, convertme_denomenator;
			std::string searchforme = "/";
			
			
			size_t found;
			std::string crap;
			
			
			mpreal numerator, denomenator;
			int digits;
			
			
			if (this->numtype==2)
			{
				
				converter << patch_equation[ii][jj];
				converter >> tmpstr_real;
				converter >> tmpstr_imag;
				converter.clear();
				converter.str("");
				
				found=tmpstr_real.find(searchforme);
				
				
				convertme_numerator.resize(found);
				convertme_denomenator.resize(tmpstr_real.size()-found-1);
				for (int mm = 0; mm<int(found)+1; mm++) {
					convertme_numerator[mm] = tmpstr_real[mm];
				}
				
				for (int mm = found+1; mm<int(tmpstr_real.size()); mm++) {
					convertme_denomenator[mm-found-1] = tmpstr_real[mm];
				}
				
				
				digits = std::max(convertme_denomenator.size(),convertme_numerator.size());
				mpreal::set_default_prec(mpfr::digits2bits(digits));
				std::cout.setf (std::ios_base::scientific);
				
				
				
				converter << convertme_numerator;
				converter >> numerator;
				converter.clear();
				converter.str("");
				
				converter << convertme_denomenator;
				converter >> denomenator;
				converter.clear();
				converter.str("");
				

				OUT.precision(digits);
				OUT << numerator/denomenator << " ";
				
				
				
				
				found=tmpstr_imag.find(searchforme);
				convertme_numerator.resize(found);
				convertme_denomenator.resize(tmpstr_imag.size()-found-1);
				for (int mm = 0; mm<int(found)+1; mm++) {
					convertme_numerator[mm] = tmpstr_imag[mm];
				}
				
				for (int mm = found+1; mm<int(tmpstr_imag.size()); mm++) {
					convertme_denomenator[mm-found-1] = tmpstr_imag[mm];
				}
				
	
				digits = std::max(convertme_denomenator.size(),convertme_numerator.size());
				mpreal::set_default_prec(mpfr::digits2bits(digits));

				converter << convertme_numerator;
				converter >> numerator;
				converter.clear();
				converter.str("");
				
				converter << convertme_denomenator;
				converter >> denomenator;
				converter.clear();
				converter.str("");
				
				
				OUT.precision(digits);
				OUT << numerator/denomenator << std::endl;
				
				
				
			}
			else
			{ // simply write the strings, no need to convert data types or set precision.
				converter << patch_equation[ii][jj];
				converter >> tmpstr_real;
				converter >> tmpstr_imag;
				converter.clear();
				converter.str("");
				OUT << tmpstr_real << " " << tmpstr_imag << std::endl;
			}
			
			
		} // re: for jj
		OUT << std::endl; //  write the line which separates linears
	} // re: for ii
	
	
	//finally close the file
	OUT.close();
	
}; //re: dehomogenize_and_write




void write_summary(std::map< std::pair < int, int >, witness_set> gathered_data,
									 int component_counter[],
									 int codim_indicator[],
									 int num_nonempty_codims,
									 const std::string directoryName);







#endif



