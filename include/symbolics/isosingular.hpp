#ifndef _ISOSINGULAR_H
#define _ISOSINGULAR_H


/**
 \file isosingular.hpp

 \brief contains core of the functions for isosingular deflation, which are not found in Bertini.
 */


#include "io/fileops.hpp"
#include "symbolics/derivative_systems.hpp"


/**
 \brief main method for isosingular deflation

 Feed in a point with respect to which to deflate, and a system to deflate, and it will automatically deflate using Matlab.

 \return boolean integer of whether it was successful.  Successful==1, other==0.
 \param num_deflations the number of deflations it took.
 \param deflation_sequence The sequence of coranks of jacobian matrices.
 \param program_options The current state of Bertini_real
 \param inputFile the name of the starting file.
 \param witness_point_filename Name of file containing file with point in it, which we deflate.
 \param output_name The desired output file name.
 \param max_deflations The maximum number of allowable deflations iterations.
 */
int isosingular_deflation(int *num_deflations, int **deflation_sequence,
						  BertiniRealConfig & program_options,
						  boost::filesystem::path inputFile,
						  boost::filesystem::path witness_point_filename,
						  boost::filesystem::path output_name,
						  int max_deflations);



/**
\brief Create a Matlab .m file which calls `deflate_no_subst` from the matlab_codes folder.

this method of deflating prevents the substitution of subfunctions into the functions while deflating, which can be helpful in preventing expression swell.

\param OUT the file to write to.  ideally named 'matlab_deflate.m'.
\param filename The name of the input file to deflate
\param deflation_number The sequential number of deflations having been run.  Increment by 1 between calls to prevent duplicate symbols.
\param minorSize the size of the minors of the jacobian of which to take determinants.
\param degrees An array of integers containing the degrees of the functions.  Likely obtained from the bertini parser.
\param numFuncs The number of functions in the system being deflated.
*/
bool createMatlabDeflationNoSubst(std::string const& filename, int deflation_number, int minorSize, int* degrees, int numFuncs, boost::filesystem::path inputOutputName);




/**
 \brief Create a matlab .m file to perform symbolic determinants of minors.

 Write a .m file, which must be called by matlab, to perform symbolic deflation of a system.  Computes minors of the Jacobian matrix, and write them into a new system.

 \param OUT an open file to which to write
 \param numVars the number of variables
 \param vars The names of the variables
 \param lineVars the lines on which the variables occur.
 \param numConstants the number of constants
 \param consts the names of the constants
 \param lineConstants the lines on which the constants appear
 \param numFuncs the number of functions in the input file being deflated
 \param funcs the functions being deflated
 \param lineFuncs the lines on which the functions appear.
 \param IN open file from which to read.
 \param minorSize the size of the minor matrices to append determinants of.
 \param degrees the degrees of the functions
 \param deflation_number the integer index of the deflation iteration.
 */
bool createMatlabDeflation(int numVars, char **vars, int *lineVars,
						   int numConstants, char **consts, int *lineConstants,
						   int numFuncs, char **funcs, int *lineFuncs,
						   FILE *IN,
						   int minorSize, int *degrees, int deflation_number);



/**
 \brief Create a python .py file to perform symbolic determinants of minors.

 Write a .py file, which must be called by python, to perform symbolic deflation of a system.  Computes minors of the Jacobian matrix, and write them into a new system.

 \param OUT an open file to which to write
 \param numVars the number of variables
 \param vars The names of the variables
 \param lineVars the lines on which the variables occur.
 \param numConstants the number of constants
 \param consts the names of the constants
 \param lineConstants the lines on which the constants appear
 \param numFuncs the number of functions in the input file being deflated
 \param funcs the functions being deflated
 \param lineFuncs the lines on which the functions appear.
 \param IN open file from which to read.
 \param minorSize the size of the minor matrices to append determinants of.
 \param degrees the degrees of the functions
 \param deflation_number the integer index of the deflation iteration.
 */
bool createPythonDeflation(int numVars, char **vars, int *lineVars,
						   int numConstants, char **consts, int *lineConstants,
						   int numFuncs, char **funcs, int *lineFuncs,
						   FILE *IN,
						   int minorSize, int *degrees, int deflation_number);



/**
\brief turn the two files `deflation_polynomials` and `deflation_polynomials_declaration`, together with the original input file, into the final deflated input file.


*/
void DeflPolyDeclAndPolyToFinal(boost::filesystem::path inputOutputName, FILE* IN,
				int *declarations);

/**
 \brief setup input file for one deflation iteration

 \param declarations The numbers of various types of items appearing in the input file
 \param inputOutputName the name of the file
 \param matlab_command command to run Matlab
 \param nullSpaceDim used to calculate the size of the minor matrices.
 \param deflation_number integer counter
 */
void isosingular_deflation_iteration(int *declarations,
									 boost::filesystem::path inputOutputName,
									 BertiniRealConfig & program_options, int nullSpaceDim, int deflation_number);





/**
 \brief setup input file to test for stabilization of isosingular deflation

 \param outputFile The name of the file to write
 \param funcInput The name of the file containing the functions (etc) to parse out.
 \param configInput The name of the file containing the CONFIG part of the Bertini input file.
 */
void stabilization_input_file(boost::filesystem::path outputFile,
							  boost::filesystem::path funcInput,
							  boost::filesystem::path configInput);


/**
 \brief A sanity check for the numbers of various types of items in a Bertini input file.

 \param declarations The counters of stuffs, created by partitionParse
 */
void check_declarations(int *declarations);

#endif
