#ifndef _FILEOPS_H
#define _FILEOPS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include <time.h>
#include <float.h>
#include <limits.h>

#include <gmp.h>
#include <mpfr.h>
#include <mpf2mpfr.h>

#include <sys/stat.h>         /* declare the 'stat' structure */
#include <sys/types.h>


#include <dirent.h>

#include <string>
#include <sstream>
#include <set>



#include "bertini_headers.hpp"


extern "C" {
#include "partitionParse.h"
}


#include "boost/filesystem.hpp"


int partition_parse(int **declarations,
					boost::filesystem::path input_filename,
					boost::filesystem::path functions_filename,
					boost::filesystem::path config_filename,
					int not_sc_flag);


/**
 renames the arr.out deg.out num.out config and preproc_data files to *.bak.  these files restored by restore_bertini_files_dotbak().
 */
void rename_bertini_files_dotbak();


/**
 restores the .bak files which were renamed to *.bak by move_bertini_files_dotbak
 */
void restore_bertini_files_dotbak();

/**
 purges an entire directory of all files except those which start with a period (.)
 */
void purge_previous_directory(char *directoryName);

/**
 opens a file to read, giving fatal error if cannot.
 */
FILE *safe_fopen_read(boost::filesystem::path filename);

/**
 opens a file to write, giving fatal error if cannot.
 */
FILE *safe_fopen_write(boost::filesystem::path filename);

/**
 opens a file to append, giving fatal error if cannot.
 */
FILE *safe_fopen_append(boost::filesystem::path filename);

/**
 copies a file character by character.
 */


void copyfile(boost::filesystem::path input_file, boost::filesystem::path OUTfile);

void copyfile(FILE *IN,FILE *OUT);

void read_matrix(boost::filesystem::path INfile, mat_mp matrix);


void br_exit(int errorCode);

void deliberate_segfault();


/**
 Have the user input a value untl it's an integer, return that value.
 \return int - The integer the user inputted.
 */
int getInteger();




/**
 Parse a string that has an integer value in string form.
 \param text  - The integer value as a string.
 \param results - The value to set as a string.
 \return bool - A boolean to indicate whether the parsing was successful or not.
 */
bool parseInteger( std::string const& text, int& results );

/**
 Display a menu option to the user and ask for an integer input within the specified range.
 \param display_string - The menu as a string.
 \param min_value - The minimum value allowed.
 \param max_value - The maximum value allowed.
 \return int - The integer the user specified.
 */
int get_int_choice(std::string display_string,int min_value,int max_value);


/**
 Display a menu option to the user and ask for an integer input within the specified range.
 \param display_string - The menu as a string.
 \param valid_values.  A std::set of valid integer values.  all others will be rejected.
 \return int - The integer the user specified.
 */
int get_int_choice(std::string display_string,const std::set<int> & valid_values);



#endif

