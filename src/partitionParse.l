%{
/* This will take an input file and partition it into a configurations file and a function file */

	#include "polysolve.h"


// settings to switch between config & function
int configFoundPartition = 0, endofConfigPartition = 1, endofFuncPartition = 0;
// number of times the statement appears
int num_variable_group = 0, num_hom_variable_group = 0, num_variable = 0, num_pathvariable = 0, num_parameter = 0, num_definedSubfunction = 0, num_random = 0, num_random_real = 0, num_constant = 0, num_function = 0;
// files for config & function
FILE *CONFIGFILEPARTITION = NULL, *FUNCFILEPARTITION = NULL;

/* %option nounput - do not warn about unput */
%}
%option nounput 

config "CONFIG"
input "INPUT"
end ("END"[ \t\r\n]+)
ends "END;"
homvargp "hom_variable_group"
vargp "variable_group"
pathvar "pathvariable"
var "variable"
param "parameter"
const "constant"
func "function"
subfunc "definedSubfunction"
rand "random"
randReal "random_real"
number ([0-9]+[\.]?([0-9]*)?([eE][-+]?[0-9]+)?)|([0-9]*\.[0-9]+([eE][-+]?[0-9]+)?) 
I [I]
name [a-zA-Z][a-zA-Z0-9_\[\]]*
comments [%].+[\n]
white [ \t\r\n]
eol [;]
colon [:]
other .

%%

{end}|{ends} { // end of something
             if (!endofConfigPartition)
             { // end config
               endofConfigPartition = 1;
             }
             else if (!endofFuncPartition)
             { // end func
               endofFuncPartition = 1;
             }
           }
{comments}|{white} { // ignore comments and whitespace
           }
{config}   { // start of config found - print there until reach END or INPUT
             configFoundPartition = 1;
             endofConfigPartition = 0;
           }
{input}    { // start of input found - print there until reach END or EOF
             endofConfigPartition = 1;
             endofFuncPartition = 0;
           }
{param}    { // declare parameter
             if (endofConfigPartition && !endofFuncPartition)
             { // print to FUNCFILEPARTITION
               num_parameter++;
               fprintf(FUNCFILEPARTITION, "%s ", partitionParsetext);
             }
           }
{homvargp} { // declare hom_variable_group
             if (endofConfigPartition && !endofFuncPartition)
             { // print to FUNCFILEPARTITION
               num_hom_variable_group++;
               fprintf(FUNCFILEPARTITION, "%s ", partitionParsetext);
             }
           }
{vargp}    { // declare variable_group
             if (endofConfigPartition && !endofFuncPartition)
             { // print to FUNCFILEPARTITION
               num_variable_group++;
               fprintf(FUNCFILEPARTITION, "%s ", partitionParsetext);
             }
           }
{var}      { // declare variable
             if (endofConfigPartition && !endofFuncPartition)
             { // print to FUNCFILEPARTITION
               num_variable++;
               fprintf(FUNCFILEPARTITION, "%s ", partitionParsetext);
             }
           }
{pathvar}  { // declare pathvariable
             if (endofConfigPartition && !endofFuncPartition)
             { // print to FUNCFILEPARTITION
               num_pathvariable++;
               fprintf(FUNCFILEPARTITION, "%s ", partitionParsetext);
             }
           }
{const}    { // declare constant
             if (endofConfigPartition && !endofFuncPartition)
             { // print to FUNCFILEPARTITION
               num_constant++;
               fprintf(FUNCFILEPARTITION, "%s ", partitionParsetext);
             }
           }
{func}     { // declare function
             if (endofConfigPartition && !endofFuncPartition)
             { // print to FUNCFILEPARTITION
               num_function++;
               fprintf(FUNCFILEPARTITION, "%s ", partitionParsetext);
             }
           }
{subfunc}  { // declare definedSubfunction
             if (endofConfigPartition && !endofFuncPartition)
             { // print to FUNCFILEPARTITION
               num_definedSubfunction++;
               fprintf(FUNCFILEPARTITION, "%s ", partitionParsetext);
             }
           }
{rand}     { // declare random
             if (endofConfigPartition && !endofFuncPartition)
             { // print to FUNCFILEPARTITION
               num_random++;
               fprintf(FUNCFILEPARTITION, "%s ", partitionParsetext);
             }
           }
{randReal} { // declare random_real
             if (endofConfigPartition && !endofFuncPartition)
             { // print to FUNCFILEPARTITION
               num_random_real++;
               fprintf(FUNCFILEPARTITION, "%s ", partitionParsetext);
             }
           }
{colon}    { // a colon
             if (!endofConfigPartition)
             { // add a space after a colon in configurations
               fprintf(CONFIGFILEPARTITION, ": ");
             }
           }
{eol}      { // end of line
             if (!endofConfigPartition)
             { // print to CONFIGFILEPARTITION
  	       fprintf(CONFIGFILEPARTITION, ";\n");
             }
             else if (!endofFuncPartition)
             { // print to FUNCFILEPARTITION
               fprintf(FUNCFILEPARTITION, ";\n");
             }
           }
{number}|{I}|{name}|{other}   { // other symbol
             if (!endofConfigPartition)
             { // print to CONFIGFILEPARTITION
               fprintf(CONFIGFILEPARTITION, "%s", partitionParsetext);
             }
             else if (!endofFuncPartition)
             { // print to FUNCFILEPARTITION
               fprintf(FUNCFILEPARTITION, "%s", partitionParsetext);
             }
           }
%%

// parse a file
int partitionParse(int **declarations, FILE *fp, char *funcName, char *configName)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: partitions fp into function & configurations               *
\***************************************************************/
{
  int rV = 0;

  // setup FUNCFILEPARTITION & CONFIGFILEPARTITION
  FUNCFILEPARTITION = fopen(funcName, "w");
  CONFIGFILEPARTITION = fopen(configName, "w");
  if (FUNCFILEPARTITION == NULL)
  {
    printf("ERROR: The name %s is not valid!\n\n", funcName);
    bexit(ERROR_INPUT_SYNTAX);
  }
  if (CONFIGFILEPARTITION == NULL)
  {
    printf("ERROR: The name %s is not valid!\n\n", configName);
    bexit(ERROR_INPUT_SYNTAX);
  }

  // parse the file
  partitionParsein = fp;
  rV = partitionParselex();

  // print the number of declarations to a file
  *declarations = (int *)bmalloc(10 * sizeof(int));
  (*declarations)[0] = num_variable_group;
  (*declarations)[1] = num_hom_variable_group;
  (*declarations)[2] = num_variable;
  (*declarations)[3] = num_pathvariable;
  (*declarations)[4] = num_parameter;
  (*declarations)[5] = num_definedSubfunction;
  (*declarations)[6] = num_random;
  (*declarations)[7] = num_random_real;
  (*declarations)[8] = num_constant;
  (*declarations)[9] = num_function;

  // close FUNCFILEPARTITION & CONFIGFILEPARTITION
  if (!endofFuncPartition)
    fprintf(FUNCFILEPARTITION, "\nEND;\n\n");

  fclose(FUNCFILEPARTITION);
  fclose(CONFIGFILEPARTITION);
  FUNCFILEPARTITION = NULL;
  CONFIGFILEPARTITION = NULL;

  return rV;
}

int partitionParsewrap()
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: called when parser reaches EOF                         *
\***************************************************************/
{ // parser is complete
  return 1;
}

int partitionParseerror(char *s)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: called when parser reaches an error-should never occur *
\***************************************************************/
{ // error in parser
  printf("ERROR: There was an error (%s) near '%s'.\n", s, partitionParsetext);
  bexit(ERROR_INPUT_SYNTAX);

  return 0;
}




