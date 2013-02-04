%{
/* This will take an input file and partition it into a configurations file and a function file */
#include "bertini_real.h"

// settings to switch between config & function
int configFound = 0, endOfConfig = 1, endOfFunc = 0;
// number of times the statement appears
int num_variable_group = 0, num_hom_variable_group = 0, num_variable = 0, num_pathvariable = 0, num_parameter = 0, num_definedSubfunction = 0, num_random = 0, num_random_real = 0, num_constant = 0, num_function = 0;
// files for config & function
FILE *CONFIGFILE = NULL, *FUNCFILE = NULL;

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
             if (!endOfConfig)
             { // end config
               endOfConfig = 1;
             }
             else if (!endOfFunc)
             { // end func
               endOfFunc = 1;
             }
           }
{comments}|{white} { // ignore comments and whitespace
           }
{config}   { // start of config found - print there until reach END or INPUT
             configFound = 1;
             endOfConfig = 0;
           }
{input}    { // start of input found - print there until reach END or EOF
             endOfConfig = 1;
             endOfFunc = 0;
           }
{param}    { // declare parameter
             if (endOfConfig && !endOfFunc)
             { // print to FUNCFILE
               num_parameter++;
               fprintf(FUNCFILE, "%s ", partitionParsetext);
             }
           }
{homvargp} { // declare hom_variable_group
             if (endOfConfig && !endOfFunc)
             { // print to FUNCFILE
               num_hom_variable_group++;
               fprintf(FUNCFILE, "%s ", partitionParsetext);
             }
           }
{vargp}    { // declare variable_group
             if (endOfConfig && !endOfFunc)
             { // print to FUNCFILE
               num_variable_group++;
               fprintf(FUNCFILE, "%s ", partitionParsetext);
             }
           }
{var}      { // declare variable
             if (endOfConfig && !endOfFunc)
             { // print to FUNCFILE
               num_variable++;
               fprintf(FUNCFILE, "%s ", partitionParsetext);
             }
           }
{pathvar}  { // declare pathvariable
             if (endOfConfig && !endOfFunc)
             { // print to FUNCFILE
               num_pathvariable++;
               fprintf(FUNCFILE, "%s ", partitionParsetext);
             }
           }
{const}    { // declare constant
             if (endOfConfig && !endOfFunc)
             { // print to FUNCFILE
               num_constant++;
               fprintf(FUNCFILE, "%s ", partitionParsetext);
             }
           }
{func}     { // declare function
             if (endOfConfig && !endOfFunc)
             { // print to FUNCFILE
               num_function++;
               fprintf(FUNCFILE, "%s ", partitionParsetext);
             }
           }
{subfunc}  { // declare definedSubfunction
             if (endOfConfig && !endOfFunc)
             { // print to FUNCFILE
               num_definedSubfunction++;
               fprintf(FUNCFILE, "%s ", partitionParsetext);
             }
           }
{rand}     { // declare random
             if (endOfConfig && !endOfFunc)
             { // print to FUNCFILE
               num_random++;
               fprintf(FUNCFILE, "%s ", partitionParsetext);
             }
           }
{randReal} { // declare random_real
             if (endOfConfig && !endOfFunc)
             { // print to FUNCFILE
               num_random_real++;
               fprintf(FUNCFILE, "%s ", partitionParsetext);
             }
           }
{colon}    { // a colon
             if (!endOfConfig)
             { // add a space after a colon in configurations
               fprintf(CONFIGFILE, ": ");
             }
           }
{eol}      { // end of line
             if (!endOfConfig)
             { // print to CONFIGFILE
  	       fprintf(CONFIGFILE, ";\n");
             }
             else if (!endOfFunc)
             { // print to FUNCFILE
               fprintf(FUNCFILE, ";\n");
             }
           }
{number}|{I}|{name}|{other}   { // other symbol
             if (!endOfConfig)
             { // print to CONFIGFILE
               fprintf(CONFIGFILE, "%s", partitionParsetext);
             }
             else if (!endOfFunc)
             { // print to FUNCFILE
               fprintf(FUNCFILE, "%s", partitionParsetext);
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

  // setup FUNCFILE & CONFIGFILE
  FUNCFILE = fopen(funcName, "w");
  CONFIGFILE = fopen(configName, "w");
  if (FUNCFILE == NULL)
  {
    printf("ERROR: The name %s is not valid!\n\n", funcName);
    bexit(ERROR_INPUT_SYNTAX);
  }
  if (CONFIGFILE == NULL)
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

  // close FUNCFILE & CONFIGFILE
  if (!endOfFunc)
    fprintf(FUNCFILE, "\nEND;\n\n");

  fclose(FUNCFILE);
  fclose(CONFIGFILE);
  FUNCFILE = NULL;
  CONFIGFILE = NULL;

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




