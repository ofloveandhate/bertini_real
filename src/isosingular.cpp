#include "isosingular.hpp"



int isosingular_deflation(int *num_deflations, int **deflation_sequence, char *inputFile, char *witness_point_filename, char *bertini_command, char *matlab_command, int max_deflations,
													int dim, int component_number)
/***************************************************************\
* USAGE: Perform isosingular deflation for system at given point*
* ARGUMENTS: name of file for polynomial system, point, command *
*  to call Bertini & Matlab, and maximum number of deflations   *
* RETURN VALUES: success or failure, number of deflations, and  *
*  deflation sequence                                           *
* NOTES:                                                        *
\***************************************************************/
{
  int ii, success = 0, strLength = 0, nullSpace = 0, *declarations = NULL;
  char ch, *strStabilizationTest = NULL;
  FILE *IN = NULL, *OUT = NULL;

	// create the string for running system()
  strLength = 1 + snprintf(NULL, 0, "%s input_stabilization_test %s", bertini_command, witness_point_filename);
  strStabilizationTest = (char *)bmalloc(strLength * sizeof(char));
  sprintf(strStabilizationTest, "%s input_stabilization_test %s", bertini_command, witness_point_filename);
	
	
	// remove previous files.

	
	//open the input file.
	IN = safe_fopen_read(inputFile);	
  partitionParse(&declarations, IN, "func_input_real", "config_real",0); // the 0 means not self conjugate mode
//	check_declarations(declarations);
  fclose(IN);


  // setup input file to test for stabilization
  stabilization_input_file("input_stabilization_test", "func_input_real", "config_real");

	
  printf("\nPerforming a stabilization test\n");

	

	remove("isosingular_summary");
		
	//print the command to the screen
	printf("running system command '%s'\n",strStabilizationTest);
	
	// perform stabilization test
	if (system(strStabilizationTest)!=0){
		bexit(ERROR_CONFIGURATION);
	}


	
  // read in the file	
	
	IN = safe_fopen_read("isosingular_summary");
  fscanf(IN, "%d%d", &nullSpace, &success);     // set the success flag
	fclose(IN);
  // setup the first entry in the deflation sequence
  *num_deflations = 0;
  *deflation_sequence = (int *)bmalloc((*num_deflations + 1) * sizeof(int));
  (*deflation_sequence)[*num_deflations] = nullSpace;
 
  // loop until successful or run out of iterations
  while (*num_deflations < max_deflations && !success)
  { // create input file for deflation
    isosingular_deflation_iteration(declarations, "func_input_real", matlab_command, nullSpace, *num_deflations + 1);

    // setup input file to test for stabilization
    stabilization_input_file("input_stabilization_test", "func_input_real", "config_real");


    // perform stabilization test
    printf("\nPerforming a stabilization test\n");
    if (system(strStabilizationTest)!=0){
			bexit(ERROR_CONFIGURATION);
		}
		
		
    // read in the file
		IN = safe_fopen_read("isosingular_summary");
    fscanf(IN, "%d%d", &nullSpace, &success); // get the success indicator

    // setup the next entry in the deflation sequence
    (*num_deflations)++;
    *deflation_sequence = (int *)brealloc(*deflation_sequence, (*num_deflations + 1) * sizeof(int));
    (*deflation_sequence)[*num_deflations] = nullSpace;

    if ((*deflation_sequence)[*num_deflations] > (*deflation_sequence)[*num_deflations - 1])
    {
      printf("ERROR: The deflation sequence must be a nonincreasing sequence of nonnegative integers!\n\n");
      success = 0; 
      break;
    }
  }

  if (success)
  { // print deflation sequence
    printf("\nIsosingular deflation was successful!\n\n");
    printf("Number of deflations: %d\n", *num_deflations);
    printf("Deflation sequence: ");
    for (ii = 0; ii < *num_deflations; ii++)
      printf("%d, ", (*deflation_sequence)[ii]);
    printf("%d, ...\n\n", (*deflation_sequence)[*num_deflations]);

    // create deflated system
    strLength = 1 + snprintf(NULL, 0, "%s_dim_%d_comp_%d_deflated", inputFile, dim, component_number);
    strStabilizationTest = (char *)bmalloc(strLength * sizeof(char));
    sprintf(strStabilizationTest, "%s_dim_%d_comp_%d_deflated", inputFile, dim, component_number);
    OUT = fopen(strStabilizationTest, "w");
    fprintf(OUT, "CONFIG\n");
		IN = safe_fopen_read("config_real");
    while ((ch = fgetc(IN)) != EOF)
      fprintf(OUT, "%c", ch);
    fclose(IN);
    fprintf(OUT, "END;\nINPUT\n");
		IN = safe_fopen_read("func_input_real");
    while ((ch = fgetc(IN)) != EOF)
      fprintf(OUT, "%c", ch);
    fclose(IN);
    fprintf(OUT, "END;\n\n");
    fclose(OUT);

    printf("Deflated system printed to '%s'.\n\n", strStabilizationTest);

  }
  else if (*num_deflations >= max_deflations)
  {
    printf("ERROR: The maximum number of deflations (%d) was exceeded without stabilization!\n\n", max_deflations);
  }
  // delete temporary files
  remove("func_input_real");
  remove("config_real");

  // clear memory
  free(strStabilizationTest);
  free(declarations);

  return success;
}










void isosingular_deflation_iteration(int *declarations, char *inputOutputName, char *matlab_command, int nullSpaceDim, int deflation_number)
/***************************************************************\
* USAGE: setup input file for one deflation iteration           *
* ARGUMENTS: number of declaration statments, name of file,     *
*  command to run Matlab, and dimension of nullspace            *
* RETURN VALUES: creates new file                               *
* NOTES:                                                        *
\***************************************************************/
{
  int ii, numVars = 0, numFuncs = 0, numConstants = 0, minorSize = 0, strLength = 0;
  int *lineVars = NULL, *lineFuncs = NULL, *lineConstants = NULL, *degrees = NULL;
  char ch, *str = NULL, **vars = NULL, **funcs = NULL, **consts = NULL;
  FILE *IN = NULL, *OUT = NULL;

  // test for existence of input file
	IN = safe_fopen_read(inputOutputName);
  fclose(IN);
   
  // move the file & open it
  rename(inputOutputName, "deflation_input_file");
	IN = safe_fopen_read("deflation_input_file");
  // setup variables
  if (declarations[0] > 0)
  { // using variable_group
    parse_names(&numVars, &vars, &lineVars, IN, "variable_group", declarations[0]);
  }
  else
  { // using hom_variable_group
    parse_names(&numVars, &vars, &lineVars, IN, "hom_variable_group", declarations[1]);
  }

  // setup constants 
  rewind(IN);
  parse_names(&numConstants, &consts, &lineConstants, IN, "constant", declarations[8]);

  // setup functions
  rewind(IN);
  parse_names(&numFuncs, &funcs, &lineFuncs, IN, "function", declarations[9]);

  // read in the degrees
  degrees = (int *)bmalloc(numFuncs * sizeof(int));
	OUT = safe_fopen_read("deg.out");
  for (ii = 0; ii < numFuncs; ii++)
    fscanf(OUT, "%d", &degrees[ii]);
  fclose(OUT);;

  // setup Matlab script
  rewind(IN);
//  OUT = fopen("matlab_deflate.m", "w");
	OUT = safe_fopen_write("matlab_deflate.m");
  minorSize = numVars - declarations[1] - nullSpaceDim + 1;
  createMatlabDeflation(OUT, numVars, vars, lineVars, numConstants, consts, lineConstants, numFuncs, funcs, lineFuncs, IN, minorSize, degrees, deflation_number);
  fclose(OUT);  

  // run Matlab script
  printf("\nPerforming an isosingular deflation\n");
  strLength = 1 + snprintf(NULL, 0, "%s < matlab_deflate.m", matlab_command);
  str = (char *)bmalloc(strLength * sizeof(char));
  sprintf(str, "%s < matlab_deflate.m", matlab_command);
  system(str);
  printf("\n");

	
  // setup new file
  OUT = safe_fopen_write(inputOutputName);
  rewind(IN);
  while ((ch = fgetc(IN)) != EOF) 
    fprintf(OUT, "%c", ch);
  fclose(IN);

  
	IN = safe_fopen_read("deflation_polynomials_declaration");
  while ((ch = fgetc(IN)) != EOF) 
    fprintf(OUT, "%c", ch);
  fclose(IN);
	
	
	
	IN = safe_fopen_read("deflation_polynomials");
  while ((ch = fgetc(IN)) != EOF)
    fprintf(OUT, "%c", ch);
  fclose(IN);
  fclose(OUT);

  // increment the number of function declarations
  declarations[9]++;

  // clear memory
  free(str);
  free(degrees);
  for (ii = 0; ii < numVars; ii++)
    free(vars[ii]);
  free(vars);
  free(lineVars);
  for (ii = 0; ii < numConstants; ii++)
    free(consts[ii]);
  free(consts);
  free(lineConstants);
  for (ii = 0; ii < numFuncs; ii++)
    free(funcs[ii]);
  free(funcs);
  free(lineFuncs);

  return;
}

void createMatlabDeflation(FILE *OUT, int numVars, char **vars, int *lineVars, int numConstants, char **consts, int *lineConstants, int numFuncs, char **funcs, int *lineFuncs, FILE *IN, int minorSize, int *degrees, int deflation_number)
/***************************************************************\
* USAGE: setup a Matlab script to perform the deflation         *
* ARGUMENTS: input file, declaration name, and number of lines  *
* RETURN VALUES: Matlab script                                  *
* NOTES:                                                        *
\***************************************************************/
{
  int ii, lineNumber = 1, cont = 1, declares = 0, strLength = 0, strSize = 1;
  char ch;
  char *str = (char *)bmalloc(strSize * sizeof(char));

  // setup Bertini constants in Matlab
  fprintf(OUT, "syms I Pi;\n");

  // setup variables
  fprintf(OUT, "syms");
  for (ii = 0; ii < numVars; ii++)
    fprintf(OUT, " %s", vars[ii]);
  fprintf(OUT, ";\nX = [");
  for (ii = 0; ii < numVars; ii++)
    fprintf(OUT, " %s", vars[ii]);
  fprintf(OUT, "];\n");

  // setup constants
  if (numConstants > 0)
  {
    fprintf(OUT, "syms");
    for (ii = 0; ii < numConstants; ii++)
      fprintf(OUT, " %s", consts[ii]);
    fprintf(OUT, ";\n");
  }

  // setup degrees
  fprintf(OUT, "deg = [");
  for (ii = 0; ii < numFuncs; ii++)
    fprintf(OUT, " %d", degrees[ii]);
  fprintf(OUT, "];\n");

  // copy lines which do not declare items or define constants (keep these as symbolic objects)
  while (cont)
  { // see if this line number declares items
    declares = 0;
    for (ii = 0; ii < numVars; ii++)
      if (lineNumber == lineVars[ii])
        declares = 1;
    for (ii = 0; ii < numConstants; ii++)
      if (lineNumber == lineConstants[ii])
        declares = 1;
    for (ii = 0; ii < numFuncs; ii++)
      if (lineNumber == lineFuncs[ii])
        declares = 1;

    if (declares)
    { // move past this line
      do
      { // read in character
        ch = fgetc(IN);
      } while (ch != '\n' && ch != EOF);
    }
    else
    { // check to see if this defines a constant - line must be of the form '[NAME]=[EXPRESSION];' OR EOF
      ch = fgetc(IN);
      if (ch != EOF)
      { // definition line
        strLength = 0;
        do
        { // add to str
          if (strLength + 1 == strSize)
          { // increase strSize
            strSize *= 2;
            str = (char *)brealloc(str, strSize * sizeof(char));
          }
          str[strLength] = ch;
          strLength++;
        } while ((ch = fgetc(IN)) != '=');
        str[strLength] = '\0';
        strLength++;

        // compare against constants
        declares = 0;
        for (ii = 0; ii < numConstants; ii++)
          if (strcmp(str, consts[ii]) == 0)
            declares = 1;

        if (declares)
        { // move past this line
          do
          { // read in character
            ch = fgetc(IN);
          } while (ch != '\n' && ch != EOF);
        } 
        else
        { // print line
          fprintf(OUT, "%s=", str);
          do
          { // read in character & print it
            ch = fgetc(IN);
            if (ch != EOF)
              fprintf(OUT, "%c", ch);
          } while (ch != '\n' && ch != EOF);
        }
      }
    }

    // increment lineNumber
    lineNumber++;

    // test for EOF
    if (ch == EOF)
      cont = 0;
  }  

  // setup functions
  fprintf(OUT, "\nF = [");
  for (ii = 0; ii < numFuncs; ii++)
    fprintf(OUT, " %s", funcs[ii]);
  fprintf(OUT, "];\n");

  // compute the jacobian
  fprintf(OUT, "J = jacobian(F,X);\n");

  // find the combinations of the rows & columns
  fprintf(OUT, "R = nchoosek(1:%d,%d);\nC = nchoosek(1:%d,%d);\n", numFuncs, minorSize, numVars, minorSize); // changed to nchoosek april 16, 2013 DAB.
  fprintf(OUT, "r = size(R,1);\nc = size(C,1);\n");

  // loop over rows & columns printing the nonzero minors to a file
  fprintf(OUT, "count = 0;\n");
  fprintf(OUT, "OUT = fopen('deflation_polynomials','w');\n");
  fprintf(OUT, "for j = 1:r\n");
  fprintf(OUT, "  for k = 1:c\n");
  fprintf(OUT, "    A = simplify(det(J(R(j,:),C(k,:))));\n");
  fprintf(OUT, "    if A ~= 0\n");
  fprintf(OUT, "      count = count + 1;\n");
  fprintf(OUT, "      fprintf(OUT, 'f_%d_%cd = %cs;%cn', count, char(A/prod(deg(R(j,:)))));\n", deflation_number, '%', '%', '\\');
  fprintf(OUT, "    end;\n");
  fprintf(OUT, "  end;\n");
  fprintf(OUT, "end;\n");
  fprintf(OUT, "fclose(OUT);\n");

  // print number of minors to a file
  fprintf(OUT, "OUT = fopen('deflation_polynomials_declaration','w');\n");
  fprintf(OUT, "fprintf(OUT, 'function ');\n");
  fprintf(OUT, "for j = 1:count\n");
  fprintf(OUT, "  fprintf(OUT, 'f_%d_%cd', j);\n", deflation_number, '%');
  fprintf(OUT, "  if j == count");
  fprintf(OUT, "    fprintf(OUT, ';%cn');\n", '\\');
  fprintf(OUT, "  else");
  fprintf(OUT, "    fprintf(OUT, ',');\n");
  fprintf(OUT, "  end;\n");
  fprintf(OUT, "end;\n");

  // clear memory
  free(str);

  return;
}

void parse_names(int *numItems, char ***itemNames, int **itemLines, FILE *IN, char *name, int num_declarations)
/***************************************************************\
* USAGE: find the lines where items are declared and find names *
* ARGUMENTS: input file, declaration name, and number of lines  *
* RETURN VALUES: number of items, names, and line numbers       *
* NOTES:                                                        *
\***************************************************************/
{
  int lineNumber = 1, lengthName = strlen(name), strLength = 0, strSize = 1;
  char ch;
  char *str = (char *)bmalloc(strSize * sizeof(char));

  // move through the file looking for the items
  while ((ch = fgetc(IN)) != EOF)
  { // add to str
    if (strLength + 1 == strSize)
    { // increase strSize
      strSize *= 2;
      str = (char *)brealloc(str, strSize * sizeof(char));
    }
    str[strLength] = ch;
    strLength++;

    if (strLength == lengthName)
    { // compare against name
      str[strLength] = '\0';
      if (strcmp(str, name) == 0)
      { // we have a good line -- add the items
        ch = fgetc(IN); // skip the space
        addItems(numItems, itemNames, itemLines, IN, lineNumber);
      }
    }

    if (ch == '\n')
    { // new line - reset data
      lineNumber++;
      strLength = 0;
    }
  }

  // clear memory
  free(str);

  return;
}

void addItems(int *numItems, char ***itemNames, int **itemLines, FILE *IN, int lineNumber)
/***************************************************************\
* USAGE: add the items in the current line                      *
* ARGUMENTS: input file and current line number                 *
* RETURN VALUES: number of items, names, and line numbers       *
* NOTES:                                                        *
\***************************************************************/
{
  int strLength = 0, strSize = 1, cont = 1;
  char ch;
  char *str = (char *)bmalloc(strSize * sizeof(char));

  // initialize ch
  ch = fgetc(IN);
  do
  { // read in next character
    if (ch == ';')
    { // end loop
      cont = 0;
    }
    else
    { // we have a new item to add
      (*numItems)++;
      *itemNames = (char **)brealloc(*itemNames, *numItems * sizeof(char *));
      *itemLines = (int *)brealloc(*itemLines, *numItems * sizeof(int));

      // read in the name
      do
      { // save the character 
        if (strLength + 1 == strSize)
        { // increase strSize
          strSize *= 2;
          str = (char *)brealloc(str, strSize * sizeof(char));
        }
        str[strLength] = ch;
        strLength++;

        // read in the next character
        ch = fgetc(IN);

      } while (ch != ',' && ch != ';');

      // save the information
      str[strLength] = '\0';
      strLength++;
      (*itemNames)[*numItems - 1] = (char *)bmalloc(strLength * sizeof(char));
      strcpy((*itemNames)[*numItems - 1], str);
      (*itemLines)[*numItems - 1] = lineNumber;

      // reset strLength
      strLength = 0;

      if (ch == ',')
      { // read in the next character
        ch = fgetc(IN);
      }
    }
  } while (cont);  

  free(str);

  return;
}

void stabilization_input_file(char *outputFile, char *funcInput, char *configInput)
/***************************************************************\
* USAGE: setup input file to test for stabilization             *
* ARGUMENTS: name of output file, function & configuration input*
* RETURN VALUES: none                                           *
* NOTES:                                                        *
\***************************************************************/
{
  char ch;
  FILE *OUT = fopen(outputFile, "w"), *IN = NULL;
  if (OUT == NULL)
  {
    printf("\n\nERROR: '%s' is an improper name of a file!!\n\n\n", outputFile);
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
  fprintf(OUT, "TrackType: 6;\nReducedOnly: 1;\nDeleteTempFiles: 0;\nEND;\nINPUT\n");

  // setup system in OUT
//  IN = fopen(funcInput, "r");
	IN = safe_fopen_read(funcInput);
//  if (IN == NULL)
//  {
//    fclose(OUT);
//    printf("\n\nERROR: '%s' does not exist!!!\n\n\n", funcInput);
//    bexit(ERROR_FILE_NOT_EXIST);
//  }
  while ((ch = fgetc(IN)) != EOF)
    fprintf(OUT, "%c", ch);
  fclose(IN);
  fprintf(OUT, "END;\n");
  fclose(OUT);

  return;
}





void check_declarations(int *declarations){
	
  // error check based on declarations
  if ( (declarations[0] + declarations[1]) > 1)
  { // multiple variable groups
    printf("\n\n%d variable_groups, %d hom_variable_groups\nERROR: Please use either one 'variable_group' or one 'hom_variable_group'!\n\n",declarations[0],declarations[1]);
    bexit(ERROR_CONFIGURATION);
  }
  else if (declarations[2] > 0)
  { // 'variable'
    printf("\n\nERROR: The system needs to be defined using either 'variable_group' or 'hom_variable_group'.\n\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (declarations[3] > 0)
  { // 'pathvariable'
    printf("\n\nERROR: The system should not depend on pathvariables!\n\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (declarations[4] > 0)
  { // 'parameter'
    printf("\n\nERROR: The system should not depend on parameters!\n\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (declarations[5] > 0)
  { // 'definedSubfunction'
    printf("\n\nERROR: The system should not depend on defined subfunctions!\n\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (declarations[6] > 0)
  { // 'random'
    printf("\n\nERROR: The system should not involve 'random' numbers!\n\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (declarations[7] > 0)
  { // 'random_real'
    printf("\n\nERROR: The system should not involve 'random_real' numbers!\n\n");
    bexit(ERROR_CONFIGURATION);
  }
	
}

