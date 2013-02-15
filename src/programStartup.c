#include "programStartup.h"



int startup(int argC, char *args[], char **inputName, char **startName)
/***************************************************************\
 * USAGE:    prepares the variables inputname and startname 
 *      for use later in the program
 * ARGUMENTS:                                                    *
 * RETURN VALUES:           an integer                           *
 * NOTES:                                                        *
 \***************************************************************/
{

	printf("\n BertiniReal(TM) v%s\n\n", BERTINI_REAL_VERSION_STRING);
  printf(" D.J. Bates, D. Brake,\n W. Hao, J.D. Hauenstein,\n A.J. Sommese, C.W. Wampler\n\n");
  printf("(using GMP v%d.%d.%d, MPFR v%s)\n\n", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());
	
  // check for write privilege
  if (checkWritePrivilege())
  {
    printf("ERROR: BertiniReal does not have write privileges!\n");
    bexit(ERROR_WRITE_PRIVILEGE);
  }
	
  if (argC > 1 && args[1] != NULL && (!strcmp(args[1], "--help") || !strcmp(args[1], "-help"))) // help
  { // print information about Bertini
    printf("\nThis is BertiniReal v%s, developed by\nDan J. Bates, Daniel Brake,\nWenrui Hao, Jonathan D. Hauenstein,\nAndrew J. Sommmese, and Charles W. Wampler.\n\n", BERTINI_REAL_VERSION_STRING);
    printf("See ??? for details about BertiniReal.\n\n");
  }
  else if (argC > 1 && args[1] != NULL && (!strcmp(args[1], "--version") || !strcmp(args[1], "-version"))) // version
  { // simply exit
  }
  else
  { 
		
    // setup inputName & startName
    if (argC >= 2 && args[1] != NULL)
    { // inputName is args[1]
      *inputName = (char *)bmalloc((strlen(args[1]) + 1) * sizeof(char));
      strcpy(*inputName, args[1]);
			
      // setup startName
      if (argC >= 3 && args[2] != NULL)
      { // startName is args[2]
        *startName = (char *)bmalloc((strlen(args[2]) + 1) * sizeof(char));
        strcpy(*startName, args[2]);
      }
      else
      { // default to 'start'
        *startName = (char *)bmalloc(12 * sizeof(char));
        strcpy(*startName, "witness_set");
      }
    }
    else
    { // default to 'input' & 'start'
      *inputName = (char *)bmalloc(6 * sizeof(char));
      strcpy(*inputName, "input");
      *startName = (char *)bmalloc(12 * sizeof(char));
      strcpy(*startName, "witness_set");
    }

	}
		return 0;
}

