#include "output.hpp"













/********************************************************/
void print_matrix_to_file_mp(FILE *OUT, int digits, mat_mp M)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES:                                                        *
 \***************************************************************/
{
  int i, j, rows = M->rows, cols = M->cols;
	
  for (i = 0; i < rows; i++)
  {
    for (j = 0; j < cols; j++)
    {
      print_mp(OUT, digits, &M->entry[i][j]);
      fprintf(OUT, "\n");
    }
    fprintf(OUT, "\n");
  }
  fprintf(OUT, "\n");
	
  return;
}

