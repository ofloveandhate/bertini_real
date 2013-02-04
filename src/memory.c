#include "memory.h"

void bexit(int errorCode)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: exits Bertini - either standard or using MPI           *
 \***************************************************************/
{
  if (errorCode == 0)
    errorCode = ERROR_OTHER;
	
  printf("%s\n", BERTINI_QUIT_MESSAGE);
#ifdef _HAVE_MPI
  MPI_Abort(MPI_COMM_WORLD, errorCode);
#else
  exit(errorCode);
#endif
}

void *bmalloc(size_t size)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: does malloc with error checking                        *
 \***************************************************************/
{
  if (size <= 0)
  { // nothing to allocate
    return NULL;
  }
  else
  { // try to allocate memory
    void *x = malloc(size);
    if (x == NULL)
    {
      printf("ERROR: malloc was unable to allocate memory (%d)!\n", (int) size);
      bexit(ERROR_MEMORY_ALLOCATION);
    }
    return x;
  }
}

void *brealloc(void *ptr, size_t size)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: does realloc with error checking                       *
 \***************************************************************/
{
  if (size <= 0)
  { // nothing to allocate - free memory and return NULL
    free(ptr);
    ptr = NULL;
  }
  else
  { // try to reallocate memory
    ptr = realloc(ptr, size);
    if (ptr == NULL)
    {
      printf("ERROR: realloc was unable to allocate memory!\n");
      bexit(ERROR_MEMORY_ALLOCATION);
    }
  }
  return ptr;
}
