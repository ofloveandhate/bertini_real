#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#ifndef _MEMORY_H
#define _MEMORY_H

// EXIT
#define BERTINI_QUIT_MESSAGE "BertiniReal will now exit due to this error.\n"
#define ERROR_OTHER 1
#define ERROR_WRITE_PRIVILEGE 2   // write privilege problems - e.g. not able to create and write to files
#define ERROR_FILE_NOT_EXIST 3    // file exist problems - e.g. expected files do not exist
#define ERROR_INVALID_SIZE 4      // size problems - e.g. expected size not the same as the current size
#define ERROR_MEMORY_ALLOCATION 5 // memory problems - e.g. unable to allocate memory
#define ERROR_CONFIGURATION 6     // configuration problems - e.g. function calls called with wrong input values
#define ERROR_INPUT_SYSTEM 7      // input errors - e.g. degenerate system
#define ERROR_INPUT_SYNTAX 8      // syntax errors for input system
#define ERROR_LOOP_FAILURE 9      // loop failed to exit properly - e.g. fail-safe for 'infinite loops'



void bexit(int errorCode); // Bertini's exit function


void *bmalloc(size_t size);
void *brealloc(void *ptr, size_t size);

#endif
