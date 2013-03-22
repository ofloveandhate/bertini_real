#include "fileops.h"




FILE *safe_fopen_read(char * filename)
{
	
	FILE* IN;
	
	struct stat stat_p;
	stat (filename, &stat_p);
	
	if (S_ISDIR(stat_p.st_mode)){
		printf("trying to open directory %s as a file!\n",filename);
		bexit(ERROR_FILE_NOT_EXIST);
	}
	
	
	IN = fopen(filename,"r");
	
	
	if (IN == NULL) {
		printf("unable to open specified file to read: %s\n",filename);
		bexit(ERROR_FILE_NOT_EXIST);
	}
	else{
	}
	
	return IN;
	
}


FILE *safe_fopen_write(char * filename)
{
	
	FILE* OUT;
	
	
	OUT = fopen(filename,"w");
	
	
	if (OUT == NULL) {
		printf("unable to open specified file to write: %s\n",filename);
		bexit(ERROR_FILE_NOT_EXIST);
	}
	else{
	}
	
	return OUT;
	
}


FILE *safe_fopen_append(char * filename)
{
	
	FILE* OUT;
	
	
	
	struct stat stat_p;
	stat (filename, &stat_p);
	
	if (S_ISDIR(stat_p.st_mode)){
		printf("trying to open directory %s as a file!\n",filename);
		bexit(ERROR_FILE_NOT_EXIST);
	}
	
	
	OUT = fopen(filename,"a");
	
	
	if (OUT == NULL) {
		printf("unable to open specified file to write: %s\n",filename);
		bexit(ERROR_FILE_NOT_EXIST);
	}
	else{
	}
	
	return OUT;
	
}



