#include "fileops.h"

//will remove all files which do not start with a . from directory.
void purge_previous_directory(char *directoryName)
{
	DIR *dp;
	struct dirent *ep;
	char tempname[1000];
	
	
	dp = opendir (directoryName);
	if (dp != NULL)
	{
		while ( (ep = readdir (dp)) )
			if (ep->d_name[0] != '.')
	    {
				sprintf(tempname,"%s/%s",directoryName, ep->d_name);
				remove(tempname);
			}
		
		
		(void) closedir (dp);
	}
	return;
}



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



