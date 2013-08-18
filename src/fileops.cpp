#include "fileops.hpp"



int partition_parse(int **declarations,
										boost::filesystem::path input_filename,
										boost::filesystem::path functions_filename,
										boost::filesystem::path config_filename,
										int not_sc_flag)
{
	
	FILE *IN = safe_fopen_read(input_filename.c_str());
	
	int retval = partitionParse(declarations, IN,
															const_cast<char *> (functions_filename.c_str()),
															const_cast<char *> (config_filename.c_str()),
															not_sc_flag); // the 0 means not self conjugate.
	fclose(IN);
	
	return retval;
}



void rename_bertini_files_dotbak(){
	
	rename("arr.out","arr.out.bak");
	rename("num.out","num.out.bak");
	rename("deg.out","deg.out.bak");
	rename("config","config.bak");
	rename("preproc_data","preproc_data.bak");
}

void restore_bertini_files_dotbak(){
	rename("config.bak","config");
	rename("arr.out.bak","arr.out");
	rename("num.out.bak","num.out");
	rename("deg.out.bak","deg.out");
	rename("preproc_data.bak","preproc_data");
	
}

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



FILE *safe_fopen_read(boost::filesystem::path filename)
{
	
	FILE* IN;
	
	if (boost::filesystem::is_directory(filename)){
		std::cerr << "trying to open directory " << filename.string() << " as a file!" << std::endl;
		deliberate_segfault();
		br_exit(ERROR_FILE_NOT_EXIST);
	}
	
	if (!boost::filesystem::exists(filename)){
		std::cerr << "unable to open specified file to read: " << filename.string() << std::endl;
		std::cerr << boost::filesystem::current_path().string() << " is the current directory" << std::endl;
		deliberate_segfault();
		br_exit(ERROR_FILE_NOT_EXIST);
	}
	
	
	IN = fopen(filename.c_str(),"r");

	if (IN == NULL) {
		std::cerr << "failed to open file: " << filename.string() << std::endl;
		fclose(IN);
		deliberate_segfault();
		br_exit(ERROR_FILE_NOT_EXIST);
	}
	else{
	}
	
	return IN;
	
}


FILE *safe_fopen_write(boost::filesystem::path filename)
{
	
	FILE* OUT;
	
	
	OUT = fopen(filename.c_str(),"w");
	
	
	if (OUT == NULL) {
		printf("unable to open specified file to write: %s\n",filename.c_str());
		fclose(OUT);
		deliberate_segfault();
//		br_exit(ERROR_FILE_NOT_EXIST);
	}
	else{
	}
	
	return OUT;
	
}


FILE *safe_fopen_append(boost::filesystem::path filename)
{
	FILE* OUT;

	struct stat stat_p;
	stat (filename.c_str(), &stat_p);
	
	if (S_ISDIR(stat_p.st_mode)){
		printf("trying to open directory %s as a file!\n",filename.c_str());
		deliberate_segfault();
//		br_exit(ERROR_FILE_NOT_EXIST);
	}
	
	
	OUT = fopen(filename.c_str(),"a");
	
	
	if (OUT == NULL) {
		printf("unable to open specified file to write: %s\n",filename.c_str());
		fclose(OUT);
		br_exit(ERROR_FILE_NOT_EXIST);
	}
	else{
	}
	
	return OUT;
}


void copyfile(boost::filesystem::path input_file, boost::filesystem::path OUTfile)
{
	char ch;
	FILE *IN,*OUT;
	
	IN  = safe_fopen_read(const_cast< char*> (input_file.c_str()));
	OUT = safe_fopen_write(const_cast< char*> (OUTfile.c_str()));
	
	while ((ch = fgetc(IN)) != EOF)
		fprintf(OUT, "%c", ch);
	
	fclose(IN);
	fclose(OUT);
}



void copyfile(FILE *IN,FILE *OUT)
{
	char ch;
	while ((ch = fgetc(IN)) != EOF)
		fprintf(OUT, "%c", ch);

}




void read_matrix(boost::filesystem::path INfile, mat_mp matrix)
{
	int rows,cols, precision;
	FILE *IN= safe_fopen_read(INfile);
	fscanf(IN,"%d %d", &rows, &cols); scanRestOfLine(IN);
	fscanf(IN,"%d", &precision); scanRestOfLine(IN);
	init_mat_mp2(matrix,rows,cols,precision);
	matrix->rows=rows; matrix->cols=cols;
	
	for (int ii=0;ii<matrix->rows;ii++)
		for (int jj=0;jj<matrix->cols;jj++)
		{
			mpf_inp_str(matrix->entry[ii][jj].r, IN, 10);
			mpf_inp_str(matrix->entry[ii][jj].i, IN, 10);
			scanRestOfLine(IN);
		}
	
	fclose(IN);
}



