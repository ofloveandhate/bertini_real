#include "output.h"
void Output_Main(char *inputName, char *input_deflated_Name,int *deflation_sequence,int num_vars, curveDecomp_d C, int MPType)
{
	int  strLength = 0;
	char *directoryName=NULL,tmp_file[1000];
	FILE *OUT;
	strLength = 1 + snprintf(NULL, 0, "%s_comp%d_curve", inputName,deflation_sequence[0]);
	directoryName = (char *)bmalloc(strLength * sizeof(char));
	sprintf(directoryName,  "%s_comp%d_curve", inputName,deflation_sequence[0]);
	
	mkdir(directoryName,0777);
	purge_previous_directory(directoryName);
	
	//Do we need witness_data?
	//sprintf(tmp_file,  "%s/witness_data", directoryName);
	//copyfile("witness_data",tmp_file);
	printf("before witness_set\n");fflush(stdout);
	
	sprintf(tmp_file,  "%s/witness_set", directoryName);
	copyfile("witness_set",tmp_file);
	
	sprintf(tmp_file,  "%s/%s", directoryName,input_deflated_Name);
	copyfile(input_deflated_Name,tmp_file);
	printf("after input_deflated_Name\n");fflush(stdout);
	sprintf(tmp_file,  "%s/Rand_Matrix", directoryName);
	copyfile("Rand_matrix",tmp_file);
	sprintf(tmp_file,  "%s/V0.vert", directoryName);
	print_vertices(C.V0,C.num_V0,num_vars,tmp_file, MPType);
	sprintf(tmp_file,  "%s/V1.vert", directoryName);
	print_vertices(C.V1,C.num_V1,num_vars,tmp_file, MPType);
	sprintf(tmp_file,  "%s/E.edge", directoryName);
	print_edges(C.E,C.num_E,num_vars,input_deflated_Name,tmp_file,MPType);
	
	OUT = fopen("Dir_Name", "w");
	fprintf(OUT,"%d\n",strLength);
	fprintf(OUT,"%s\n",directoryName);
	fprintf(OUT,"%d\n",MPType);
	
	fclose(OUT);
	free(directoryName);
}

void print_edges(edge_d *E, int num_E, int num_vars,char *input_deflated_Name,char *outputfile,int MPType)
/**Output edge structure as follows:
 # variables
 # edges
 name of input file
 edge 1
 
 edge 2
 .
 .
 .
 
 for each edge, output the following information:
 index to left vertex in V1
 index to right vertex in V1
 
 midpoint
 
 projection coefficients
 **/
{
	int i;
	FILE *OUT = fopen(outputfile, "w");
	
	if (OUT == NULL)
	{
		printf("\n\nERROR: '%s' is an improper name of a file!!\n\n\n", outputfile);
		bexit(ERROR_FILE_NOT_EXIST);
	}
	
	// output the number of vertices
	fprintf(OUT,"%d\n",num_vars);
	fprintf(OUT,"%d\n",num_E);
	fprintf(OUT,"%d\n",strlen(input_deflated_Name));
	fprintf(OUT,"%s\n",input_deflated_Name);
	for(i=0;i<num_E;i++)
		print_each_edge(E[i],num_vars,OUT,MPType);
	
	fclose(OUT);
	
}

void print_each_edge(edge_d E,int num_vars,FILE *OUT, int MPType)
/**
 for each edge, output the following information:
 index to left vertex in V1
 index to right vertex in V1
 
 midpoint
 
 projection coefficients
 
 **/
{
	int i;
	fprintf(OUT,"%d\n",E.left);
	fprintf(OUT,"%d\n\n",E.right);
	
	if(MPType==0)
	{
		for(i=0;i<num_vars;i++)
		{
			print_d(OUT, 0, &(E.midpt->coord[i]));
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
		for(i=0;i<num_vars;i++)
		{
			print_d(OUT, 0, &(E.pi->coord[i]));
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
	}
	else
	{
		for(i=0;i<num_vars;i++)
		{
			print_mp(OUT, 0, &(E.midpt_mp->coord[i]));
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
		for(i=0;i<num_vars;i++)
		{
			print_mp(OUT, 0, &(E.pi_mp->coord[i]));
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
		
	}
}

void print_vertices(vertex_d *V, int num_V, int num_vars,char *outputfile, int MPType)
/**Output vertex structure as follows:
 # pts
 pt.1
 
 pt.2
 
 .
 .
 .
 **/
{
	FILE *OUT = fopen(outputfile, "w");
	int i,j;
	
	if (OUT == NULL)
	{
		printf("\n\nERROR: '%s' is an improper name of a file!!\n\n\n", outputfile);
		bexit(ERROR_FILE_NOT_EXIST);
	}
	
	// output the number of vertices
	fprintf(OUT,"%d\n\n",num_V);
	
	for (i = 0; i < num_V; i++)
	{ // output points
		if(MPType==0)
		{
			for(j=0;j<num_vars;j++)
			{
				print_d(OUT, 0, &(V[i].pt->coord[j]));
				fprintf(OUT,"\n");
			}
		}
		else
		{
			for(j=0;j<num_vars;j++)
			{
				print_mp(OUT, 0, &(V[i].pt_mp->coord[j]));
				fprintf(OUT,"\n");
			}
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT);
	
}

void copyfile(char *INfile,char *OUTfile)
{
	char ch;
	FILE *IN,*OUT;
	IN = safe_fopen_read(INfile);
	OUT = fopen(OUTfile, "w");
	if (OUT == NULL)
	{
		printf("\n\nERROR: '%s' is an improper name of a file!!\n\n\n", OUTfile);
		bexit(ERROR_FILE_NOT_EXIST);
	}
	
	while ((ch = fgetc(IN)) != EOF)
		fprintf(OUT, "%c", ch);
	
	fclose(IN);
	fclose(OUT);
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


