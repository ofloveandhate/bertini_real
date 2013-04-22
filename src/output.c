#include "output.h"




void Output_Main(char *outputName, char *input_deflated_Name,int component_number,int num_vars, curveDecomp_d C, int MPType)
{
	int  strLength = 0;
	char *directoryName=NULL,tmp_file[1000];
	FILE *OUT;
	strLength = 1 + snprintf(NULL, 0, "%s_comp%d_curve", outputName, component_number);
	directoryName = (char *)bmalloc(strLength * sizeof(char));
	sprintf(directoryName,  "%s_comp%d_curve", outputName, component_number);
	
	mkdir(directoryName,0777);
	purge_previous_directory(directoryName);
	
	//Do we need witness_data?
	
	sprintf(tmp_file,  "%s/witness_data", directoryName);
	copyfile("witness_data",tmp_file);
//	printf("before witness_set\n");fflush(stdout);
	
	sprintf(tmp_file,  "%s/witness_set", directoryName);
	copyfile("witness_set",tmp_file);
	
	sprintf(tmp_file,  "%s/%s", directoryName,input_deflated_Name);
	copyfile(input_deflated_Name,tmp_file);
//	printf("after input_deflated_Name\n");fflush(stdout);
	sprintf(tmp_file,  "%s/Rand_Matrix", directoryName);
	copyfile("Rand_matrix",tmp_file);
	
	sprintf(tmp_file,  "%s/V0.vert", directoryName);
	print_vertices(C.V0,C.num_V0,num_vars,tmp_file, MPType);
	sprintf(tmp_file,  "%s/V1.vert", directoryName);
	print_vertices(C.V1,C.num_V1,num_vars,tmp_file, MPType);
	sprintf(tmp_file,  "%s/E.edge", directoryName);
	print_edges(C.edges,C.num_edges,num_vars,input_deflated_Name,tmp_file,MPType);
	
	OUT = safe_fopen_write("Dir_Name");
	fprintf(OUT,"%d\n",strLength);
	fprintf(OUT,"%s\n",directoryName);
	fprintf(OUT,"%d\n",MPType);
	fclose(OUT);
	
	
	free(directoryName);
}





void print_edges(edge_d *set_of_edges, int num_edges, int num_vars,char *input_deflated_Name,char *outputfile,int MPType)
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
	int ii;
	FILE *OUT = safe_fopen_write(outputfile);
	
	// output the number of vertices
	fprintf(OUT,"%d\n",num_vars);
	fprintf(OUT,"%d\n",num_edges);
	fprintf(OUT,"%d\n",strlen(input_deflated_Name));
	fprintf(OUT,"%s\n",input_deflated_Name);
	for(ii=0;ii<num_edges;ii++)
		print_individual_edge(set_of_edges[ii],num_vars,OUT,MPType);
	
	fclose(OUT);
	
}

void print_individual_edge(edge_d curr_edge,int num_vars,FILE *OUT, int MPType)
/**
 for each edge, output the following information:
 index to left vertex in V1
 index to right vertex in V1
 
 midpoint
 
 projection coefficients
 
 **/
{
	int ii;
	fprintf(OUT,"%d\n",curr_edge.left);
	fprintf(OUT,"%d\n\n",curr_edge.right);
	
	if(MPType==0)
	{
		for(ii=0;ii<num_vars;ii++)
		{
			print_d(OUT, 0, &(curr_edge.midpt->coord[ii]));
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
		for(ii=0;ii<num_vars;ii++)
		{
			print_d(OUT, 0, &(curr_edge.pi->coord[ii]));
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
	}
	else
	{
		//print the point.
		for(ii=0;ii<num_vars;ii++)
		{
			print_mp(OUT, 0, &(curr_edge.midpt_mp->coord[ii]));
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
		
		//print the projection
		for(ii=0;ii<num_vars;ii++)
		{
			print_mp(OUT, 0, &(curr_edge.pi_mp->coord[ii]));
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
		
	}
}

void print_vertices(vertex_d *current_vertex, int num_V, int num_vars,char *outputfile, int MPType)
/**Output vertex structure as follows:
 # pts
 pt.1
 
 pt.2
 
 .
 .
 .
 **/
{
	
	FILE *OUT = safe_fopen_write(outputfile);
	int ii,jj;
	

	
	// output the number of vertices
	fprintf(OUT,"%d\n\n",num_V);
	
	for (ii = 0; ii < num_V; ii++)
	{ // output points
		if(MPType==0)
		{
			for(jj=0;jj<num_vars;jj++)
			{
				print_d(OUT, 0, &(current_vertex[ii].pt->coord[jj]));
				fprintf(OUT,"\n");
			}
		}
		else
		{
			for(jj=0;jj<num_vars;jj++)
			{
				print_mp(OUT, 0, &(current_vertex[ii].pt_mp->coord[jj]));
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
	
	IN  = safe_fopen_read(INfile);
	OUT = safe_fopen_write(OUTfile);

	while ((ch = fgetc(IN)) != EOF)
		fprintf(OUT, "%c", ch);
	
	fclose(IN);
	fclose(OUT);
}

