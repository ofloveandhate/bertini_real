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
	
//	sprintf(tmp_file,  "%s/V0.vert", directoryName);
//	print_vertices(C.V0,C.num_V0,num_vars,tmp_file, MPType);
	sprintf(tmp_file,  "%s/V.vert", directoryName);
	print_vertices(C.vertices,C.num_vertices,num_vars,tmp_file, MPType);
	sprintf(tmp_file,  "%s/E.edge", directoryName);
	print_edges(C.edges,C.num_edges,num_vars,input_deflated_Name,tmp_file,MPType);
	
	sprintf(tmp_file,  "%s/C.curve", directoryName);
	print_curve(C, num_vars, tmp_file, MPType);
	
	OUT = safe_fopen_write("Dir_Name");
	fprintf(OUT,"%d\n",strLength);
	fprintf(OUT,"%s\n",directoryName);
	fprintf(OUT,"%d\n",MPType);
	fclose(OUT);
	
	
	free(directoryName);
}





void print_edges(edge *set_of_edges, int num_edges, int num_vars,char *input_deflated_Name,char *outputfile,int MPType)
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
 index to left vertex in vertices
 index to right vertex in vertices
 index to midpoint vertex in vertices
 
 **/
{
	int ii;
	FILE *OUT = safe_fopen_write(outputfile);
	
	// output the number of vertices
	fprintf(OUT,"%d\n",num_vars);
	fprintf(OUT,"%d\n",num_edges);
	fprintf(OUT,"%zu\n",strlen(input_deflated_Name));
	fprintf(OUT,"%s\n",input_deflated_Name);
	for(ii=0;ii<num_edges;ii++)
		print_individual_edge(set_of_edges[ii],num_vars,OUT,MPType);
	
	fclose(OUT);
	
}

void print_individual_edge(edge curr_edge,int num_vars,FILE *OUT, int MPType)
/**
 for each edge, output the following information:
 index to left vertex in vertices
 index to right vertex in vertices
 index to midpoint in vertices
 **/
{
	fprintf(OUT,"%d\n%d\n%d\n\n",curr_edge.left,curr_edge.midpt,curr_edge.right);
}

void print_vertices(vertex *current_vertex, int num_V, int num_vars,char *outputfile, int MPType)
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
			print_d(OUT, 0, current_vertex[ii].projVal);
			fprintf(OUT,"\n");
			
		}
		else
		{
			for(jj=0;jj<num_vars;jj++)
			{
				print_mp(OUT, 0, &(current_vertex[ii].pt_mp->coord[jj]));
				fprintf(OUT,"\n");
			}
			print_mp(OUT, 0, current_vertex[ii].projVal_mp);
			fprintf(OUT,"\n");
		}
		
		fprintf(OUT,"%d\n",current_vertex[ii].type);
		
		
		fprintf(OUT,"\n"); // the line after the vertex is written
	}
	
	fclose(OUT);
	
}


/**
 Output curve overall info as follows:
 
 **/
void print_curve(curveDecomp_d C, int num_vars, char *outputfile, int MPType)

{
	int ii;
	
	FILE *OUT = safe_fopen_write(outputfile);
	fprintf(OUT,"%d %d %d\n",num_vars, C.num_vertices, C.num_edges);
	fprintf(OUT,"%d %d %d %d\n",C.num_V0, C.num_V1, C.num_midpts, C.num_new);
	

//	int			*V0_indices;
	for (ii=0; ii<C.num_V0; ii++) {
		fprintf(OUT,"%d\n",C.V0_indices[ii]);
	}
	fprintf(OUT,"\n");
	
	for (ii=0; ii<C.num_V1; ii++) {
		fprintf(OUT,"%d\n",C.V1_indices[ii]);
	}
	fprintf(OUT,"\n");
	
	for (ii=0; ii<C.num_midpts; ii++) {
		fprintf(OUT,"%d\n",C.midpt_indices[ii]);
	}
	fprintf(OUT,"\n");
	
	for (ii=0; ii<C.num_new; ii++) {
		fprintf(OUT,"%d\n",C.new_indices[ii]);
	}
	fprintf(OUT,"\n");

	// and finally the projection
	if (MPType==0) {
		for(ii=0;ii<num_vars;ii++)
		{
			print_d(OUT, 0, &(C.pi_d->coord[ii]));
			fprintf(OUT,"\n");
		}
	}
	else{
		for(ii=0;ii<num_vars;ii++)
		{
			print_mp(OUT, 0, &(C.pi->coord[ii]));
			fprintf(OUT,"\n");
		}
	}
	
	fclose(OUT);
	
}


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

