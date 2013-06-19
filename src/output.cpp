#include "output.hpp"




void Output_Main(program_configuration program_options, witness_set W, curveDecomp_d C, vertex_set V)
{
	int  strLength = 0;
	char *directoryName=NULL,tmp_file[1000];
	FILE *OUT;
	strLength = 1 + snprintf(NULL, 0, "%s_comp%d_curve", program_options.output_basename.c_str(), W.comp_num);
	directoryName = (char *)bmalloc(strLength * sizeof(char));
	sprintf(directoryName,  "%s_comp%d_curve", program_options.output_basename.c_str(), W.comp_num);
	
	mkdir(directoryName,0777);
	purge_previous_directory(directoryName);
	

	sprintf(tmp_file,  "%s/witness_data", directoryName);
	copyfile("witness_data",tmp_file);
//	printf("before witness_set\n");fflush(stdout);
	
	sprintf(tmp_file,  "%s/witness_set", directoryName);
	copyfile(program_options.witness_set_filename,tmp_file);
	
	sprintf(tmp_file,  "%s/%s", directoryName,program_options.input_deflated_filename.c_str());
	copyfile(program_options.input_deflated_filename,tmp_file);

	sprintf(tmp_file,  "%s/Rand_Matrix", directoryName);
	copyfile("Rand_matrix",tmp_file);
	
	sprintf(tmp_file,  "%s/V.vertex", directoryName);
	print_vertices(V, W.num_variables,tmp_file, program_options.MPType);
	sprintf(tmp_file,  "%s/E.edge", directoryName);
	print_edges(C.edges,C.num_edges,W.num_variables,
							tmp_file,program_options.MPType);
	
	sprintf(tmp_file,  "%s/C.curve", directoryName);
	print_curve(C, W.num_variables, program_options.input_deflated_filename.c_str(), tmp_file, program_options.MPType);
	
	OUT = safe_fopen_write("Dir_Name");
	fprintf(OUT,"%d\n",strLength);
	fprintf(OUT,"%s\n",directoryName);
	fprintf(OUT,"%d\n",program_options.MPType);
	fclose(OUT);
	
	
	free(directoryName);
}





void print_edges(edge *set_of_edges, int num_edges, int num_vars,char *outputfile,int MPType)
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
	fprintf(OUT,"%d\n\n",num_edges);

	for(ii=0;ii<num_edges;ii++)
		fprintf(OUT,"%d %d %d \n",set_of_edges[ii].left,
															set_of_edges[ii].midpt,
															set_of_edges[ii].right);
	
	

	fclose(OUT);
	
}

\

void print_vertices(vertex_set V, int num_vars,
										char *outputfile, int MPType)
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
	fprintf(OUT,"%d %d\n\n",V.num_vertices,num_vars);
	
	for (ii = 0; ii < V.num_vertices; ii++)
	{ // output points
			for(jj=0;jj<num_vars;jj++) {
				print_mp(OUT, 0, &(V.vertices[ii].pt_mp->coord[jj]));
				fprintf(OUT,"\n");
			}
			print_mp(OUT, 0, V.vertices[ii].projVal_mp);
			fprintf(OUT,"\n");
		
		
		
		fprintf(OUT,"%d\n",V.vertices[ii].type);
		
		
		fprintf(OUT,"\n"); // the line after the vertex is written
	}
	
	fclose(OUT);
	
}


/**
 Output curve overall info as follows:
 
 **/
void print_curve(curveDecomp_d C, int num_vars, boost::filesystem::path input_deflated_Name, boost::filesystem::path outputfile, int MPType)

{
	int ii;
	
	FILE *OUT = safe_fopen_write(outputfile.c_str());
	fprintf(OUT,"%zu\n",strlen(input_deflated_Name.c_str()));
	fprintf(OUT,"%s\n",input_deflated_Name.c_str());
	
	fprintf(OUT,"%d %d\n",num_vars, C.num_edges);
	fprintf(OUT,"%d %d %d %d %d\n",C.num_V0, C.num_V1, C.num_midpts, C.num_new, C.num_isolated);
	
	
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

	for(ii=0;ii<num_vars;ii++)
	{
		print_mp(OUT, 0, &(C.pi_mp->coord[ii]));
		fprintf(OUT,"\n");
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

