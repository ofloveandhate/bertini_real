#include "output.hpp"




void Output_Main(BR_configuration program_options, witness_set & W, curve_decomposition & C, vertex_set & V)
{

	FILE *OUT;
	std::stringstream converter;
	converter << "_comp" << W.comp_num << "_curve";
	boost::filesystem::path base = program_options.output_dir;
	base += converter.str(); converter.clear(); converter.str("");

	
//	purge_previous_directory(const_cast<char *>(base.c_str()));
	boost::filesystem::remove_all(base);
	boost::filesystem::create_directory(base);
	

	copyfile("witness_data",base / "witness_data");
	
	copyfile(program_options.witness_set_filename, base / "witness_set"); 
// TODO:  this should be a write call, not a copy
	copyfile(program_options.input_deflated_filename, base / program_options.input_deflated_filename.filename());

	copyfile("Rand_matrix", base / "Rand_Matrix");
	
	V.print(base / "V.vertex");
	
	C.print_edges(base / "E.edge");
	
	C.print(program_options.input_deflated_filename, base / "C.curve");
	
//	W.print_patches(base / "patches");
	OUT = safe_fopen_write("Dir_Name");
	fprintf(OUT,"%s\n",base.c_str());
	fprintf(OUT,"%d\n",program_options.MPType);
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

