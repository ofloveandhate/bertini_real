// brings in necessary modules, functions, and parameters
// this file is used in `br_surf_complete.scad`
//
// depends on `plugs_and_sockets.scad`.

module import_piece(indices){
   filename = str("br_piece_smooth_",indices[0],"-",indices[1],"-",indices[2],".stl");
   echo(filename);
   import(filename);
}

include<plugs_and_sockets.scad>


