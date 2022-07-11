// brings in necessary modules, functions, and parameters
// this file is used in `br_surf_complete.scad`
//
// depends on `plugs_and_sockets.scad`.

module import_piece(filename){
   echo("importing",filename);
   import(filename);
}

include<plugs_and_sockets.scad>


