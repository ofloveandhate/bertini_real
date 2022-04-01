// this file puts together the singularly-separated pieces
// of a nodally singular algebraic surface computed in
// bertini_real.
//
// depends on:
// * `br_surf_piece_data.scad` -- automatically generated using the Python library for bertini_real
// * `br_modules.scad` -- should be copied into your surface folder via py bertini_real.
// * `plugs_and_sockets.scad` -- should be copied into your surface folder.  no parameters in here.
// * `plug_socket_params.scad` -- should be copied into your surface folder.  change as you need to effect the shape and size of the plugs and sockets.
// * `br_piece_smooth_A_B_C.stl` -- generated stl files for each piece of the surface.  A,B,C are integer indices into the faces of the surface.  each piece is a union of faces -- A,B,C are the first three.




// silviana amethyst
// university of wisconsin-eau claire
// spring 2022

include<br_surf_piece_data.scad>; // AUTOMATICALLY GENERATED





include<br_modules.scad>; // not automatically generated, but copied in per-surface



// a few parameters
//
// adjust to make pieces skewed from each other.  0 is in-place as on the surface, so probably not suitable for printing.
piece_offset_x = 0.;
piece_offset_y = 0.;
piece_offset_z = 0.;
cut = true; // cutouts help understand the inside of the surface



difference(){ // the difference is for the cutout.
for (ii=[0:len(piece_indices)-1]){

    indices = piece_indices[ii]; // this is a Piece property

    sings = singularities_on_pieces[ii];

    translate([ii*piece_offset_x,piece_offset_y,ii*piece_offset_z]) // to offset pieces from each other
    Piece(indices, ii, parities, sings, sing_directions, sing_locations, conn_size);

}
    if (cut){
        translate([0,0,-50])
        cube([100,100,100]);
    }
}

module Piece(indices, piece_index, parities, sings, directions, locations, conn_size){
    echo(indices)


    difference(){
        import_piece(indices);

        for (ii=[0:len(sings)-1]){

        sing_index = sings[ii];

        // for each singularity on the piece
        direction = directions[sing_index];
        location = locations[sing_index];
        parity_this_sing = parities[sing_index][piece_index];

        if (parity_this_sing==-1){
            socket_neg(direction, location, conn_size);
        }
        else if (parity_this_sing==1){
            plug_neg(direction, location, conn_size);
        }

        }
    }

    for (ii=[0:len(sings)-1]){

        sing_index = sings[ii];
        echo(sing_index);
        // for each singularity on the piece
        direction = directions[sing_index];
        location = locations[sing_index];
        parity_this_sing = parities[sing_index][piece_index];

        if (parity_this_sing==-1){
            socket_pos(direction, location, conn_size);
        }
        else if (parity_this_sing==1){
            plug_pos(direction, location, conn_size);
        }

       }

}
