



///////////////////////////////




electonics_parts_toggle = true;




wire_hole_dia = 7;


socket_cyl_dia = 15;
socket_wall_thickness = 2;
socket_length = 10;

plug_body_overlap = 10; // how much the connection extends into the body.  this is because the points are there...
connection_play = 0.2; // a gap between the plug and the socket

plug_neg_factor = 3;
socket_neg_factor = 3;


plug_taper_factor = .7;
plug_taper_length = 5;


plug_length_overage = 1; // a manual adjustment for making the plug long enough to snap in... this really should be automatically computed..

plug_wedge_h = 2; // how much the tab protrudes above the cylinder.
plug_wedge_w = 2; // how wide the tab is.

plug_tab_depth = 4; // also how long the tab wedge is.  this distance produces overage on how long the plug cylinder is, how long it extends beyond the socket.  
plug_tab_thickness = 1.5;
plug_tab_cutout_depth = 10; // includes the plug_tab_depth
plug_tab_cutout_thickness = 1.2; // the amount of material missing. the gap between the tab and the cylinder remnants




ps_fn = 30;

// preview only, does not affect render
socket_color = "red";
plug_color = "blue";