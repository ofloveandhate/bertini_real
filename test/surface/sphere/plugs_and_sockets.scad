include<plug_socket_params.scad>;


//d = [0,0,1];
//
//size = 1;
//
//p = 20;
//
//socket_pos(d,[0,0,0],size);
//
//translate([0,p,0])
//socket_neg(d,[0,0,0],size);
//
//translate([p,0,0])
//plug_pos(d,[0,0,0],size);
//
//translate([p,p,0])
//plug_neg(d,[0,0,0],size);



 

//////////////////////////////




// we next set up some helper functions
eps = 0.01; // helps avoid coplanar faces, for doing booleans

function theta(direction) = atan2(direction[1], direction[0]); 
function phi(direction) = acos(direction[2] / norm(direction)); 

module put_sp_in_place(direction, location, size){
    translate(location)
    rotate(theta(direction), [0,0,1])
    rotate(phi(direction), [0,1,0])
    scale([size,size,size])
        children();
}








// now for the top-level code which actually makes plugs and sockets


module plug_neg(direction, location, size, ps_params){

    color(plug_color)
    put_sp_in_place(direction, location, size)
    {
        cutting_r = wire_hole_dia/2;
        h = plug_neg_factor*(socket_length+plug_length_overage+plug_body_overlap);

        translate([0,0,-eps])
            cylinder(r=cutting_r, h=h,center=true,$fn=ps_fn);
    }
    
}





module plug_pos(direction, location, size, ps_params){


    
    color(plug_color)
    put_sp_in_place(direction, location, size)
    {
  
  plug_r = socket_cyl_dia/2-socket_wall_thickness - connection_play;
	
  tapered_r = plug_taper_factor * plug_r;
  


  translate([0,0,-(socket_length+plug_length_overage)]) 
  
    union(){
      
      difference()
      {
        union() // the total length of this should be socket_length + plug_body_overlap + plug_length_overage
        {
          translate([0,0,plug_taper_length])
          cylinder(r=plug_r, h=socket_length+plug_body_overlap - plug_taper_length+plug_length_overage,$fn=ps_fn); // the main body of the plug
          cylinder(r2=plug_r, r1= tapered_r, h=plug_taper_length,$fn=ps_fn); // the tapered tip of the plug
        }
        
        
        if (electonics_parts_toggle)
        {
            cutting_r = wire_hole_dia/2;
            if (cutting_r >= plug_r)
              echo("check your dims, cutting radius", cutting_r, "larger than plug radius",plug_r);
            plug_neg(direction, location, size);
        }
        
        
        
        // the tab cutouts
        ell = plug_r-plug_tab_thickness - plug_tab_cutout_thickness/2;
        for (ii=[-1,1]) // do it to both sides!
        {
          translate([0,ii*ell,plug_tab_cutout_depth/2])
            cube([plug_r*2,plug_tab_cutout_thickness,plug_tab_cutout_depth+eps],center=true);
        } // for
      } // diff
      
      // add wedges for snaps, to retain
      for (ii=[-1,1]){
        translate([0,ii*( plug_r-plug_tab_thickness ),plug_taper_length])
        rotate(ii*90,[0,0,1])
            wedge(l=plug_tab_depth,w=plug_wedge_w,h=plug_wedge_h+plug_tab_thickness);
      }
      
      
    } // union after translating into place
    
    
    } // put in place
} // plug_pos








module socket_neg(direction, location, size, params){

       
    color(socket_color)
    
    put_sp_in_place(direction, location, size)
    rotate(180,[0,1,0])
    {
        cylinder(r=socket_cyl_dia/2-socket_wall_thickness, h=socket_neg_factor*socket_length+eps,$fn=ps_fn,center=true);
    }
}


module socket_pos(direction, location, size, params){

     color(socket_color)
    put_sp_in_place(direction, location, size)
    rotate(180,[0,1,0])
    {
        difference(){
            cylinder(r=socket_cyl_dia/2, h=socket_length, $fn=ps_fn);
         translate([0,0,-eps])
            cylinder(r=socket_cyl_dia/2-socket_wall_thickness, h=socket_length+2*eps, $fn=ps_fn);
          } // diff
    }
}









// and some sub-pieces, like for tabs, etc


module wedge(l,w,h)
{
	
	rotate(-90,[1,0,0])
  translate([0,0,-w/2])
  linear_extrude(height=w)
  polygon(
	points = [ [0,0], [0,l], [h,0] ], 
	paths = [ [0,1,2,0]], 
	convexity = 1);
}

