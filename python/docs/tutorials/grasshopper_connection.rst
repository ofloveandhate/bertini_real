Snap together pieces in Grasshopper 🧩
===========================================================================

Attaching plug and sockets to pieces of nodally singular algebraic surfaces.

Exporting complete surfaces
*****************************

After decomposing a surface, you can run the following python codes to export and solidify OBJ in the Python Shell.
First, you should navigate to your desired surface directory via terminal. 
Then using an interactive python shell, execute the following lines.
If you have previously run *"gather()"* on the surface and have a **BRdata.pkl** file for the surface you can skip line 2

::

    import bertini_real as br

    br.data.gather() # do this once after decomposing and sampling the surface. If 

    surface = br.data.read_most_recent()

    pieces = surface.separate_into_nonsingular_pieces()

Next we need to write to a file all the information required to orient the connectors. 
::

    surf.write_piece_data()

The surf.write_piece_data function positions and orients the plug and socket and writes that data to a .scad file and a .json file. This is accomplished by:

1. Finding each singularities that have two connecting pieces
2. for each of the connected pieces found. determine the centroid
3. the centroid and determines a vector to place the plug and socket along
4. write the following data to **br_surf_piece_data.json** and **br_surf_piece_data.scad**

    * piece indices
    * singularities on a piece
    * direction of the vector
    * singularity coordinates
    * and the parity of each piece(whether it gets a plug or socket)

The .json and .scad files are very similar and contain only syntax difference. The .json file is used for adding connectors in grasshopper and we can ignore the .scad file for now.

Awesome! You are now ready to *hop* over to grasshopper

Viewing Connectors in Rhino3D and Grasshopper
*****************************

Launch Rhino3D and navigate to a Grasshopper editor. When creating connectors for a surface there are 6 primary componets required. They can be found in the plug_socket_library.

    1. Socket Hole
    2. Positive Socket
    3. Plug Hole
    4. Positive Plug
    5. Transform Plug and Socket
    6. PieceJson

It is reccommended that you add them as User Components to have easy access.
The following steps are all done in plug_sock_library.gh and if you have access to the file, I reccomend copying all the slider variables.
Components 1-4 produce single connector parts as brep geometries. *Hole* components are used to take the difference of the surface to create space for their corresponding plug or socket.
Positve components add the plug or socket to the surface.
They taken in an assortment of variables, primarily numbers. By using sliders as inputs you can easily make tweaks to their sizes, positions, and other details.

*base* and *eps* are used by all connector part components. *base* is a plane object set at the origin
Place these 4 components and connect all their inputs. Some input variables are used for more than one component.
Go ahead and try messing with the sliders and figuring out what they change
(Need to insert image)
.. image:: snap_together_pictures/hold
    :width: 300

Next place 2 Transform Plug and Socket. Click and drag the geometry output wire from the Positive Socket Componet to the **pos** input of one of the Transform Plug and Socket components
Drag the wire of from the output geometry of Socket Hole to the **neg** input of the same Transform componet that you just connected the Positive Socket to.
Repeat this for the second Transform component, connecting the geometries of Positive Plug, and Plug Hole to the **pos** and **neg** input of the transform respectively.

.. image:: snap_together_pictures/hold
    :width: 300

Transform Plug and Socket takes in the positive and negative parts of a connector and moves them according to inputed directions, locations, and sizes. 
Create a slider and connect it to the size of both Transform components.
If it recieves a list of directions and locations then it will create a connector for each item in the list.

The PieceJson component reads in **br_surf_complete.json** and sends out locations and direction information for each singularity. 
create a Path component, right click it, and select *"Select one existing file"* navigate to the surface that we created the .json file for earlier.
Select **br_surf_piece_data.json**
Then drag the wire of the path component into "filename" input of the PieceJson component.

.. image:: snap_together_pictures/hold
    :width: 300

Finally, connect the outputs of the PieceJson, to the Transform component. 
Connect the Socket_dir output to the direction input for the transform component that is connected to the socket geometries

.. image:: snap_together_pictures/hold
    :width: 300

Amazing! You can now take a look at Rhino3d, you will see an arrangment of plugs and sockets place were the singularites would be!
You may need to disable *preview* on each of the modules, except for the ones wth the final geometry output. 

.. image:: snap_together_pictures/hold
    :width: 300

To do work 🚧
*****************************

    1. Connect with components with surface
    2. Have a better way of organizing the inputs
    3. Investigate Blocks


:Author:
	Caden Joergens

:Version: 1.0 2023/10/12
