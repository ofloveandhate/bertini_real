Snap together pieces
===========================================================================

Attaching plug and sockets to pieces of nodally singular algebraic surfaces.

Exporting complete surfaces
*****************************

After decomposing a surface, you can run the following python codes to export and solidify OBJ in the Python Shell.
We are using a surface **"octdong"** in this example .

::

    import bertini_real as br

    br.data.gather() # do this once after decomposing and sampling the surface.

    surface = br.data.read_most_recent()

    pieces = surface.separate_into_nonsingular_pieces()
