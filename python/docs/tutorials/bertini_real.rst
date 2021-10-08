
Introduction to Bertini_real's Python Visualization Suite
==========================================================

Bertini_real is a software for real algebraic sets. There are two options available to decompose and visualize algebraic curves & surfaces either using (1) Matlab or (2) Python. In this documentation, we will focus on the Python visualization suite.

Decomposition
**************
There are two types of child classes for Decomposition (parent class): (1) Curve and (2) Surface. Bertini_real will return either a Curve or Surface object decomposition data based on the given algebraic equations.

Gather and save decomposition
******************************

Use an interactive Python shell, such as iPython.  Start the instance of iPython from the terminal, when you are already in the folder containing the data you want to process.  For example, if you have a folder called `data/whitney/`, and in that folder is your bertini_real decomposition, then move to that folder and THEN launch iPython.

One you're in, you can run the following python codes to gather and save a decomposition of a curve/surface.

::

    import bertini_real

    bertini_real.data.gather_and_save()

It is only necessary to `import bertini_real` once per session.

**Caution** Using tab completion in a shell will not place `()` at the end of functions, make sure to include them.

If you have multiple BRdata*.pkl files, you can retrieve the most recent. This can be done by typing:

::

    data = bertini_real.data.read_most_recent()

Plot curves and surfaces
*************************

There are two plotting modules (Plotter and `GlumpyPlotter <glumpy.html>`_) in Bertini_real Python visualization suite. In this example, we are going to demonstrate how to plot it using Plotter object with matplotlib.

::

    bertini_real.plot.plot()

Here are some of the plotted curve/surface examples

**Alpha Curve** (Curve)

.. image:: bertini_real_pictures/alphacurve_plotter.PNG
   :width: 600

**Dingdong** (Surface)

.. image:: bertini_real_pictures/dingdong_plotter.PNG
   :width: 600

Checking **Smooth STL** will generate and export an 3D object to your current folder. Checking **Raw STL** will do the same.
**Caution** The plot pop-up screen must be closed to use the shell again.

Separate surfaces into pieces
******************************
We are working on the solidification feature for exporting singular algebraic surfaces STL in Bertini_real. We created a Piece object in Fall 2019 to separate surfaces into nonsingular pieces. In this example, we are going to separate surface **Dingdong** into nonsingular pieces.

If you have not yet done so in the current session you must import bertini_real and retrieve the most recent BRdata*.pkl:
::

    import bertini_real

    data = bertini_real.data.read_most_recent()


To separate into non singular pieces:
::

    pieces = data.separate_into_nonsingular_pieces()



We can print out the piece and it should return 2 pieces for **Dingdong**  with its corresponding lists of indices:

::

    print(pieces)


will output:

::

    [piece with indices:[0, 1, 2, 3, 5, 6]
    , piece with indices:[4, 7, 8]
    ]

We can access each piece by specifying their indices. For example,

::

    print(pieces[0])

will output

::

    piece with indices:[0, 1, 2, 3, 5, 6]


Properties and functions of a Piece object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are a few properties and callable functions in a Piece object. You can run this command to generate a list of properties/functions for a Piece by running
::

    dir(pieces[0])

you will get the following output:
::

    ['__doc__',
    '__init__',
    '__module__',
    '__repr__',
    '__str__',
    'indices',
    'is_compact',
    'point_singularities',
    'surface']


To access the indices of a Piece object, type
::

    pieces[0].indices

to output:

::

    [0, 1, 2, 3, 5, 6]


To check whether a Piece object is compact, type
::

  pieces[0].is_compact()

and it'll output:

::

    True


To retrieve the list of point singularities from a Piece object, type:
::

    pieces[0].point_singularities()

and we get:

::

    [0]



There are three modules used to plot surfaces & export stereolithography and 3d animations. To learn more, check out the following modules:

* `Anaglypy <anaglypy.html>`_ (A module that exports 3d anaglyph/non-anaglyph animations of algebraic surfaces)
* `GlumpyPlotter <glumpy.html>`_ (A module that plot curves/surfaces using Glumpy)
* `Tmesh <tmesh.html>`_ (A module that export stereolithography of surfaces for 3d printing using Trimesh)

:Author:
	Foong Min Wong, Caden Joergens

:Version: 1.2 2021-10-08
