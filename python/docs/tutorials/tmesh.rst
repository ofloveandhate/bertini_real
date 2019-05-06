Introduction to TMesh with bertini_real
========================================

The TMesh class in bertini_real allows user to export and solidify (raw & smooth surfaces) stereolithography (STL) for 3D printing.

Python Scripting
*****************

After decomposing a surface, you can run the following python codes to export and solidify STL in the Python Shell. 
We are using a surface **"Whitney"** in this example ☂️ .

::

    import bertini_real

    bertini_real.data.gather()

    surface = bertini_real.data.ReadMostRecent()

To export the raw version of STL, you can type:

::

	bertini_real.tmesh.stl_raw(surface)

And you will see an output like this:

::

	bertini_real.tmesh.stl_raw(surface)

	Generating raw STL surface...

	Export stl_raw_whitney.stl successfully

You can display STL using any 3D STL viewer:

.. image:: tmesh_pictures/stl_raw_whitney.PNG
   :width: 300

To export the smooth version of STL, you can type:

::

	bertini_real.tmesh.stl_smooth(surface)

And you will see an output like this:

::

	bertini_real.tmesh.stl_smooth(surface)

	Generating smooth STL surface...
	Export stl_smooth_whitney.stl successfully

You can display STL using any 3D STL viewer:

.. image:: tmesh_pictures/stl_smooth_whitney.PNG
   :width: 300

To solidify the raw version of STL, you can type:

::

	bertini_real.tmesh.solidify_raw(surface)

And you will see an output like this:

::

	bertini_real.tmesh.solidify_raw(surface)

	Solidiying raw STL surface...
	Export solidify_raw_whitney.stl successfully

You can display STL using any 3D STL viewer:

.. image:: tmesh_pictures/solidify_raw_whitney.PNG
   :width: 300

To solidify the smooth version of STL, you can type:

::

	bertini_real.tmesh.solidify_smooth(surface)

And you will see an output like this:

::

	bertini_real.tmesh.solidify_smooth(surface)

	Solidiying smooth STL surface...
	Export solidify_smooth_whitney.stl successfully

You can display STL using any 3D STL viewer:

.. image:: tmesh_pictures/solidify_smooth_whitney.PNG
   :width: 300

:Author:
	Foong Min Wong

:Version: 1.0 2019/04/22
