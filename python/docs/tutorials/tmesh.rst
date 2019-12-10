Introduction to TMesh with bertini_real
========================================

The TMesh class in bertini_real allows user to export and solidify (raw & smooth surfaces) OBJ filesfor 3D printing.

Python Scripting
*****************

After decomposing a surface, you can run the following python codes to export and solidify OBJ in the Python Shell. 
We are using a surface **"Whitney"** in this example ☂️ .

::

    import bertini_real

    bertini_real.data.gather()

    surface = bertini_real.data.read_most_recent()

To export the raw version of OBJ, you can type:

::

	bertini_real.tmesh.obj_raw(surface)

And you will see an output like this:

::

	bertini_real.tmesh.obj_raw(surface)

	Generating raw OBJ surface...

	Export obj_raw_whitney.obj successfully

You can display OBJ using any 3D OBJ viewer:

.. image:: tmesh_pictures/obj_raw_whitney.PNG
   :width: 300

To export the smooth version of OBJ, you can type:

::

	bertini_real.tmesh.obj_smooth(surface)

And you will see an output like this:

::

	bertini_real.tmesh.obj_smooth(surface)

	Generating smooth OBJ surface...
	Export obj_smooth_whitney.obj successfully

You can display OBJ using any 3D OBJ viewer:

.. image:: tmesh_pictures/obj_smooth_whitney.PNG
   :width: 300

To solidify the raw version of OBJ, you can type:

::

	bertini_real.tmesh.solidify_raw(surface)

And you will see an output like this:

::

	bertini_real.tmesh.solidify_raw(surface)

	Solidiying raw OBJ surface...
	Export solidify_raw_whitney.obj successfully

You can display OBJ using any 3D OBJ viewer:

.. image:: tmesh_pictures/solidify_raw_whitney.PNG
   :width: 300

To solidify the smooth version of OBJ, you can type:

::

	bertini_real.tmesh.solidify_smooth(surface)

And you will see an output like this:

::

	bertini_real.tmesh.solidify_smooth(surface)

	Solidiying smooth OBJ surface...
	Export solidify_smooth_whitney.obj successfully

You can display OBJ using any 3D OBJ viewer:

.. image:: tmesh_pictures/solidify_smooth_whitney.PNG
   :width: 300

:Author:
	Foong Min Wong

:Version: 1.1 2019/09/17
.. :Version: 1.0 2019/04/22
