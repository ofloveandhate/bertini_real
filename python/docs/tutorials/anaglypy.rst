Introduction to Anaglypy with bertini_real
===========================================

The Anaglypy class in bertini_real allows user to render 3d anaglyph movies of algebraic surfaces using Blender/Python API. 

Python Scripting
*****************

After decomposing a surface, you can either run the rendering process manually from command-line or automate it through a shell script, ``anaglypy.sh`` (can be found in bertini_real's python ``anaglypy`` folder) to export 3d stereoscopic movies. Please make sure you are using Python 3.7 (matches with the python version of Blender), have already installed Blender 2.8 and added it to the environmental path.

We are using a surface **"Crixxi"** in this example.

Running from command-line
**************************

Copy ``anaglypy.py`` into your surface folder (e.g., ``crixxi/``) folder. The directory structure will look like the following:

::

    crixxi/
      ├── anaglypy.py
      ├── BRdata0.pkl
      ├── input
      ├── witness_data
      └── ...

Then, run this command inside your surface folder:

::

    $ blender -b P anaglypy.py

You can choose different options from the command-line menu for the video:

::

    Please choose:
    1) Raw
    2) Smooth
    3) Both
    Enter number: 

    reading from file BRdata0.pkl

    Please choose:
    1) Rotate Z
    2) Rotate XYZ
    3) Spin
    4) Multi Rotate
    Enter number: 

    Please choose:
    1) Anaglyph 3D
    2) Non-Anaglyph 3D

It will begin the rendering process, for example:

::

    Created avi: render/crixxi_rotate_z_both_anaglyph0000-0100_L.avi

    Fra:0 Mem:25.35M (0.00M, Peak 25.37M) | Time:00:00.00 | Mem:0.00M, Peak:0.00M | Scene, View Layer, left | Synchronizing object | crixxi

    Fra:0 Mem:25.71M (0.00M, Peak 25.71M) | Time:00:00.01 | Mem:0.00M, Peak:0.00M | Scene, View Layer, left | Initializing

    ....

Finally, you can find the video and blender file in a newly created folder, ``render/``.

::

    Export Anaglyph 3D crixxi_rotate_z_both_anaglyph.avi successfully

Here's an example of video output:

.. image:: anaglypy_pictures/crixxi_rotate_z_both_anaglyph.gif
   :width: 50 %
   :align:   center

Using Shell Scripting
**********************

First, create a folder (e.g., ``data/`` ) containing all subfolders of surfaces that you have already decomposed by bertini_real.

Copy ``anaglypy.py`` and ``anaglypy.sh`` located in bertini_real's ``anaglypy`` python folder to the ``data/`` folder you just created.

Remember to change permission of the shell script:

::

    $ chmod 755 anaglypy.sh


You will need **four** files to automate the video rendering process:

1. anaglypy.py (A Blender Python API script)
2. anaglypy.sh (An automation video rendering shell script)
3. surfaces.txt (A texfile containing the names of surfaces)
4. options.txt (A textfile specifying the rendering options)

Example: ``surfaces.txt``
++++++++++++++++++++++++++
::

    crixxi
    daisy


Options choices
++++++++++++++++
The first number specifies the style of surfaces you want to render:

1) Raw
2) Smooth
3) Both

The second number indicates the animation modes: 

1) Rotate Z
2) Rotate XYZ
3) Spin
4) Multi-Rotate (Only available if the first number is 1 or 2)

The third number turns on/off the stereographic 3D mode: 

1) Anaglyph 3D
2) Non-Anaglyph 3D

Example: ``options.txt``
++++++++++++++++++++++++++
::

    3
    1
    2


Place your files in the following standard structure, for example:

::

    data/
      ├── anaglypy.py
      ├── anaglypy.sh
      ├── surfaces.txt
      ├── options.txt
      ├── crixxi/
      |     ├── BRdata0.pkl
      |     └── ...
      └── daisy/
            ├── BRdata0.pkl
            └── ...

Run this command in the terminal to automate the rendering process:
::

	$ ./anaglypy.sh

It will begin multiple rendering processes accordingly based on the ``surfaces.txt``:

::

    Created avi: render/crixxi_rotate_z_both_non_anaglyph0000-0100_L.avi

    Fra:0 Mem:25.35M (0.00M, Peak 25.37M) | Time:00:00.00 | Mem:0.00M, Peak:0.00M | Scene, View Layer, left | Synchronizing object | crixxi

    ....

    Export Anaglyph 3D crixxi_rotate_z_both_non_anaglyph.avi  successfully

    Created avi: render/daisy_rotate_z_both_non_anaglyph0000-0100_L.avi

    Fra:0 Mem:25.35M (0.00M, Peak 25.37M) | Time:00:00.00 | Mem:0.00M, Peak:0.00M | Scene, View Layer, left | Synchronizing object | crixxi

    ....

    Export Anaglyph 3D daisy_rotate_z_both_non_anaglyph.avi  successfully

These are the exported videos from this example:

.. image:: anaglypy_pictures/crixxi_rotate_z_both_non_anaglyph.gif
   :width: 49 %

.. image:: anaglypy_pictures/daisy_rotate_z_both_non_anaglyph.gif
   :width: 49 %

Finally, you can find all videos and blender files in a newly created folder ``render/`` in each surfaces subfolder.

Change object colors
*********************
You can modify the rgb values from this line of code in ``anaglypy.py``:

::

    r,g,b = 1.0, 0.0, 0.2

You can convert videos into gifs and here's a gallery of 3d anaglyph and non-anaglyph animations!

.. image:: anaglypy_pictures/both_multi.gif
   :width: 100 %

.. image:: anaglypy_pictures/schneeflocke_raw_smooth.gif
   :width: 100 %


:Author:
	Foong Min Wong

:Version: 1.0 2019/07/18
