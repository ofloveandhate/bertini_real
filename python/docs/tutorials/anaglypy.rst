3d anaglyph movies in Blender -- bertinireal.anaglypy
=======================================================

The Anaglypy class in bertini_real allows user to render 3d anaglyph movies of algebraic surfaces using Blender/Python API. 


Preliminaries
**************

1. Please make sure you are using Python 3.7 (matches with the python version of Blender), have already installed Blender 2.8 and added it to the environmental path.

* MacOS: Install python 3.7.0 from python.org
* Linux (Ubuntu): ``sudo apt install python3.7``

2. Install Blender 2.8 (latest version which uses Python 3.7 interpreter)

* Download and install `Blender 2.8 <https://www.blender.org/2-8/>`_

3. Add Blender to PATH:

* MacOS: ``echo "alias blender=/Applications/blender2.8/blender.app" >> ~/.bashrc``

* Linux (Ubuntu): ``export PATH=/path/to/blender2.8/blender:$PATH``

4. Install python packages required by Bertini_real into Blender Python Modules Folder in your local terminal: (use the latest version pip for Python 3.7)

* MacOS: ``pip3 install -t /Applications/blender2.8/blender.app/Contents/Resources/2.80/scripts/modules <python-package-name> --upgrade``

* Linux (Ubuntu): ``pip3 install -t /path/to/blender2.8/2.80/scripts/modules <python-package-name> --upgrade``

* +------------+------------+-----------+-----------+-----------+
  | List of python package names                                |
  +============+============+===========+===========+===========+
  |     dill   | matplotlib |    sympy  | scipy     | numpy     |
  +------------+------------+-----------+-----------+-----------+
  | algopy     | mpmath     | trimesh   | shapely   | cython    |
  +------------+------------+-----------+-----------+-----------+
  | pyopengl   | triangle   | glumpy    | pillow    | encodings |
  +------------+------------+-----------+-----------+-----------+

5. Other issues:

* If the commands don’t work, add option ``--system``

* For Linux users, if you get error: “command 'x86_64-linux-gnu-gcc' failed with exit status 1”, try running this command: ``sudo apt-get install python3.7-dev`` or try this `solution <https://github.com/scrapy/scrapy/issues/2115>`_

*   If you get errors such as “pickle.load  EOF Error: Ran out of input” related to pickle, try running “pip 3 install dill pillow --upgrade”

MacOS Configuration
++++++++++++++++++++

1. Make environment modifications

::

    export BLENDER_SYSTEM_PYTHON=/usr/local/lib/python3.7/site-packages

2. Solve the finding of the scripts problem

::

    export BLENDER_SYSTEM_SCRIPTS=/Applications/Blender.app/Contents/Resources/2.80/scripts

3. Add the entire python path to the one in the shell. The paths to add come from getting the system path from the python console **INSIDE** blender. Type the following commands **in the blender python console**

::

    import sys

    print(sys.path)

4. Clean up the formatting ``sys.path`` and replace with colons

5. Then, in the local terminal, export the path from Step 4 (note the escaped space in that line!):

::

    export PYTHONPATH=/Applications/Blender.app/Contents/Resources/2.80/scripts/addons_contrib:/Applications/Blender.app/Contents/Resources/2.80/scripts/addons:/Applications/Blender.app/Contents/Resources/2.80/scripts/startup:/Applications/Blender.app/Contents/Resources/2.80/scripts/modules:/Applications/Blender.app/Contents/Resources/2.80/python/lib/python37.zip:/Applications/Blender.app/Contents/Resources/2.80/python/lib/python3.7:/Applications/Blender.app/Contents/Resources/2.80/python/lib/python3.7/lib-dynload:/Applications/Blender.app/Contents/Resources/2.80/python/lib/python3.7/site-packages:/Applications/Blender.app/Contents/Resources/2.80/scripts/freestyle/modules:/Applications/Blender.app/Contents/Resources/2.80/scripts/addons/modules:/Users/brakeda/Library/Application\ Support/Blender/2.80/scripts/addons/modules


Python Scripting
*****************

After done configuring and installing Blender, you can decompose a surface, and  automate video rendering process of the surface through a shell script, ``anaglypy.sh`` (can be found in bertini_real's python ``anaglypy`` folder) to export 3d stereoscopic movies. 

We are using surface **Crixxi** and **Daisy** in this example.

Using Shell Scripting
**********************

First, create a folder (e.g., ``data/`` ) containing all subfolders of surfaces that you have already decomposed by bertini_real.

Copy ``anaglpy.json``, ``anaglypy.py`` and ``anaglypy.sh`` located in bertini_real's ``anaglypy`` python folder to the ``data/`` folder you just created.

Remember to change permission of the shell script:

::

    $ chmod 755 anaglypy.sh


You will need **four** files to automate the video rendering process:

1. anaglypy.json (A JSON file specifying object and video settings)
2. anaglypy.py (A Blender Python API script)
3. anaglypy.sh (An automation video rendering shell script)
4. surfaces.txt (A texfile containing the names of surfaces)

Explanation of ``anaglypy.json``
+++++++++++++++++++++++++++++++++

::

  {

    "anaglyph": 0 | 1 # 0: Non-Anaglyph 3D Effect, 1: Anaglyph 3D Effect
    "raw_smooth": "Raw" | "Smooth" | "Both"
    "animation": "Rotate_Z" | "Rotate_XYZ" | "Spin" | "Multi_Rotate" # Multi_Rotate is only available if raw_smooth = "Raw" or "Smooth")]

    "color": {
      "r": 0 - 255 # red
      "g": 0 - 255 # green
      "b": 0 - 255 # blue
      "alpha": 0 - 1 # opacity of color
    },

    "video_settings":{
      "num_frames": 100 # Number of frames 
      "film_exposure": 3.5 # Film exposure
      "file_format": "AVI_JPEG" | "AVI_RAW | "FFMPEG" # Blender Movie File Format 
      "color_mode": "RGB" # Color mode
      "views_format": "STEREO_3D" # For stereocopy movies only
    },

    "resolution":{
      "percentage":100 # Percentage scale of render resolution
      "x": 800 # Resolution X
      "y": 700 # Resolution Y
    },

    "object":{ # Single Object
      "dimension": 1.5
      "inflation":0.15
    },

    "object_both":{ # Two Objects
      "dimension": 1.15
      "inflation":0.15
      "location": 1.15
    },

    "object_multi":{ # Three Objects
      "dimension": 1.0
      "inflation":0.1
      "location": 1.65
    },

    "convergence_distance": 11 # Converge point for the stereo cameras
    "interocular_distance":1.5 # Set the distance between the eyes - the stereo plane distance 
  }



Example: ``surfaces.txt``
++++++++++++++++++++++++++
::

    crixxi
    daisy


Place your files in the following standard structure, for example:

::

    data/
      ├── anaglypy.json
      ├── anaglypy.py
      ├── anaglypy.sh
      ├── surfaces.txt
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

Examples of video output
************************

Here's a gallery of 3d anaglyph and non-anaglyph algebraic surface animations!

.. image:: anaglypy_pictures/both_multi.gif
   :width: 100 %

.. image:: anaglypy_pictures/schneeflocke_raw_smooth.gif
   :width: 100 %


:Author:
	Foong Min Wong

:Version: 1.0 2019/07/18
