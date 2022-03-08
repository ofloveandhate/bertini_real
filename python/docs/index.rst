.. Test documentation master file, created by
   sphinx-quickstart on Mon Feb 25 16:34:57 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome!
**************************************************************

This is the documentation for the `Bertini_real <https://github.com/ofloveandhate/bertini_real>`_ python library for working with numerical cellular decompositions of algebraic curves and surfaces in any number of variables. 

Follow the tutorials below to learn how to plot surfaces 🌈 & export files suitable for 3d printing 🧱!

.. image:: pictures/croissant_collage.png
   :width: 300

As a super-quick-start, after successfully decomposing a surface in Bertini_real:

::

    import bertini_real

    bertini_real.data.gather_and_save()
    decomposition = bertini_real.data.read_most_recent()
    decomposition.plot()


Tutorials ✏️
============

.. toctree::
   :maxdepth: 1
   :caption: Contents:
   
   tutorials/bertini_real
   tutorials/matplotlib
   tutorials/glumpy
   tutorials/mesh_export
   tutorials/anaglypy


Details 📝
================

.. toctree::
   :maxdepth: 14
   :caption: Contents:

   anaglypy
   curve
   data
   decomposition
   dehomogenize
   parse
   plot
   surface
   util
   vertex
   vertextype
   glumpy

Implementation notes
======================

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   todo

Indices and tables 📋
==========================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Contact 📧
=================

If you have any questions or problems, please `submit an issue on github <https://github.com/ofloveandhate/bertini_real/issues>`_!

.. image:: pictures/3dprints_collage2.png
   :width: 300


