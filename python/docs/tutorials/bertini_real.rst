
Introduction to Bertini_real's Python Visualization Suite
==========================================================

Bertini_real is a software for real algebraic sets. There are two options available to decompose and visualize algebraic curves & surfaces either using (1) Matlab or (2) Python. In this documentation, we will focus on the Python visualization suite.

Python Scripting
*****************

Use an interactive Python shell, such as IPython, and you can run the following python codes to gather and save a decomposition of a curve/surface.

::

    import bertini_real

    bertini_real.data.gather_and_save()

    surface = bertini_real.data.read_most_recent()

There are several modules in Bertini_real's Python visualization suite. To learn more, feel free to check out the following modules:

* `Anaglypy <anaglypy.html>`_
* `GlumpyPlotter <glumpy.html>`_
* `Tmesh <tmesh.html>`_




:Author:
	Foong Min Wong

:Version: 1.1 2019/12/08
.. :Version: 1.0 2019/04/22
