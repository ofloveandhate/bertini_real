Introduction to Matplotlib with bertini_real
===========================================

The plot class in bertini_real uses matplotlib that allows you to render the decompositions
computed by bertini_real.

Python Scripting
******************
After decomposing a surface, you can run the following python codes to plot curve/surface in the Python Shell. 
We are plotting a curve **Alpha Curve** and a surface  **"Whitney"** in this example ☂️ .


Example: Alpha Curve
*********************
::

    import bertini_real

    bertini_real.data.gather_and_save()

    curve = bertini_real.data.read_most_recent()

    bertini_real.plot.plot()


Example: Whitney
*****************
::

    import bertini_real

    bertini_real.data.gather_and_save()

    surface = bertini_real.data.read_most_recent()

    bertini_real.plot.plot()






:Author:
	Foong Min Wong

:Version: 1.0 2020/03/01
