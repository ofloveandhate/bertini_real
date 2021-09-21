.. Test documentation master file, created by
   sphinx-quickstart on Mon Feb 25 16:34:57 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

bertinireal.vertextype
================================

.. automodule:: vertextype
   :members:

There are eleven VertexTypes of algebraic curves/surfaces in this module, which are:

* unset = 0
* critical = 1
* semicritical = 2
* midpoint = 4
* isolated = 8
* new = 16
* curve_sample_point = 32
* surface_sample_point = 64
* removed = 128
* problematic = 256
* singular = 512

We use bit operations to indicate the various types, so that a vertex can be multiple types at once.