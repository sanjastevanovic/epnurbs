.. epnurbs documentation master file, created by
   sphinx-quickstart on Thu May 10 22:02:25 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to epnurbs's documentation!
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

epnurbs package serves to create shading surfaces with nurbs outer edge in EnergyPlus .idf/.imf files.
It relies on *eppy* package to handle .idf/.imf files and 
*NURBS-Python* to create points on a nurbs curve.
It can be installed with the usual

  ``pip install epnurbs``

Consider first an example of using epnurbs::

  ###################################################
  # jEPlus supplies the following arguments:
  #   sys.argv[1]  -  project's base folder where the project files are located
  #   sys.argv[2]  -  output folder of the project where the RunTimes.csv is located
  #   sys.argv[3]  -  Other arguments specified in the parameter definition. They are passed in as a ',' delimitted string
  #   sys.argv[4]  -  folder of the binary files of the simulation program, e.g. the location of Energy+.idd

  import sys
  import os

  # path to E+ idd file, required by eppy
  idd_filename = os.path.join(sys.argv[4], 'Energy+.idd')

  # path to idf input file within each simulated folder
  idf_filename = os.path.join(sys.argv[2], 'in.idf')

  ###################################################
  # idf template for shading elements
  shading_str_top  = 'Shading:Zone:Detailed, ShadingTop<IDX>, <BASESURFACE>, , , <VERTICES>;'

  # alternative values for NURBS control points
  start_values = [ [0.3, 0, 0.6], [0.3, 0, 1.1], [0.3, 0, 1.6], [0.3, 0, 2.1], [0.3, 0, 2.6],
                 [0.8, 0, 2.6], [1.3, 0, 2.6], [1.8, 0, 2.6], [2.3, 0, 2.6], [2.8, 0, 2.6],
                 [3.3, 0, 2.6], [3.3, 0, 2.1], [3.3, 0, 1.6], [3.3, 0, 1.1], [3.3, 0, 0.6] ]
  ctrl_point_alternatives = [ [ [start_values[i][0], start_values[i][1] - 0.25*j, start_values[i][2]]
                                 for j in range(9)]
                               for i in range(len(start_values))]

  # parse the jEPlus parameters from sys.argv[3]
  param = sys.argv[3].split(',')

  # param[0], param[1], ... are strings representing control point indices from ctrl_point_alternatives list
  indices = [int(param[i]) for i in range(len(param))]
  actual_ctrl_points_top  = [ctrl_point_alternatives[i][indices[i]] for i in range(0, 7)]

  ###################################################
  # call the appropriate epnurbs method
  import epnurbs

  epnurbs.createnurbsshading(idd_filename, idf_filename, 'ZidJug', shading_str_top,  actual_ctrl_points_top, evaluated_points=15)

It is assumed here that one is using jEPlus for …
