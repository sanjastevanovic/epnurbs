# epnurbs.createshading module contains a method for creating approximation of a NURBS shading in EnergyPlus files

from eppy import modeleditor
from eppy.modeleditor import IDF

from .helper_methods import crossproduct, normalize, distance, foot

# the method for creating approximation of a NURBS shading
def createnurbsshading(idd_filename, idf_filename, base_surface, shading_str, ctrl_points, evaluated_points=20):
    ##########################################
    # loading base surface data from idf file 
    ##########################################
    
    # check if IDD has been already set
    try:
        IDF.setiddname(idd_filename)
    except modeleditor.IDDAlreadySetError as e:
        pass

    # load idf file into idf collection
    idf_collection = IDF(idf_filename)

    # find the named base_surface in idf collection
    walls = idf_collection.idfobjects['BuildingSurface:Detailed'.upper()]
    for wall in walls:
        if wall.Name == base_surface:
            # named base_surface is now contained in wall
            break
    else:
        # named base_surface was not found in idf file
        print('epnurbs.createshading: unable to find the base surface', base_surface, 'in', idf_filename)
        return

    #################################
    # calculating NURBS curve points
    #################################
    from geomdl import NURBS
    from geomdl import utilities

    curve = NURBS.Curve()

    # 1/delta corresponds to the number of trapezoids used in approximation of NURBS shading
    curve.delta = 1/evaluated_points
    curve.degree = 3

    # add weight 1 to each control point
    # unless it has been already weighted
    for cpt in ctrl_points:
        if len(cpt)<4:
            cpt.append(1.0)

    # sets curve control points
    curve.set_ctrlpts(ctrl_points)

    # sets curve knot vector
    curve.knotvector = utilities.generate_knot_vector(curve.degree, len(curve.ctrlpts))

    # evaluates curve points
    curve.evaluate()

    # make a local copy of evaluated curve points
    crv_points = curve.curvepts

    ###################################################################################
    # feet of perpendiculars from the remaining NURBS curve points to the base surface
    ###################################################################################

    # getting the found base_surface's coordinates
    coord = wall.coords
    ulc = coord[0]
    blc = coord[1]
    brc = coord[2]

    # find the base_surface's plane normal and normalize it
    N = crossproduct( (blc[0]-ulc[0], blc[1]-ulc[1], blc[2]-ulc[2]), \
                      (brc[0]-ulc[0], brc[1]-ulc[1], brc[2]-ulc[2]) )
    N = normalize(N)

    # calculate feet of perpendiculars
    feet_points = [foot(crv_points[i], ulc, N) for i in range(len(crv_points))]

    #################################################################################
    # create string of idf definitions for trapezoids that approximate NURBS shading
    #################################################################################
    idf_total_shading_def = ""
    for i in range(1, len(crv_points)):
        # the width of a trapezoid must be at least 0.01
        if distance(crv_points[i-1], crv_points[i])>=0.01 and \
           distance(feet_points[i-1], feet_points[i])>=0.01:
            # are trapezoid arms at least 0.01 or do we have a triangle?
            if distance(crv_points[i-1], feet_points[i-1])>=0.01:
                if distance(crv_points[i], feet_points[i])>=0.01:
                    # both arms are at least 0.01, so we have a trapezoid
                    vertices_str = "{:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}".format(
                                   feet_points[i-1][0], feet_points[i-1][1], feet_points[i-1][2],
                                   crv_points[i-1][0],  crv_points[i-1][1],  crv_points[i-1][2],
                                   crv_points[i][0],    crv_points[i][1],    crv_points[i][2],
                                   feet_points[i][0],   feet_points[i][1],   feet_points[i][2])
                    countervertices_str = "{:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}".format(
                                          feet_points[i][0],   feet_points[i][1],   feet_points[i][2],
                                          crv_points[i][0],    crv_points[i][1],    crv_points[i][2],
                                          crv_points[i-1][0],  crv_points[i-1][1],  crv_points[i-1][2],
                                          feet_points[i-1][0], feet_points[i-1][1], feet_points[i-1][2])
                    single_shading_def = shading_str.replace('<IDX>', str(i)).replace('<BASESURFACE>', base_surface).replace('<VERTICES>', vertices_str).replace('<COUNTERVERTICES>', countervertices_str)
                    idf_total_shading_def = idf_total_shading_def + single_shading_def
                else: 
                    # arm i is less than 0.01, so we have a triangle
                    vertices_str = "{:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}".format(
                                   feet_points[i-1][0], feet_points[i-1][1], feet_points[i-1][2],
                                   crv_points[i-1][0],  crv_points[i-1][1],  crv_points[i-1][2],
                                   feet_points[i][0],   feet_points[i][1],   feet_points[i][2])
                    countervertices_str = "{:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}".format(
                                          feet_points[i][0],   feet_points[i][1],   feet_points[i][2], 
                                          crv_points[i-1][0],  crv_points[i-1][1],  crv_points[i-1][2],
                                          feet_points[i-1][0], feet_points[i-1][1], feet_points[i-1][2])               
                    single_shading_def = shading_str.replace('<IDX>', str(i)).replace('<BASESURFACE>', base_surface).replace('<VERTICES>', vertices_str).replace('<COUNTERVERTICES>', countervertices_str)
                    idf_total_shading_def = idf_total_shading_def + single_shading_def
            else:
                # arm i-1 is less than 0.01, but do we still have a triangle?
                if distance(crv_points[i], feet_points[i])>=0.01:
                    # we have a triangle
                    vertices_str = "{:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}".format(
                                   feet_points[i-1][0], feet_points[i-1][1], feet_points[i-1][2],
                                   crv_points[i][0],    crv_points[i][1],    crv_points[i][2],
                                   feet_points[i][0],   feet_points[i][1],   feet_points[i][2])
                    countervertices_str = "{:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}".format(
                                          feet_points[i][0],   feet_points[i][1],   feet_points[i][2],
                                          crv_points[i][0],    crv_points[i][1],    crv_points[i][2],
                                          feet_points[i-1][0], feet_points[i-1][1], feet_points[i-1][2])              
                    single_shading_def = shading_str.replace('<IDX>', str(i)).replace('<BASESURFACE>', base_surface).replace('<VERTICES>', vertices_str).replace('<COUNTERVERTICES>', countervertices_str)
                    idf_total_shading_def = idf_total_shading_def + single_shading_def
                else:
                    # we do not have a shading element in this case
                    pass
    
    # create idf shading objects from the string containing shading definitions   
    from io import StringIO
    idf_shading = IDF(StringIO(idf_total_shading_def))

    # copy idf shading objects to the existing idf file
    shadings = idf_shading.idfobjects["SHADING:ZONE:DETAILED"]
    for shading in shadings:
        idf_collection.copyidfobject(shading)

    ###############################
    # at the end, save the changes
    ###############################
    idf_collection.save()
