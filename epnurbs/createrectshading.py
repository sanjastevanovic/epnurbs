# epnurbs.createrectshading module contains a method for creating a sequence of rectangular shadings in EnergyPlus files

from eppy import modeleditor
from eppy.modeleditor import IDF

from .helper_methods import crossproduct, normalize, distance, foot

# the method for creating a sequence of rectangular shadings
# between the start and the end points with depths given in the list depths
# the start point and the end point do not need to actually belong to the base_surface
# as they are projected on it first
def createrectshading(idd_filename, idf_filename, base_surface, shading_str, start_point, end_point, depths):
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
        print('epnurbs.createrectshading: unable to find the base surface', base_surface, 'in', idf_filename)
        return

    ##############################################################################
    # feet of perpendiculars from the start and the end point to the base surface
    ##############################################################################

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
    start_foot = foot(start_point, ulc, N)
    end_foot = foot(end_point, ulc, N)

    # steps between the start and the end point
    num = len(depths)
    step = [(end_foot[0]-start_foot[0])/num, (end_foot[1]-start_foot[1])/num, (end_foot[2]-start_foot[2])/num]

    ########################################################################
    # create string of idf definitions for a sequence of shading rectangles
    ########################################################################
    idf_total_shading_def = ""
    for i in range(num):
        if depths[i]>0.01:
            # two base_surface vertices are start_foot + i*step and start_foot + (i+1)*step
            # the other two vertices are at depth[i] distance from the base surface
            v1 = [start_foot[0]+i*step[0], start_foot[1]+i*step[1], start_foot[2]+i*step[2]]
            v2 = [v1[0]+depths[i]*N[0], v1[1]+depths[i]*N[1], v1[2]+depths[i]*N[2]]

            v4 = [start_foot[0]+(i+1)*step[0], start_foot[1]+(i+1)*step[1], start_foot[2]+(i+1)*step[2]] 
            v3 = [v4[0]+depths[i]*N[0], v4[1]+depths[i]*N[1], v4[2]+depths[i]*N[2]]

            vertices_str = "{:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}".format(
                            v1[0], v1[1], v1[2],
                            v2[0], v2[1], v2[2],
                            v3[0], v3[1], v3[2],
                            v4[0], v4[1], v4[2])
            countervertices_str = "{:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}".format(
                                    v4[0], v4[1], v4[2],
                                    v3[0], v3[1], v3[2],
                                    v2[0], v2[1], v2[2],
                                    v1[0], v1[1], v1[2])
            single_shading_def = shading_str.replace('<IDX>', str(i)).replace('<BASESURFACE>', base_surface).replace('<VERTICES>', vertices_str).replace('<COUNTERVERTICES>', countervertices_str)
            idf_total_shading_def = idf_total_shading_def + single_shading_def
        else:
            # we do not have a shading rectangle if it is not deep enough
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
