# epnurbs.createopening module contains a method for creating approximation of a NURBS opening in EnergyPlus files

from eppy import modeleditor
from eppy.modeleditor import IDF

from .helper_methods import subtract, dotproduct, crossproduct, length, normalize, distance, foot

# Modified geomdl.utilities.check_knot_vector method to
# check if the input knot vector follows the mathematical rules
def new_check_knot_vector(degree=0, knot_vector=(), control_points_size=0, tol=0.001):
    """ Checks if the input knot vector follows the mathematical rules. """    
    if not knot_vector:
        raise ValueError("Input knot vector cannot be empty")

    # Check the formula; m = p + n + 1
    if len(knot_vector) is not degree + control_points_size + 1:
        return False

    # Set up a return value
    ret_val = True

    # Check ascending order
    if ret_val:
        prev_knot = knot_vector[0]
        for knot in knot_vector:
            if prev_knot > knot:
                ret_val = False
                break
            prev_knot = knot

    return ret_val


# the method for creating approximation of a NURBS opening
def createnurbsopening(idd_filename, idf_filename, base_surface, opening_str, ctrl_points, evaluated_points=50):
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
        print('epnurbs.createopening: unable to find the base surface', base_surface, 'in', idf_filename)
        return

    #################################
    # calculating NURBS curve points
    #################################
    from geomdl import NURBS
    from geomdl import utilities

    # one has to exclude checks for leading zeros and trailing ones in geomdl.utilities.check_knot_vector
    # in order for this knotvector to be allowed!
    utilities.check_knot_vector = new_check_knot_vector
    
    curve = NURBS.Curve()

    # 1/delta corresponds to the number of evaluated points of a NURBS curve that defines the opening
    curve.delta = 1/evaluated_points
    curve.degree = 3

    # add weight 1 to each control point
    for cpt in ctrl_points:
        if len(cpt)<4:
            cpt.append(1.0)

    # add copy of the first degree control points in order to close NURBS curve
    extra_points = ctrl_points[:curve.degree]
    ctrl_points.extend(extra_points)

    # sets curve control points
    curve.set_ctrlpts(ctrl_points)

    # sets knot vector for closed NURBS curve
    # note that the length of ctrl_points has increased for curve.degree due to extend operation
    knotstep = 1/(len(ctrl_points) + curve.degree)
    curve.knotvector = [i*knotstep for i in range(len(ctrl_points)+curve.degree+1)]

    # note that for the closed NURBS curve
    # we have to use only the domain [knotvector[curve.degree], knotvector[len(ctrl_points)]]!!!
    # evaluates curve points
    curve.evaluate(curve.degree*knotstep, len(ctrl_points)*knotstep)

    # make a local copy of evaluated curve points
    crv_points = curve.curvepts

    #########################################################################
    # feet of perpendiculars from the NURBS curve points to the base surface,
    # just to make sure we are not getting out of its plane
    #########################################################################

    # getting the found base_surface's coordinates
    coord = wall.coords
    ulc = coord[0]
    blc = coord[1]
    brc = coord[2]

    # find the base_surface's plane normal and normalize it
    u = [ blc[0]-ulc[0], blc[1]-ulc[1], blc[2]-ulc[2] ]
    v = [ brc[0]-blc[0], brc[1]-blc[1], brc[2]-blc[2] ]
    N = crossproduct(u,v)
    N = normalize(N)

    # calculate feet of perpendiculars
    feet_points = [foot(crv_points[i], ulc, N) for i in range(len(crv_points))]

    #######################################################################################
    # create a list of 0.1m squares that partition the whole wall surface
    # !!! it is assumed that the wall surface is a rectangle
    # !!! three of whose vertices are ulc, blc and brc (and the fourth one is ulc+brc-blc)
    #######################################################################################

    # a new orthonormal system for the wall surface
    u = normalize(u)
    v = crossproduct(N,u)
    v = normalize(v)

    # transform feet_point coordinates to the uv system
    # since N, u, v is an orthonormal system we have that
    # for each vector X = <X,N>N + <X,u>u + <X,v>v
    # this will be applied to vectors feet_points-ulc which belong to the wall surface plane
    feet_points_uv = [ [dotproduct(subtract(fp, ulc), u),
                        dotproduct(subtract(fp, ulc), v)] for fp in feet_points]

    # form a list of squares and a list of their centers in uv coordinates
    squaresize = 0.1
    maxi = int(length(subtract(blc, ulc))/squaresize)-1
    maxj = int(length(subtract(brc, blc))/squaresize)-1

    squarelist = [ [ [ ((i+0.5)*squaresize, (j+0.5)*squaresize), ((i+1.5)*squaresize, (j+0.5)*squaresize),
                       ((i+1.5)*squaresize, (j+1.5)*squaresize), ((i+0.5)*squaresize, (j+1.5)*squaresize) ]
                   for j in range(maxj) ]
                 for i in range(maxi)]
    squarecenters = [ ( (i+1)*squaresize, (j+1)*squaresize ) for i in range(maxi) for j in range(maxj) ]

    # for each square check whether its center belongs to the polygon defined by NURBS curve points
    import matplotlib.path as mplPath
    path = mplPath.Path(feet_points_uv)
    inside = path.contains_points(squarecenters)

    # select only those squares whose centers are inside the NURBS curve
    insideindices = [ [ j for j in range(maxj) if inside[i*maxj+j] ] for i in range(maxi) ]

    # now insideindices[i] contains all j such that squarelist[i][j] is inside a NURBS curve

    #################################################################################
    # create string of idf definitions for rectangles that approximate NURBS opening
    #################################################################################
    idf_total_opening_def = ""

    for i in range(maxi):
        # go through all inside squares within i-th column, if there are any
        if len(insideindices[i])>0:
            # pn and kn contain the beginning and the end of
            # subsequences of consecutive numbers within the insideindices[i] list
            pn = insideindices[i][0]
            kn = insideindices[i][0]

            for k in range(1, len(insideindices[i])+1):
                if k<len(insideindices[i]) and insideindices[i][k]==insideindices[i][k-1]+1:
                    # consecutive subsequence goes on
                    kn = insideindices[i][k]
                else:
                    # consecutive subsequence had ended, so create an opening element
                    # squarelist[i][pn][0],[1] contain its first two vertices in uv system,
                    # squarelist[i][kn][2],[3] contain its last two vertices in uv system
                    # recalculate vertices in xyz system first
                    vertex1 = [ ulc[0] + squarelist[i][pn][0][0] * u[0] + squarelist[i][pn][0][1] * v[0],
                                ulc[1] + squarelist[i][pn][0][0] * u[1] + squarelist[i][pn][0][1] * v[1],
                                ulc[2] + squarelist[i][pn][0][0] * u[2] + squarelist[i][pn][0][1] * v[2] ]
                    vertex2 = [ ulc[0] + squarelist[i][pn][1][0] * u[0] + squarelist[i][pn][1][1] * v[0],
                                ulc[1] + squarelist[i][pn][1][0] * u[1] + squarelist[i][pn][1][1] * v[1],
                                ulc[2] + squarelist[i][pn][1][0] * u[2] + squarelist[i][pn][1][1] * v[2] ]
                    vertex3 = [ ulc[0] + squarelist[i][kn][2][0] * u[0] + squarelist[i][kn][2][1] * v[0],
                                ulc[1] + squarelist[i][kn][2][0] * u[1] + squarelist[i][kn][2][1] * v[1],
                                ulc[2] + squarelist[i][kn][2][0] * u[2] + squarelist[i][kn][2][1] * v[2] ]
                    vertex4 = [ ulc[0] + squarelist[i][kn][3][0] * u[0] + squarelist[i][kn][3][1] * v[0],
                                ulc[1] + squarelist[i][kn][3][0] * u[1] + squarelist[i][kn][3][1] * v[1],
                                ulc[2] + squarelist[i][kn][3][0] * u[2] + squarelist[i][kn][3][1] * v[2] ]
                    
                    vert_str = "{:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}".format(
                               vertex1[0], vertex1[1], vertex1[2],
                               vertex2[0], vertex2[1], vertex2[2],
                               vertex3[0], vertex3[1], vertex3[2],
                               vertex4[0], vertex4[1], vertex4[2])
                    countervert_str = "{:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}, {:f}".format(
                                      vertex1[0], vertex1[1], vertex1[2],
                                      vertex4[0], vertex4[1], vertex4[2],
                                      vertex3[0], vertex3[1], vertex3[2],
                                      vertex2[0], vertex2[1], vertex2[2])
                    
                    # fill out the entries in the opening string definition
                    single_opening_def = opening_str.replace('<IDX>', str(i)+'_'+str(pn)).replace('<BASESURFACE>', base_surface).replace('<VERTICES>', vert_str).replace('<COUNTERVERTICES>', countervert_str)
                    idf_total_opening_def = idf_total_opening_def + single_opening_def

                    # a new consecutive subsequence has just begun provided that k<len(insideindices[i])
                    if k<len(insideindices[i]):
                        pn = insideindices[i][k]
                        kn = insideindices[i][k]
     
    # create idf opening objects from the string containing opening definitions   
    from io import StringIO
    idf_opening = IDF(StringIO(idf_total_opening_def))

    # copy idf shading objects to the existing idf file
    openings = idf_opening.idfobjects["FENESTRATIONSURFACE:DETAILED"]
    for opening in openings:
        idf_collection.copyidfobject(opening)
    
    ###############################
    # at the end, save the changes
    ###############################
    idf_collection.save()
