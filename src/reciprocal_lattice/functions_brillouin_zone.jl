################################################################################
#
#   COLLECTION OF FUNCTIONS FOR
#   BRILLOUIN ZONES
#
#   FILE CONTAINS
#   - mathematical helper functions for planes, lines and points
#   - construction of BZ parts (corners, faces, etc.)
#   - construction of BZ objects from reciprocal unitcells
#
################################################################################








################################################################################
#
#   MATHEMATICAL HELPER FUNCTIONS
#
################################################################################

# function to find out if two vectors are parallel
function parallelVectors(
            v1 :: Vector{<:Real},
            v2 :: Vector{<:Real}
        ) :: Bool

    # check if they can be scaled into one another (i.e. the rank of the matrix)
    r = rank([v1'; v2'])
    if r == 1
        return true
    elseif r == 2
        return false
    else
        error("Got " * string(r) * " from calculating the rank of two vectors:\n" *
        "v1 = " * string(v1) * ", and v2 = " * string(v2))
    end
end

# Find out if a point point_p is on a plane given by (plane_p,plane_n)
function pointOnPlane(
            point_p :: Vector{<:Real},
            plane_p :: Vector{<:Real},
            plane_n :: Vector{<:Real},
            error_allowed :: Real = 1e-15
        ) :: Bool

    # check if the plane normal is perpendicular to the connecting vector between plane basis point and point
    return dot(plane_p .- point_p, plane_n) < error_allowed
end

# Find out if two planes (p1,n2) and (p2,n2) are identical
function identicalPlanes(
            p1 :: Vector{<:Real},
            n1 :: Vector{<:Real},
            p2 :: Vector{<:Real},
            n2 :: Vector{<:Real},
            error_allowed :: Real = 1e-15
        ) :: Bool

    # check if p2 is on plane 1 and the normals are parallel
    return pointOnPlane(p2, p1,n1, error_allowed) && parallelVectors(n1,n2)
end

# Find out if three planes (pi,ni) have one intersecting point
function haveSingleIntersection(
            p1 :: Vector{<:Real},
            n1 :: Vector{<:Real},
            p2 :: Vector{<:Real},
            n2 :: Vector{<:Real},
            p3 :: Vector{<:Real},
            n3 :: Vector{<:Real},
            error_allowed :: Real = 1e-15
        ) :: Bool

    # calculate the determinant
    determinant = dot(n3, cross(n1, n2))
    # check if they intersect in a point
    if abs(determinant) < error_allowed
        # there will not be a single point of intersection
        return false
    else
        # they will intersect in one point
        return true
    end
end

# Get the intersecting point of three planes
function getSingleIntersection(
            p1 :: Vector{<:Real},
            n1 :: Vector{<:Real},
            p2 :: Vector{<:Real},
            n2 :: Vector{<:Real},
            p3 :: Vector{<:Real},
            n3 :: Vector{<:Real},
            error_allowed :: Real = 1e-15
        ) :: Vector{Float64}

    # calculate the main determinant
    determinant = dot(n3, cross(n1, n2))

    # calculate the other relevant determinants
    D = Float64[-dot(n1,p1), -dot(n2,p2), -dot(n3,p3)]
    det_x = det([
        D[1] n1[2] n1[3]
        D[2] n2[2] n2[3]
        D[3] n3[2] n3[3]
    ])
    det_y = det([
        n1[1] D[1] n1[3]
        n2[1] D[2] n2[3]
        n3[1] D[3] n3[3]
    ])
    det_z = det([
        n1[1] n1[2] D[1]
        n2[1] n2[2] D[2]
        n3[1] n3[2] D[3]
    ])
    # get the interesection point
    return [det_x, det_y, det_z] ./ determinant
end






################################################################################
#
#   COMPUTING BRILLOUIN ZONE OBJECTS
#   (GENERAL FUNCTIONS)
#
################################################################################

# obtain the Brillouin zone for a given reciprocal unitcell (L=2,N=2)
function getBrillouinZone(
            ::Type{BZ},
            reciprocal_unitcell :: RU
            ;
            max_r :: Real = 1.001
        ) :: BZ where {L,P<:AbstractReciprocalPoint{2},B<:AbstractBond{L,2},RU<:AbstractReciprocalUnitcell{P,B},BZ<:AbstractBrillouinZone{RU}}


    ##########
    # STEP 1 - Compute the reciprocal lattice points that are needed for the calculation
    ##########

    # build a test lattice
    reciprocal_lattice = getLatticePeriodic(reciprocal_unitcell, 2)

    # find the minimal distance to a the Gamma point
    max_distance = maximum([dot(k,k) for k in point.(sites(reciprocal_lattice))])
    min_distance = minimum([dot(k,k) for k in point.(sites(reciprocal_lattice))])

    # reconstruct the reciprocal lattice
    reciprocal_lattice = getLatticeInSphere(reciprocal_unitcell, sqrt(max_distance)*max_r)

    # list of mid points
    mid_points = Vector{Float64}[
        k.*0.5  for k in point.(sites(reciprocal_lattice)) if sum(k.*k) > 1e-5
    ]

    # list of normals
    normals = Vector{Float64}[
        [k[2], -k[1]] ./ (k[1]*k[1] + k[2]*k[2]) for k in mid_points
    ]



    ##########
    # STEP 2 - Compute the corners of the BZ
    ##########


    ############
    # STEP 2.1 - Compute insersection points
    ############

    # list of intersection points
    intersections = Vector{Float64}[]

    # find the intersections
    for i in 1:length(mid_points)
    for j in 1:length(mid_points)
        # continue if same index
        if i==j
            continue
        end
        # get data for point 1
        p1 = mid_points[i]
        d1 = normals[i]
        # get data for point 2
        p2 = mid_points[j]
        d2 = normals[j]
        # check the interesection point
        delta_1 =  d1[2]*p1[1] - d1[1]*p1[2]
        delta_2 =  d2[2]*p2[1] - d2[1]*p2[2]
        delta   = -d1[2]*d2[1] + d2[2]*d1[1]
        # push to the intersections (if not devided by zero)
        if abs(delta) > 1e-8
            push!(intersections,[
                (d1[1]*delta_2 - d2[1]*delta_1) / delta,
                (d1[2]*delta_2 - d2[2]*delta_1) / delta
            ])
        end
    end
    end


    ############
    # STEP 2.2 - Find closest intersection points --> points of BZ
    ############

    # list of distances
    intersection_distances = Float64[
        sum(isec.*isec) for isec in intersections
    ]

    # find the minimal intersection distance
    min_distance = minimum(intersection_distances)

    # list of points which are only the minimal intersection distance away (=corners)
    corners_BZ = Array{Float64, 1}[]

    # check all points
    for i in 1:length(intersections)
        # check the relative deviation to the calculated minimal distance
        if intersection_distances[i] < min_distance * (1.0 + 1e-8)
            # check if not inserted already
            inserted_already = false
            for p in corners_BZ
                if sum((p.-intersections[i]).*(p.-intersections[i])) < 1e-8
                    inserted_already = true
                end
            end
            # found a new point
            if !inserted_already
                push!(corners_BZ, intersections[i])
            end
        end
    end



    ##########
    # STEP 3 - Compute the faces of the BZ (based on corners)
    ##########

    # list of points if they are contained in the loop
    in_loop = Bool[false for p in corners_BZ]

    # the list of the loop (to bring the points in correct order)
    loop = Int64[]

    # start with the first point
    push!(loop, 1)
    in_loop[1] = true

    # look for next point (which is closest to the given point in all remaing ones)
    while length(loop) < length(corners_BZ)
        # check for the closest
        closest_distance = sum(corners_BZ[1].*corners_BZ[1]) * 1000
        closest_point = -1
        # check all other points
        for i in 1:length(corners_BZ)
            # if already in loop, continue
            if in_loop[i]
                continue
            end
            # compare to the nearest one
            if sum((corners_BZ[i].-corners_BZ[loop[end]]).*(corners_BZ[i].-corners_BZ[loop[end]])) < closest_distance
                closest_point    = i
                closest_distance = sum((corners_BZ[i].-corners_BZ[loop[end]]).*(corners_BZ[i].-corners_BZ[loop[end]]))
            end
        end
        # insert the closest point as the next in line
        push!(loop, closest_point)
        in_loop[closest_point] = true
    end

    # push the starting point into the loop to wrap up and make a closed line
    #push!(loop, loop[1])

    # list of edges (only contains this one loop)
    faces_BZ = Vector{Int64}[loop]


    ##########
    # STEP 4 - Construct a new object and return it
    ##########

    # return the finished BZ
    return newBrillouinZone(
        BZ,
        reciprocal_unitcell,
        corners_BZ,
        faces_BZ
    )
end


# obtain the Brillouin zone for a given reciprocal unitcell (L=3,N=3)
function getBrillouinZone(
            ::Type{BZ},
            reciprocal_unitcell :: RU
            ;
            max_r :: Real = 1.001
        ) :: BZ where {L,P<:AbstractReciprocalPoint{3},B<:AbstractBond{L,3},RU<:AbstractReciprocalUnitcell{P,B},BZ<:AbstractBrillouinZone{RU}}


    ##########
    # STEP 1 - Compute the reciprocal lattice points that are needed for the calculation
    ##########

    # build a test lattice
    reciprocal_lattice = getLatticePeriodic(reciprocal_unitcell, 2)

    # find the minimal distance to a the Gamma point
    max_distance = maximum([dot(k,k) for k in point.(sites(reciprocal_lattice))])
    min_distance = minimum([dot(k,k) for k in point.(sites(reciprocal_lattice))])

    # reconstruct the reciprocal lattice
    reciprocal_lattice = getLatticeInSphere(reciprocal_unitcell, sqrt(max_distance)*max_r)


    # k points
    k_points = point.(sites(reciprocal_lattice))

    # list of mid points
    mid_points = Vector{Float64}[
        k.*0.5  for k in point.(sites(reciprocal_lattice)) if sum(k.*k) > 1e-5
    ]

    # list of normals of planes
    normals = [m ./ dot(m,m) for m in mid_points]




    ##########
    # STEP 2 - Find all intersection points of 3 planes and check that they are closest to the Gamma point
    ##########

    # list of points which are only the minimal intersection distance to Gamma point
    points = Array{Float64, 1}[]

    # temp point
    point_tmp = zeros(Float64, 3)


    # lists for all planes that aided in the construction
    planes_n = Array{Float64,1}[]
    planes_p = Array{Float64,1}[]



    # find the intersections
    for i in 1:length(mid_points)
        # get data for plane 1
        p1 = mid_points[i]
        n1 = normals[i]
        # get a second plane
        for j in i+1:length(mid_points)
            # get data for plane 2
            p2 = mid_points[j]
            n2 = normals[j]
            # continue if they are parallel
            if dot(cross(n1,n2), cross(n1,n2)) < 1e-8
                continue
            end
            # get a third plane
            for l in j+1:length(mid_points)
                # get data for plane 3
                p3 = mid_points[l]
                n3 = normals[l]
                # calculate the determinant
                determinant = dot(n3, cross(n1, n2))
                # check if they intersect in a point
                if abs(determinant) < 1e-8
                    # there will not be one point of intersection
                    continue
                end
                # calculate the other relevant determinants
                #A = Float64[n1[1], n2[1], n3[1]]
                #B = Float64[n1[2], n2[2], n3[2]]
                #C = Float64[n1[3], n2[3], n3[3]]
                D = Float64[-dot(n1,p1), -dot(n2,p2), -dot(n3,p3)]
                #det_x = dot(D, cross(B, C))
                #det_y = dot(A, cross(D, C))
                #det_z = dot(A, cross(B, D))
                det_x = det([
                    D[1] n1[2] n1[3]
                    D[2] n2[2] n2[3]
                    D[3] n3[2] n3[3]
                ])
                det_y = det([
                    n1[1] D[1] n1[3]
                    n2[1] D[2] n2[3]
                    n3[1] D[3] n3[3]
                ])
                det_z = det([
                    n1[1] n1[2] D[1]
                    n2[1] n2[2] D[2]
                    n3[1] n3[2] D[3]
                ])
                # get the interesection point
                point_tmp[1] = det_x / determinant
                point_tmp[2] = det_y / determinant
                point_tmp[3] = det_z / determinant
                # check if the distance to any other k point in the lattice is shorter than to the Gamma point
                shorter_distance_found = false
                distance_to_Gamma = dot(point_tmp, point_tmp)
                # check all points
                for k in k_points
                    if dot(point_tmp .- k, point_tmp .- k)*(1 + 1e-8) < distance_to_Gamma
                        shorter_distance_found = true
                        break
                    end
                end
                # insert point into list if no shorter distance found
                if !shorter_distance_found
                    # check if not inserted already
                    inserted_already = false
                    # check all points
                    for p in points
                        if dot(p.-point_tmp, p.-point_tmp) < 1e-6
                            inserted_already = true
                            break
                        end
                    end
                    # if not inserted already, insert now
                    if !inserted_already
                        # insert the point into the list
                        push!(points, copy(point_tmp))
                        # see if the planes that cross here, have to be inserted as well

                        # check plane 1
                        insert_p1 = true
                        # look through all planes that are inserted already
                        for index in 1:length(planes_p)
                            # check if the point lies in the plane
                            if abs(dot(p1 .- planes_p[index], planes_n[index])) < 1e-6
                                # point lies in this plane, check the direction of the normal
                                if dot(cross(planes_n[index], n1), cross(planes_n[index], n1)) < 1e-6
                                    # normals also match in direction
                                    insert_p1 = false
                                    break
                                end
                            end
                        end
                        # if not found, insert into lists
                        if insert_p1
                            push!(planes_p, p1)
                            push!(planes_n, n1)
                        end

                        # check plane 2
                        insert_p2 = true
                        # look through all planes that are inserted already
                        for index in 1:length(planes_p)
                            # check if the point lies in the plane
                            if abs(dot(p2 .- planes_p[index], planes_n[index])) < 1e-6
                                # point lies in this plane, check the direction of the normal
                                if dot(cross(planes_n[index], n2), cross(planes_n[index], n2)) < 1e-6
                                    # normals also match in direction
                                    insert_p2 = false
                                    break
                                end
                            end
                        end
                        # if not found, insert into lists
                        if insert_p2
                            push!(planes_p, p2)
                            push!(planes_n, n2)
                        end

                        # check plane 3
                        insert_p3 = true
                        # look through all planes that are inserted already
                        for index in 1:length(planes_p)
                            # check if the point lies in the plane
                            if abs(dot(p3 .- planes_p[index], planes_n[index])) < 1e-6
                                # point lies in this plane, check the direction of the normal
                                if dot(cross(planes_n[index], n3), cross(planes_n[index], n3)) < 1e-6
                                    # normals also match in direction
                                    insert_p3 = false
                                    break
                                end
                            end
                        end
                        # if not found, insert into lists
                        if insert_p3
                            push!(planes_p, p3)
                            push!(planes_n, n3)
                        end

                    end
                end
            end
        end
    end


    ##########
    # STEP 5 - Connect the points in planes of the BZ --> edges of BZ
    ##########

    # the list of all the loops (i.e. all edge loops of the BZ)
    edges        = Array{Int64,1}[]
    # the list of all the loops but sorted without double point for the start+end
    edges_sorted = Array{Int64,1}[]

    # find all loops of all points by constructing planes

    # iterate over all planes
    for i in 1:length(planes_n)

        # find the subset of points that lie within the plane
        in_plane = Bool[abs(dot(p .- planes_p[i], planes_n[i])) < 1e-6 for p in points]

        # use the 2D loop finding within the plane
        point_indices = collect(1:length(points))[in_plane]

        # check if there are enough points to form a loop
        if length(point_indices) < 3
            # not enough points found
            continue
        end

        # list of points if they are contained in the loop
        in_loop = [false for i in point_indices]

        # the list of the loop (to bring the points in correct order)
        # with relative indices (have to be translated back)
        loop = Int64[]

        # start with the first (relative) point
        push!(loop, 1)
        in_loop[1] = true

        # look for next point (which is closest to the given point in all remaing ones)
        while length(loop) < length(point_indices)
            # check for the closest
            closest_distance = 1e32
            closest_point = -1
            # check all other points
            for i in 1:length(point_indices)
                # if already in loop, continue
                if in_loop[i]
                    continue
                end
                # compare to the nearest one
                if dot(points[point_indices[i]].-points[point_indices[loop[end]]], points[point_indices[i]].-points[point_indices[loop[end]]]) < closest_distance
                    closest_point    = i
                    closest_distance = dot(points[point_indices[i]].-points[point_indices[loop[end]]], points[point_indices[i]].-points[point_indices[loop[end]]])
                end
            end
            # security check
            if closest_point == -1
                println()
                println(closest_distance)
                println(points[in_plane])
            end
            # insert the closest point as the next in line
            push!(loop, closest_point)
            in_loop[closest_point] = true
        end

        # get back the real indices and reconstruct loop
        loop = point_indices[loop]
        # get the sorted
        loop_sorted = sort(loop)
        # push the starting point into the loop
        #push!(loop, loop[1])

        # loop not found yet
        loop_found = false
        # check if the loop already exists
        for l in edges_sorted
            if l==loop_sorted
                loop_found = true
                break
            end
        end
        # if the loop is not found, insert the loop into the edges list
        if !loop_found
            push!(edges, loop)
        end

    end



    ##########
    # STEP 4 - Construct a new object and return it
    ##########

    # return the finished BZ
    return newBrillouinZone(
        BZ,
        reciprocal_unitcell,
        points,
        edges
    )
end







################################################################################
#
#   COMPUTING BRILLOUIN ZONE OBJECTS
#   (Wrapper)
#
################################################################################

# wrapper for all reciprocal unitcells
function getBrillouinZone(
            reciprocal_unitcell :: RU
            ;
            kwargs...
        ) :: BrillouinZone{RU} where {L,N,D,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B}}

    # call the respective function
    return getBrillouinZone(BrillouinZone{RU}, reciprocal_unitcell; kwargs...)
end

# wrapper for all real space unitcells
function getBrillouinZone(
            unitcell :: U
            ;
            kwargs...
        ) :: BrillouinZone{ReciprocalUnitcell{ReciprocalPoint{D},B}} where {L,D,S<:AbstractSite{L,D},B,U<:AbstractUnitcell{S,B}}

    # call the respective function
    return getBrillouinZone(BrillouinZone{ReciprocalUnitcell{ReciprocalPoint{D},B}}, getReciprocalUnitcell(unitcell); kwargs...)
end



# Export the function
export getBrillouinZone










################################################################################
#
#   FUNCTIONS FOR HIGH SYMMETRY MESH
#
################################################################################

# contructor for symmetric mesh of given brillouin zone
function getPointsHighSymmetryMesh(
            bz              :: BZ,
            line_resolution :: Integer
            ;
            include_bz_corners :: Bool           = true,
            include_face_centers :: Bool         = true,
            include_bz_edges :: Bool             = true,  # edges of faces (i.e. lines connecting the corners)
            include_lines_in_faces :: Bool       = true,
            include_lines_to_face_center :: Bool = true,
            include_lines_to_gamma :: Bool       = true,
            include_random_face_points :: Bool   = true,
            face_resolution :: Integer = 100
        ) :: Vector{Vector{Float64}} where {D,B, P<:AbstractReciprocalPoint{D}, R<:AbstractReciprocalUnitcell{P,B}, BZ<:AbstractBrillouinZone{R}}

    # get the gamma point
    gamma = point(gammaPoint(reciprocalUnitcell(bz)))
    # the list of points that form the mesh
    mesh = Vector{Float64}[]

    # iterate over all faces of the BZ
    for f in faces(bz)

        # obtain the relevant corners
        corners = unique([corners(bz)[i] for i in f])
        # points that connect to gamma
        points_to_gamma = Vector{Float64}[]

        # INCLUDE CORNER POINTS IN MESH
        if include_bz_corners
            for c in corners
                # put the corner into the mesh
                push!(mesh, c)
                # put the corner into the list of points that connect to the gamma point
                push!(points_to_gamma, c)
            end
        end

        # INCLUDE FACE CENTERS IN MESH
        if include_face_centers
            # push the center into the mesh
            push!(mesh, getPointCenter(corners))
            # push the center into the list of connecting points
            push!(points_to_gamma, getPointCenter(corners))
        end

        # INCLUDE EDGES IN MESH
        if include_bz_edges
            # the edge that wraps in the point list
            line = getPointsOnLine(corners[end], corners[1], line_resolution)
            for p in line
                # put the point into the mesh
                push!(mesh, p)
            end
            # put the start, end and center of line into the gamma point connecting point list
            push!(points_to_gamma, corners[end])
            push!(points_to_gamma, corners[1])
            push!(points_to_gamma, getPointCenter(corners[1], corners[end]))
            # all other edges
            for i in 1:length(corners)-1
                # construct line
                line = getPointsOnLine(corners[i], corners[i+1], line_resolution)
                for p in line
                    # put the point into the mesh
                    push!(mesh, p)
                end
                # put the start, end and center of line into the gamma point connecting point list
                push!(points_to_gamma, corners[i])
                push!(points_to_gamma, corners[i+1])
                push!(points_to_gamma, getPointCenter(corners[i], corners[i+1]))
            end
        end

        # INCLUDE LINES IN THE FACE
        if include_lines_in_faces
            # iterate over all combination of points
            for i in 1:length(corners)
            for j in i+1:length(corners)
                # build the line from i to j
                line = getPointsOnLine(corners[i], corners[j], line_resolution)
                for p in line
                    push!(mesh, p)
                end
                # put the start, end and center of line into the gamma point connecting point list
                push!(points_to_gamma, corners[i])
                push!(points_to_gamma, corners[j])
                push!(points_to_gamma, getPointCenter(corners[i], corners[j]))
            end
            end
        end

        # INCLUDE LINES TO THE FACE CENTER
        if include_lines_to_face_center
            # build the face center
            center = getPointCenter(corners)
            # iterate over all points
            for i in 1:length(corners)
                # build the line from i to j
                line = getPointsOnLine(corners[i], center, line_resolution)
                for p in line
                    push!(mesh, p)
                end
            end
        end

        # INCLUDE LINES CONNECTING TO GAMMA POINT
        if include_lines_to_gamma
            # get a unique list of connecting points
            unique!(points_to_gamma)
            # build lines from these points to the gamma point (if they are not the gamma point)
            for ptg in points_to_gamma
                if dot(ptg.-gamma, ptg.-gamma) < 1e-8
                    continue
                end
                line = getPointsOnLine(ptg, gamma, line_resolution)
                for p in line
                    push!(mesh, p)
                end
            end
        end

        # INCLUDE FACE POINTS RANDOMLY SAMPLED
        if include_random_face_points
            points = getPointsRandomInConvexHull(face_resolution, corners...)
            for p in points
                push!(mesh, p)
            end
        end

    end

    # remove redundant points
    unique!(mesh)

    # return the mesh
    return mesh
end

export getPointsHighSymmetryMesh
