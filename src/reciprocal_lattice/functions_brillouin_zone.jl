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
    reciprocal_lattice = getLatticeByBondDistance(reciprocal_unitcell, 2)

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
    push!(loop, loop[1])

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


################################################################################
#
#   COMPUTING BRILLOUIN ZONE OBJECTS
#   (Wrapper)
#
################################################################################

# wrapper for 2d unitcells
function getBrillouinZone(
            reciprocal_unitcell :: RU
            ;
            kwargs...
        ) :: BrillouinZone{RU} where {L,N,D,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B}}

    # call the respective function
    return getBrillouinZone(BrillouinZone{RU}, reciprocal_unitcell; kwargs...)
end

# Export the function
export getBrillouinZone
