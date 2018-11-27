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







#=

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





################################################################################
#
#   COMPUTING BZ PARTS FROM RECIPROCAL UNITCELLS
#
################################################################################



#-------------
# RECIPROCAL LATTICE (from reciprocal unitcell)
#-------------

# general for all N,D
function computeBrillouinZoneConstructionLattice(
            reciprocal_unitcell :: RU
        ) :: Lattice{P,B,RU} where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B}}

    ##########
    # STEP 1 - Construct the reciprocal lattice points that are relevant
    ##########

    # list of reciprocal points
    k_points = Array{Float64, 1}[]

    # build list of points
    for i in -1:1
    for j in -1:1
    for l in -1:1
        # continue if gamma point
        if i==j==l==0
            continue
        end
        # add to the list
        push!(k_points, b1.*i .+ b2.*j .+ b3.*l)
    end
    end
    end

    # find the minimal distance to a the Gamma point
    max_distance = maximum([dot(k,k) for k in k_points])
    min_distance = minimum([dot(k,k) for k in k_points])

    # print
    if verbose
        println("maximum distance: $(max_distance)")
        println("minimum distance: $(min_distance)")
    end

    # clear list of reciprocal points
    k_points = Array{Float64, 1}[]

    # get a finite patch of reciprocal lattice
    lattice = getLatticeInSphere(rec_unit, sqrt(max_distance)*max_r)

    # get the points on the reciprocal lattice
    k_points = lattice.positions

    # print how many points
    if verbose
        println("$(length(k_points)) points of reciprocal lattice considered")
    end


end




#-------------
# CORNERS (from reciprocal lattice)
#-------------

# N=?, D=? (FALLBACK)
function computeBrillouinZoneCorners(
            reciprocal_lattice :: RL
        ) :: Vector{Vector{Float64}} where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B},LS,LB,RL<:AbstractLattice{LS,LB,RU}}

    # throw an error
    error("Currently no implemented function 'computeBrillouinZoneCorners' for reciprocal unitcell type " * string(RU) *
    " within reciprocal lattice of type " * string(RL) *
    "\nwhich explicitly means D=" * string(D) * " and N=" * string(N) * " of the reciprocal lattice")
end

# N=3, D=3
function computeBrillouinZoneCorners(
            reciprocal_lattice :: RL
        ) :: Vector{Vector{Float64}} where {L,P<:AbstractReciprocalPoint{3},B<:AbstractBond{L,3},RU<:AbstractReciprocalUnitcell{P,B},
                                            LS<:AbstractSite{LSL,3} where {LSL},LB,RL<:AbstractLattice{LS,LB,RU}}


    ##########
    # STEP 1 - Construct all mid points of all lattice sites as well as directions of normals
    ##########

    # list of mid points (achors of planes)
    mid_points = Vector{Float64}[
        k.*0.5  for k in point.(sites(reciprocal_lattice)) if dot(k,k) > 1e-5
    ]

    # list of normals of planes
    normals = Vector{Float64}[m ./ dot(m,m) for m in mid_points]



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

    # print how many points
    if verbose
        println("$(length(points)) corners of Brillouin zone found")
        println("$(length(planes_p)) planes of Brillouin zone found")
    end
end




#-------------
# FACES (from corners)
#-------------

# N=?, D=? (FALLBACK)
function computeBrillouinZoneFaces()

    ##########
    # STEP 4 - Connect the points in planes of the BZ --> edges/faces of BZ
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

        # print the points within the plane
        if verbose
            print("- points in plane $(i): $(point_indices)  ---  Loop: 1")
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
            # print the closest point
            if verbose
                print("- $(closest_point)")
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

        # close the line
        if verbose
            println("")
        end

        # get back the real indices and reconstruct loop
        loop = point_indices[loop]
        # get the sorted
        loop_sorted = sort(loop)
        # push the starting point into the loop
        push!(loop, loop[1])

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

    # print how many edges
    if verbose
        println("$(length(edges)) edges of Brillouin zone found")
    end
end







################################################################################
#
#   COMPUTING BRILLOUIN ZONE OBJECTS
#
################################################################################

# obtain the Brillouin zone for a given reciprocal unitcell
function getBrillouinZone(
            ::Type{BZ},
            reciprocal_unitcell :: RU
        ) :: BZ where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B},BZ<:AbstractBrillouinZone{RU}}


    ##########
    # STEP 1 - Compute the reciprocal lattice points that are needed for the calculation
    ##########

    # call the respective function
    reciprocal_lattice = computeBrillouinZoneConstructionLattice(reciprocal_unitcell)



    ##########
    # STEP 2 - Compute the corners of the BZ
    ##########

    # call the respective function
    corners_BZ = computeBrillouinZoneCorners(reciprocal_lattice)



    ##########
    # STEP 3 - Compute the faces of the BZ (based on corners)
    ##########

    # call the respective function
    faces_BZ = computeBrillouinZoneFaces(corners_BZ)



    ##########
    # STEP 4 - Construct a new object and return it
    ##########

    # return the finished BZ
    return newBrillouinZone(
        BZ,
        reciprocal_unitcell
        corners_BZ,
        faces_BZ
    )
end
=#
