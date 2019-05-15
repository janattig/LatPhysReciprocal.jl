################################################################################
#
#   FUNCTIONS FOR CONSTRUCTION OF K SPACE POINTS
#
#   FILE CONTAINS
#   - general building of points
#   - high symmetry mesh
#
#   ALL FUNCTIONS ARE NOT EXPORTED
#
################################################################################






################################################################################
#
#   FUNCTIONS FOR GENERAL CONSTRUCTION OF K SPACE POINTS
#
################################################################################

# contruct line between vertices with a certain resolution
function getPointsOnLine(
            point1 :: Vector{<:Real},
            point2 :: Vector{<:Real},
            resolution :: Integer
        ) :: Vector{Vector{Float64}}

    # the array of points which is later returned
    line = Vector{Vector{Float64}}(undef, resolution)
    # the multipliers between the two points
    alphas = range(0.0, stop=1.0, length=resolution)

    # set all line points
    for i in 1:resolution
        line[i] = (point1.*(1-alphas[i])) .+ (point2.*alphas[i])
    end

    # return the array of points
    return line
end
export getPointsOnLine

# get random points inside a box
function getPointsRandomInBox(
            N           :: Integer,
            center      :: Vector{<:Real},
            dimensions  :: Vector{<:Real}
        )

    # the list of points to return later
    points = Vector{Vector{Float64}}(undef, N)
    # the dimension of points
    d = length(dimensions)

    # set all points
    for i in 1:N
        # set the point to a random point
        points[i] = rand(Float64,d)
        # shift the point into the box
        points[i] .*= 2
        points[i] .-= 1
        points[i] .*= dimensions
        points[i] .+= center
    end

    # return the array of points
    return points
end
export getPointsRandomInBox

# get random points inside a convex hull of some points
function getPointsRandomInConvexHull(
            N             :: Integer,
            anchor_points :: Vector{<:Real}...
        )

    # the list of points to return later
    points = Vector{Vector{Float64}}(undef, N)

    # set all points
    for i in 1:N
        # set the point to a random point
        alphas = rand(Float64,length(anchor_points))
        alphas ./= sum(alphas)
        # shift the point into the box
        points[i] = sum(anchor_points .* alphas)
    end

    # return the array of points
    return points
end
export getPointsRandomInConvexHull

# compute the center as half the diagonal from a set of vertices forming a BZ surface
function getPointCenter(
            points :: Vector{<:Vector{<:Real}}
        ) :: Vector{Float64}

    # return sum of all points and divide by number of points
    return sum(points) ./ length(points)
end
function getPointCenter(
            points :: Vector{<:Real} ...
        ) :: Vector{Float64}

    # return sum of all points and divide by number of points
    return sum(points) ./ length(points)
end
export getPointCenter
