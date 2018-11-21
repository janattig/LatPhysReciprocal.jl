################################################################################
#
#	ABSTRACT TYPE
#
#   ReciprocalPath{P}
#   --> P is the reciprocal point type (<: AbstractReciprocalPoint{D})
#       --> D is the dimension of embedding space
#
#   FILE CONTAINS
#       - abstract type definition
#       - interface definition
#       - TODO interface testing
#
################################################################################


################################################################################
#
#   ABSTRACT TYPE DEFINITION
#
################################################################################
abstract type AbstractReciprocalPath{P <: AbstractReciprocalPoint{D} where {D}} end

# export the type
export AbstractReciprocalPath






################################################################################
#
#	INTERFACING / ACCESSING RECIPROCAL PATHS
#	(functions have to be overwritten by concrete types)
#
################################################################################


# default constructor interface
# used for creation of new reciprocal paths
function newReciprocalPath(
            ::Type{RPA},
            points :: Vector{<:RPO}
        ) :: RPA where {D,RPO<:AbstractReciprocalPoint{D}, RPA<:AbstractReciprocalPath{RPO}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'newReciprocalPath' for concrete reciprocal path type " *
            string(RPA) * " with reciprocal point type " * string(RPO) *
            " in " * string(D) * " spatial dimensions" )
end

# export the newUnitcell interface
export newReciprocalPath






# access the list of points
function points(
            p :: RPA
        ) :: Vector{RPO} where {D,RPO<:AbstractReciprocalPoint{D}, RPA<:AbstractReciprocalPath{RPO}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'points' for concrete reciprocal path type " * string(RPA) )
end
# set the list of points
function points!(
            p :: RPA,
            points :: Vector{<:RPO}
        ) where {D,RPO<:AbstractReciprocalPoint{D}, RPA<:AbstractReciprocalPath{RPO}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'points!' for concrete reciprocal path type " * string(RPA) )
end


# export access to points
export points, points!




# some more beauty interface
# builds on interface defined above but can also be overwritten

# number of reciprocal points
function numPoints(
            p :: RPA
        ) :: Int64 where {D,RPO<:AbstractReciprocalPoint{D}, RPA<:AbstractReciprocalPath{RPO}}

    # return length of the point array that the reciprocal path type implements
    return length(points(p))
end

# export the function
export numPoints




# access a specific reciprocal point
function point(
            path  :: RPA,
            index :: Int64
        ) :: RPO where {D,RPO<:AbstractReciprocalPoint{D}, RPA<:AbstractReciprocalPath{RPO}}

    # return the respective reciprocal point
    return points(path)[index]
end

# export the function
export point





# splatted constructor interface
# used for creation of new reciprocal paths based on number of points
# (does not have to be overwritten)
function newReciprocalPath(
            ::Type{RPA},
            points :: RPO ...
        ) :: RPA where {D,RPO<:AbstractReciprocalPoint{D}, RPA<:AbstractReciprocalPath{RPO}}

    # return the usual constructor
    return newReciprocalPath(RPA, RPO[p for p in points])
end
# set the list of points with splatted input
function points!(
            p :: RPA,
            points :: RPO ...
        ) where {D,RPO<:AbstractReciprocalPoint{D}, RPA<:AbstractReciprocalPath{RPO}}

    # return the fitting function
    return points!(p, RPO[p for p in points])
end
