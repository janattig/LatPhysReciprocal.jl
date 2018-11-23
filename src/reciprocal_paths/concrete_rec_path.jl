################################################################################
#
#	CONCRETE TYPE
#
#   ReciprocalPath{P} <: AbstractReciprocalPath{P}
#   --> P is the reciprocal point type (<: AbstractReciprocalPoint{D})
#       --> D is the dimension of embedding space
#
#   FILE CONTAINS
#       - concrete struct definition
#       - interface implementation
#
################################################################################



################################################################################
#
#   STRUCT DEFINITION
#
################################################################################
mutable struct ReciprocalPath{P} <: AbstractReciprocalPath{P}

    # the only field is the array of points
    points :: Vector{P}

end

# export the type
export ReciprocalPath






################################################################################
#
#	INTERFACING / ACCESSING RECIPROCAL PATHS
#
################################################################################


# default constructor interface
# used for creation of new reciprocal paths
function newReciprocalPath(
            ::Type{ReciprocalPath{RPO}},
            points :: Vector{<:RPO}
        ) :: ReciprocalPath{RPO} where {D,RPO<:AbstractReciprocalPoint{D}}

    # return the newly created object
    return ReciprocalPath{RPO}(points)
end

# export the newUnitcell interface
export newReciprocalPath






# access the list of points
function points(
            p :: ReciprocalPath{RPO}
        ) :: Vector{RPO} where {D,RPO<:AbstractReciprocalPoint{D}}

    # return the respective field
    return p.points
end
# set the list of points
function points!(
            p :: ReciprocalPath{RPO},
            points :: Vector{<:RPO}
        ) where {D,RPO<:AbstractReciprocalPoint{D}}

    # set the respective field
    p.points = points
end


# export access to points
export points, points!
