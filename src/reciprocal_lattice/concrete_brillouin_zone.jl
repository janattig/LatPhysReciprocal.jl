################################################################################
#
#	CONCRETE TYPE
#
#   BrillouinZone{R}
#   --> R is the reciprocal unitcell type (<: AbstractReciprocalUnitcell{P,B})
#       --> P is the reciprocal site type (<: AbstractReciprocalPoint{D})
#           --> D is the dimension of embedding space
#       --> B is the bond type (<: AbstractBond{L,N})
#           --> N is the dimension of the Bravais lattice that the bond is located in
#           --> L is the label type of bonds
#
#   FILE CONTAINS
#       - abstract type definition
#       - interface implementation
#
################################################################################


################################################################################
#
#   ABSTRACT TYPE DEFINITION
#
################################################################################
mutable struct BrillouinZone{R} <: AbstractBrillouinZone{R}

    # reciprocal unitcell
    reciprocal_unitcell :: R

    # corners
    corners :: Vector{Vector{Float64}}

    # faces
    faces :: Vector{Vector{Int64}}

end

# export the type
export BrillouinZone









################################################################################
#
#	INTERFACING / ACCESSING BRILLOUIN ZONES
#	(functions that have to be overwritten by concrete types)
#
################################################################################


# constructor
function newBrillouinZone(
            ::Type{BrillouinZone{RU}},
            reciprocal_unitcell :: RU,
            corners             :: Vector{<:Vector{<:Real}},
            faces               :: Vector{<:Vector{<:Integer}}
        ) :: BrillouinZone{RU} where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B}}

    # call the constructor
    return BrillouinZone{RU}(reciprocal_unitcell, corners, faces)
end

# export constructor
export newBrillouinZone




# accessing the corner points of the BZ
function corners(
            bz :: BrillouinZone{RU}
        ) :: Vector{Vector{Float64}} where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B}}

    # return the corners
    return bz.corners
end

# setting the corner points of the BZ
function corners!(
            bz      :: BrillouinZone{RU},
            corners :: Vector{<:Vector{<:Real}}
        ) where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B}}

    # set the corners
    bz.corners = corners
end

# export corner interface
export corners, corners!





# accessing the faces of the BZ
function faces(
            bz :: BrillouinZone{RU}
        ) :: Vector{Vector{Int64}} where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B}}

    # return the faces
    return bz.faces
end

# setting the faces of the BZ
function faces!(
            bz    :: BrillouinZone{RU},
            faces :: Vector{<:Vector{<:Integer}}
        ) where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B}}

    # set the faces
    bz.faces = faces
end

# export faces interface
export faces, faces!





# accessing the reciprocal unitcell of the BZ
function reciprocalUnitcell(
            bz :: BrillouinZone{RU}
        ) :: RU where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B}}

    # return the reciprocal unitcell
    return bz.reciprocal_unitcell
end

# setting the reciprocal unitcell of the BZ
function reciprocalUnitcell!(
            bz  :: BrillouinZone{RU},
            ruc :: RU
        ) where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B}}

    # set the reciprocal unitcell
    bz.reciprocal_unitcell = ruc
end

# export reciprocalUnitcell interface
export reciprocalUnitcell, reciprocalUnitcell!
