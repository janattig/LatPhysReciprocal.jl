################################################################################
#
#	ABSTRACT TYPE
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
#       - interface definition
#       - TODO interface testing
#
################################################################################

################################################################################
#
#   ABSTRACT TYPE DEFINITION
#
################################################################################
abstract type AbstractBrillouinZone{
        R <: AbstractReciprocalUnitcell{P,B} where {
            P <: AbstractReciprocalPoint{D} where {D},
            B <: AbstractBond{L,N} where {L,N},
        }
    } end

# export the type
export AbstractBrillouinZone









################################################################################
#
#	INTERFACING / ACCESSING BRILLOUIN ZONES
#	(functions that have to be overwritten by concrete types)
#
################################################################################


# constructor
function newBrillouinZone(
            ::Type{BZ},
            reciprocal_unitcell :: RU,
            corners             :: Vector{<:Vector{<:Real}},
            faces               :: Vector{<:Vector{<:Integer}}
        ) :: BZ where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B},BZ<:AbstractBrillouinZone{RU}}

    # throw an error because implementation for concrete type is missing
    error(  "not implemented function 'newBrillouinZone' for concrete Brillouin zone type " *
            string(BZ) * " with reciprocal unitcell type " * string(RU) )
end

# export constructor
export newBrillouinZone




# accessing the corner points of the BZ
function corners(
            bz :: BZ
        ) :: Vector{Vector{Float64}} where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B},BZ<:AbstractBrillouinZone{RU}}

    # throw an error because implementation for concrete type is missing
    error(  "not implemented function 'corners' for concrete Brillouin zone type " *
            string(BZ) * " with reciprocal unitcell type " * string(RU) )
end

# setting the corner points of the BZ
function corners!(
            bz      :: BZ,
            corners :: Vector{<:Vector{<:Real}}
        ) where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B},BZ<:AbstractBrillouinZone{RU}}

    # throw an error because implementation for concrete type is missing
    error(  "not implemented function 'corners!' for concrete Brillouin zone type " *
            string(BZ) * " with reciprocal unitcell type " * string(RU) )
end

# export corner interface
export corners, corners!





# accessing the faces of the BZ
function faces(
            bz :: BZ
        ) :: Vector{Vector{Int64}} where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B},BZ<:AbstractBrillouinZone{RU}}

    # throw an error because implementation for concrete type is missing
    error(  "not implemented function 'faces' for concrete Brillouin zone type " *
            string(BZ) * " with reciprocal unitcell type " * string(RU) )
end

# setting the faces of the BZ
function faces!(
            bz    :: BZ,
            faces :: Vector{<:Vector{<:Integer}}
        ) where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B},BZ<:AbstractBrillouinZone{RU}}

    # throw an error because implementation for concrete type is missing
    error(  "not implemented function 'faces!' for concrete Brillouin zone type " *
            string(BZ) * " with reciprocal unitcell type " * string(RU) )
end

# export faces interface
export faces, faces!





# accessing the reciprocal unitcell of the BZ
function reciprocalUnitcell(
            bz :: BZ
        ) :: RU where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B},BZ<:AbstractBrillouinZone{RU}}

    # throw an error because implementation for concrete type is missing
    error(  "not implemented function 'reciprocalUnitcell' for concrete Brillouin zone type " *
            string(BZ) * " with reciprocal unitcell type " * string(RU) )
end

# setting the reciprocal unitcell of the BZ
function reciprocalUnitcell!(
            bz  :: BZ,
            ruc :: RU
        ) where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B},BZ<:AbstractBrillouinZone{RU}}

    # throw an error because implementation for concrete type is missing
    error(  "not implemented function 'reciprocalUnitcell!' for concrete Brillouin zone type " *
            string(BZ) * " with reciprocal unitcell type " * string(RU) )
end

# export reciprocalUnitcell interface
export reciprocalUnitcell, reciprocalUnitcell!




# CUSTOM SHOW FUNCTION

# single Brillouin Zone
function show(io::IO, bz::BZ) where {P,B,RU<:AbstractReciprocalUnitcell{P,B},BZ<:AbstractBrillouinZone{RU}}
    print(io, "Brillouin Zone object\n--> reciprocal unitcell type ", RU, "\n--> type ", BZ, "\n--> ", length(corners(bz)), " BZ corners\n--> ", length(faces(bz)), " BZ faces")
end
