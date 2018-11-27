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



# INTERFACE FUNCTIONS
#
# - corners(BZ), faces(BZ)
# - rec_unitcell(BZ)
#
