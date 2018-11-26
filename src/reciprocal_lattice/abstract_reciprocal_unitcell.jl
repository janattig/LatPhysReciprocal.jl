################################################################################
#
#	ABSTRACT TYPE
#
#   ReciprocalUnitcell{P,B} <: Unitcell{P,B}
#   --> P is the reciprocal site type (<: AbstractReciprocalPoint{D})
#       --> D is the dimension of embedding space
#   --> B is the bond type (<: AbstractBond{L,N})
#       --> N is the dimension of the Bravais lattice that the bond is located in
#       --> L is the label type of bonds
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
abstract type AbstractReciprocalUnitcell{
        P <: AbstractReciprocalPoint{D} where {D},
        B <: AbstractBond{L,N} where {L,N}
    } <: AbstractUnitcell{P,B} end

# export the type
export AbstractReciprocalUnitcell






################################################################################
#
#	INTERFACING / ACCESSING UNITCELLS
#	(functions that are overwritten by abstract reciprocal type, NOT concrete)
#
################################################################################


# default constructor interface
# used for creation of new unitcells
function newUnitcell(
            ::Type{U},
            lattice_vectors :: Vector{<:Vector{<:Real}},
            sites           :: Vector{S},
            bonds           :: Vector{B}
        ) :: U where {D,N,L,S<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},U<:AbstractReciprocalUnitcell{S,B}}

    # check if only a single site is given
    @assert length(sites) == 1
    # check if the site is located at 0
    @assert sum(abs.(point(sites[1]))) < 1e-8

    # pass to the reciprocal unitcell constructor interface
    return newReciprocalUnitcell(U, lattice_vectors, bonds)
end

# export the newUnitcell interface
export newUnitcell



# accessing a list of sites
function sites(
            unitcell :: U
        ) :: Vector{S} where {D,N,L,S<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},U<:AbstractReciprocalUnitcell{S,B}}

    # return the Gamma point inside an array
    return S[gammaPoint(unitcell),]
end
# setting a list of sites
function sites!(
            unitcell :: U,
            sites    :: Vector{S}
        ) where {D,N,L,S<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},U<:AbstractReciprocalUnitcell{S,B}}

    # check if only a single site is given
    @assert length(sites) == 1
    # check if the site is located at 0
    @assert sum(abs.(point(sites[1]))) < 1e-8

    # set the Gamma point
    gammaPoint!(unitcell, sites[1])
end

# export the sites interface
export sites, sites!



################################################################################
#
#	INTERFACING / ACCESSING RECIPROCAL UNITCELLS
#	(functions that have to be overwritten by concrete types)
#
################################################################################

# constructor
function newReciprocalUnitcell(
            ::Type{R},
            lattice_vectors :: Vector{<:Vector{<:Real}},
            bonds           :: Vector{B}
        ) :: R where {D,N,L,P<:AbstractReciprocalPoint{D}, B<:AbstractBond{N,L}, R<:AbstractReciprocalUnitcell{P,B}}

    # throw an error because implementation for concrete type is missing
    error(  "not implemented function 'newReciprocalUnitcell' for concrete reciprocal unitcell type " *
            string(R) * " with bond type " * string(B) )
end

# export the constructor
export newReciprocalUnitcell




# accessing a list of lattice vectors
function latticeVectors(
            unitcell :: U
        ) :: Vector{Vector{Float64}} where {P,B,U<:AbstractReciprocalUnitcell{P,B}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'latticeVectors' for concrete reciprocal unitcell type " * string(U) )
end
# setting a list of lattice vectors
function latticeVectors!(
            unitcell        :: U,
            lattice_vectors :: Vector{<:Vector{<:Real}}
        ) where {P,B,U<:AbstractReciprocalUnitcell{P,B}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'latticeVectors!' for concrete reciprocal unitcell type " * string(U) )
end

# export the latticeVectors interface
export latticeVectors, latticeVectors!




# accessing the gamma point
function gammaPoint(
            unitcell :: U
        ) :: S where {D,N,L,S<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},U<:AbstractReciprocalUnitcell{S,B}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'gammaPoint' for concrete reciprocal unitcell type " * string(U) )
end

# setting the gamma point
function gammaPoint!(
            unitcell :: U,
            site     :: P
        ) :: S where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},U<:AbstractReciprocalUnitcell{P,B}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'gammaPoint!' for concrete reciprocal unitcell type " * string(U) )
end

# export gamma point interface
export gammaPoint, gammaPoint!




# accessing a list of bonds
function bonds(
            unitcell :: U
        ) :: Vector{B} where {N,L,P,B<:AbstractBond{L,N},U<:AbstractReciprocalUnitcell{P,B}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'bonds' for concrete reciprocal unitcell type " *
            string(U) * " with bond type " * string(B)   )
end
# setting a list of bonds
function bonds!(
            unitcell :: U,
            bonds    :: Vector{B}
        ) where {N,L,P,B<:AbstractBond{L,N},U<:AbstractReciprocalUnitcell{P,B}}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'bonds!' for concrete reciprocal unitcell type " *
            string(U) * " with bond type " * string(B)   )
end

# export the bonds interface
export bonds, bonds!
