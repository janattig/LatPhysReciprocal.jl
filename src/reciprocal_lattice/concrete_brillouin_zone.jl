################################################################################
#
#	CONCRETE TYPE
#
#   ReciprocalUnitcell{P,B} <: AbstractReciprocalUnitcell{P,B}
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
mutable struct ReciprocalUnitcell{P,B} <: AbstractReciprocalUnitcell{P,B}

    # basis vectors of the Bravais lattice (of the reciprocal unitcell)
    lattice_vectors	:: Vector{Vector{Float64}}

    # the gamma point (only site within the reciprocal unitcell)
    gamma_point     :: P

    # list of bonds
    bonds			:: Vector{B}

end

# export the type
export ReciprocalUnitcell








################################################################################
#
#	INTERFACING / ACCESSING RECIPROCAL UNITCELLS
#	(functions that have to be overwritten by concrete types)
#
################################################################################

# constructor
function newReciprocalUnitcell(
            ::Type{ReciprocalUnitcell{P,B}},
            lattice_vectors :: Vector{<:Vector{<:Real}},
            bonds           :: Vector{B}
        ) :: ReciprocalUnitcell{P,B} where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N}}

    # call the constructor
    return ReciprocalUnitcell{P,B}(
        lattice_vectors,
        newReciprocalPoint(P, zeros(D), "Gamma", "\\Gamma"),
        bonds
    )
end

# export the constructor
export newReciprocalUnitcell




# accessing a list of lattice vectors
function latticeVectors(
            unitcell :: ReciprocalUnitcell{P,B}
        ) :: Vector{Vector{Float64}} where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N}}

    # return the lattice vectors
    return unitcell.lattice_vectors
end
# setting a list of lattice vectors
function latticeVectors!(
            unitcell        :: ReciprocalUnitcell{P,B},
            lattice_vectors :: Vector{<:Vector{<:Real}}
        ) where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N}}

    # set the lattice vectors
    unitcell.lattice_vectors = lattice_vectors
end

# export the latticeVectors interface
export latticeVectors, latticeVectors!




# accessing the gamma point
function gammaPoint(
            unitcell :: ReciprocalUnitcell{P,B}
        ) :: P where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N}}

    # return the gamma point
    return unitcell.gamma_point
end

# setting the gamma point
function gammaPoint!(
            unitcell :: ReciprocalUnitcell{P,B},
            site     :: P
        ) where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N}}

    # set the gamma point
    unitcell.gamma_point = site
end

# export gamma point interface
export gammaPoint, gammaPoint!




# accessing a list of bonds
function bonds(
            unitcell :: ReciprocalUnitcell{P,B}
        ) :: Vector{B} where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N}}

    # return the internal list of bonds
    return unitcell.bonds
end
# setting a list of bonds
function bonds!(
            unitcell :: ReciprocalUnitcell{P,B},
            bonds    :: Vector{B}
        ) where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N}}

    # set the internal list of bonds
    unitcell.bonds = bonds
end

# export the bonds interface
export bonds, bonds!
