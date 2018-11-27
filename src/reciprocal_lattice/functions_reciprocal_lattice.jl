################################################################################
#
#   COLLECTION OF FUNCTIONS FOR
#   RECIPROCAL LATTICES / UNITCELLS
#
#   FILE CONTAINS
#   - construction of reciprocal unitcells from real space unitcells
#   - testing of reciprocality between unitcell and reciprocal unitcell
#   - reciprocal point information, shifting to 1st BZ
#
################################################################################








################################################################################
#
#   RECIPROCAL LATTICE / UNITCELL CONSTRUCTION
#
################################################################################

############
# FALLBACK #
############

# N=?,D=?
function getReciprocalUnitcell(
        ::Type{R},
        unitcell :: U
    ) :: R where {D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B},
                  LBR,BR<:AbstractBond{LBR,N},P<:AbstractReciprocalPoint{D},R<:AbstractReciprocalUnitcell{P,BR}}

    # give an error
    error("Cannot build a reciprocal unitcell for realspace lattice dimensions D=" * string(D) * " and N=" * string(N) * ".")
end

###########
# WRAPPER #
###########

# N=?,D=?
function getReciprocalUnitcell(
        unitcell :: U
    ) :: ReciprocalUnitcell{ReciprocalPoint{D},Bond{Int64,N}} where {D,N,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B},
                  LBR,BR<:AbstractBond{LBR,N},P<:AbstractReciprocalPoint{D},R<:AbstractReciprocalUnitcell{P,BR}}

    # pass on to general interface
    return getReciprocalUnitcell(ReciprocalUnitcell{ReciprocalPoint{D},Bond{Int64,N}}, unitcell)
end




###########
#   1D    #
###########

# N=1,D=?
function getReciprocalUnitcell(
        ::Type{R},
        unitcell :: U
    ) :: R where {D,LS,LB,S<:AbstractSite{LS,D},B<:AbstractBond{LB,1},U<:AbstractUnitcell{S,B},
                  LBR,BR<:AbstractBond{LBR,1},P<:AbstractReciprocalPoint{D},R<:AbstractReciprocalUnitcell{P,BR}}

    # build lattice vector
    latvec_b1 = a1(unitcell) .* 2*pi/sqrt(dot(a1(unitcell), a1(unitcell)))

    # return a new reciprocal unitcell
    return newReciprocalUnitcell(
        # pass the given type
        R,
        # pass the single lattice vector
        Vector{Float64}[
            latvec_b1,
        ],
        # pass the two relevant bonds
        BR[
            newBond(BR, 1,1, getDefaultLabel(LBR), (+1,)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (-1,))
        ]
    )
end



###########
#   2D    #
###########

# N=2, D=2
function getReciprocalUnitcell(
        ::Type{R},
        unitcell :: U
    ) :: R where {LS,LB,S<:AbstractSite{LS,2},B<:AbstractBond{LB,2},U<:AbstractUnitcell{S,B},
                  LBR,BR<:AbstractBond{LBR,2},P<:AbstractReciprocalPoint{2},R<:AbstractReciprocalUnitcell{P,BR}}

    # build reciprocal lattice vectors
    latvec_b1 = [a2(unitcell)[2], -a2(unitcell)[1]]
    latvec_b2 = [a1(unitcell)[2], -a1(unitcell)[1]]

    # normalize the reciprocal vectors
    latvec_b1 .*= 2*pi/dot(latvec_b1, a1(unitcell))
    latvec_b2 .*= 2*pi/dot(latvec_b2, a2(unitcell))

    # return a new reciprocal unitcell
    return newReciprocalUnitcell(
        # pass the given type
        R,
        # pass the single lattice vector
        Vector{Float64}[
            latvec_b1,
            latvec_b2
        ],
        # pass the two relevant bonds
        BR[
            newBond(BR, 1,1, getDefaultLabel(LBR), ( 0,+1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), ( 0,-1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (+1, 0)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (-1, 0)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (+1,+1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (-1,+1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (+1,-1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (-1,-1))
        ]
    )
end

# N=2, D=3
function getReciprocalUnitcell(
        ::Type{R},
        unitcell :: U
    ) :: R where {LS,LB,S<:AbstractSite{LS,3},B<:AbstractBond{LB,2},U<:AbstractUnitcell{S,B},
                  LBR,BR<:AbstractBond{LBR,2},P<:AbstractReciprocalPoint{3},R<:AbstractReciprocalUnitcell{P,BR}}

    # build a orthogonal vector
    ortho = cross(a1(unitcell), a2(unitcell))

    # build reciprocal lattice vectors
    latvec_b1 = cross(a2(unitcell), ortho)
    latvec_b2 = cross(a1(unitcell), ortho)

    # normalize the reciprocal vectors
    latvec_b1 .*= 2*pi/dot(latvec_b1, a1(unitcell))
    latvec_b2 .*= 2*pi/dot(latvec_b2, a2(unitcell))

    # return a new reciprocal unitcell
    return newReciprocalUnitcell(
        # pass the given type
        R,
        # pass the two lattice vectors
        Vector{Float64}[
            latvec_b1,
            latvec_b2
        ],
        # pass the 8 (3^2 - 1) relevant bonds
        BR[
            newBond(BR, 1,1, getDefaultLabel(LBR), ( 0,+1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), ( 0,-1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (+1, 0)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (-1, 0)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (+1,+1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (-1,+1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (+1,-1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (-1,-1))
        ]
    )
end

# N=2, D<2 not possible
# N=2, D>3 not implemented



###########
#   3D    #
###########

# N=2, D=3
function getReciprocalUnitcell(
        ::Type{R},
        unitcell :: U
    ) :: R where {LS,LB,S<:AbstractSite{LS,3},B<:AbstractBond{LB,3},U<:AbstractUnitcell{S,B},
                  LBR,BR<:AbstractBond{LBR,3},P<:AbstractReciprocalPoint{3},R<:AbstractReciprocalUnitcell{P,BR}}

    # construct reciprocal lattice vectors
    latvec_b1 = cross(a2(unitcell), a3(unitcell))
    latvec_b2 = cross(a3(unitcell), a1(unitcell))
    latvec_b3 = cross(a1(unitcell), a2(unitcell))
    # normalize the vectors
    latvec_b1 .*= 2*pi/sum(latvec_b1.*a1(unitcell))
    latvec_b2 .*= 2*pi/sum(latvec_b2.*a2(unitcell))
    latvec_b3 .*= 2*pi/sum(latvec_b3.*a3(unitcell))

    # return a new reciprocal unitcell
    return newReciprocalUnitcell(
        # pass the given type
        R,
        # pass the three lattice vectors
        Vector{Float64}[
            latvec_b1,
            latvec_b2,
            latvec_b3
        ],
        # pass the 26 (3^3 - 1) relevant bonds
        BR[
            newBond(BR, 1,1, getDefaultLabel(LBR), ( 0, 0,+1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), ( 0, 0,-1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), ( 0,+1, 0)),
            newBond(BR, 1,1, getDefaultLabel(LBR), ( 0,-1, 0)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (-1, 0, 0)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (+1, 0, 0)),
            newBond(BR, 1,1, getDefaultLabel(LBR), ( 0,+1,+1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), ( 0,-1,+1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), ( 0,+1,-1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), ( 0,-1,-1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (+1, 0,+1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (+1, 0,-1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (-1, 0,+1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (-1, 0,-1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (+1,+1, 0)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (+1,-1, 0)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (-1,+1, 0)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (-1,-1, 0)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (+1,+1,+1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (+1,-1,+1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (+1,+1,-1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (+1,-1,-1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (-1,+1,+1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (-1,-1,+1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (-1,+1,-1)),
            newBond(BR, 1,1, getDefaultLabel(LBR), (-1,-1,-1))
        ]
    )
end

# N=3, D<3 not possible
# N=3, D>3 not implemented




# export the function
export getReciprocalUnitcell










################################################################################
#
#   RECIPROCAL LATTICE / UNITCELL TESTING
#
################################################################################

# Test function
function testReciprocality(
        uc_real :: U,
        uc_reci :: R,
        error_allowed :: Float64 = 1e-8
    ) :: Bool where {N,UL,US,UB<:AbstractBond{UL,N},U<:AbstractUnitcell{US,UB},
                     RL,RS,RB<:AbstractBond{RL,N},R<:AbstractReciprocalUnitcell{RS,RB}}

    # the complete overlap
    overlap_sum = 0.0

    # test the overlap
    for l1 in 1:length(latticeVectors(uc_real))
    for l2 in 1:length(latticeVectors(uc_reci))
        # the value that it should take
        v_correct = l1==l2 ? 2*pi : 0.0
        # calculate the value
        v_taken = dot(latticeVectors(uc_real)[l1],latticeVectors(uc_reci)[l2])
        # add to the overlap sum
        overlap_sum += abs(v_correct - v_taken)
    end
    end

    # check if the overall overlap is sufficiently small
    return overlap_sum < error_allowed ? true : false
end

# export the function
export testReciprocality











################################################################################
#
#   RECIPROCAL POINT INFORMATION / SHIFTING TO 1st BZ
#
################################################################################


# TESTING IF IN 1st BZ
function isInFirstBZ(
            reci_unitcell :: RU,
            k             :: Vector{<:Real}
        ) :: Bool where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B}}

    # build up the minimal distance to any neighbor
    min_distance_squared = 1e40

    # iterate over all neighbors and find distance
    for b in bonds(reci_unitcell)
        # find out the position of the neighbor
        pos_neighbor = Float64[
            sum([latticeVectors(reci_unitcell)[j][i] * wrap(b)[j] for j in 1:N])  for i in 1:D
        ]
        # add the k point given
        pos_neighbor .+= k
        # check the minimum squared distance
        min_distance_squared = min(min_distance_squared, dot(pos_neighbor, pos_neighbor))
    end

    # return if the minimal distance to the neighbors is smaller than the distance to the origin (@[0,0,0,..])
    return (min_distance_squared >= dot(k,k))

end

# export the test
export isInFirstBZ




# SHIFTING TO 1st BZ
function shiftToFirstBZ(
            reci_unitcell :: RU,
            k             :: Vector{<:Real}
        ) :: Vector{Float64} where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B}}

    # return the inplace function with a copied input
    return shiftToFirstBZ!(reci_unitcell, deepcopy(k))
end

# SHIFTING TO 1st BZ (inplace)
function shiftToFirstBZ!(
            reci_unitcell :: RU,
            k             :: Vector{<:Real}
        ) :: Vector{Float64} where {D,N,L,P<:AbstractReciprocalPoint{D},B<:AbstractBond{L,N},RU<:AbstractReciprocalUnitcell{P,B}}

    # find out the current distance to gamma point
    current_distance = dot(k,k)

    # create a shifted k vector
    k_shifted = k

    # boolean indicating if the point was relocated in the last iteration
    relocated = true
    # try to shift as long as there is a possiblity to decrease the distance
    while relocated
        # assume that the point is not relocated in this iteration
        relocated = false
        # iterate over all neighbors and try to shift the point
        for b in bonds(reci_unitcell)
            # get the offset from the neighbor wrap
            offset = Float64[
                sum([latticeVectors(reci_unitcell)[j][i] * wrap(b)[j] for j in 1:N])  for i in 1:D
            ]
            # see if the distance is lowered and maybe shift
            if dot(k_shifted .+ offset, k_shifted .+ offset) < current_distance
                # shift the point
                k_shifted .+= offset
                # set the current distance correctly
                current_distance = dot(k_shifted,k_shifted)
                # denote as relocated
                relocated = true
            end
        end
    end

    # return the final shifted lattice vector
    return k_shifted
end

# export the test
export shiftToFirstBZ, shiftToFirstBZ!
