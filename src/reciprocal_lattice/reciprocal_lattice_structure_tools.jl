################################################################################
#
#   RECIPROCAL LATTICE / UNITCELL CONSTRUCTION / TESTING
#
################################################################################

# CONSTRUCT THE RECIPROCAL UNITCELL FROM A UNITCELL

# 1d
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
            newBond(BR, 1,1, LBR(1), (+1,)),
            newBond(BR, 1,1, LBR(1), (-1,))
        ]
    )
end

# 2d

# 3d


# export the function
export getReciprocalUnitcell








# TEST IF TWO UNITCELLS ARE RECIPROCAL
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


# TODO LIST
# - getBrillouinZone(rec_unitcell)
# - getBrillouinZoneCorners(rec_unitcell)   (... edges, faces...)
