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
    latvec_b1 = a1(unitcell) .* 1/sqrt(dot(a1(unitcell), a1(unitcell)))

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







# export the function
export getReciprocalUnitcell




# TODO LIST
# - getReciprocalUnitcell(unitcell)
# - getBrillouinZone(rec_unitcell)
# - getBrillouinZoneCorners(rec_unitcell)   (... edges, faces...)
