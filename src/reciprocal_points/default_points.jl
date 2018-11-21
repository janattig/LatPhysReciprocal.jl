################################################################################
#
#	INTERFACE FOR DEFAULT POINTS IN RECIPROCAL SPACE
#	(implementations see below)
#
################################################################################

# general interface
function getReciprocalPoint(
            :: Type{R},
            unitcell   :: U,
            identifier :: Symbol,
            instance   :: Int64 = 1
        ) :: R where {D,N,LS,LB, R<:AbstractReciprocalPoint{D}, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B}}

    # return the specific val typed function
    return getReciprocalPoint(R, unitcell, Val(identifier), Val(instance))
end

# interface for own concrete reciprocal point type
function getReciprocalPoint(
            unitcell   :: U,
            identifier :: Symbol,
            instance   :: Int64 = 1
        ) :: ReciprocalPoint{D} where {D,N,LS,LB, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B}}

    # return the specific val typed function
    return getReciprocalPoint(ReciprocalPoint{D}, unitcell, Val(identifier), Val(instance))
end



# Fallback for all functions
function getReciprocalPoint(
            :: Type{R},
            unitcell   :: U,
            identifier :: Val{I1},
            instance   :: Val{I2}
        ) :: R where {I1,I2,D,N,LS,LB, R<:AbstractReciprocalPoint{D}, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B},V}

    # fallback / fail due to missing implementation
    error("Implementation of reciprocal point " * string(I1) * " (instance " * string(I2) * ") not implemented yet " *
    "for concrete reciprocal point type " * string(R))
end


# export the general interface function
export getReciprocalPoint
