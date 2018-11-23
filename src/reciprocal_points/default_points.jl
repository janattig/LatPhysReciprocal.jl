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
        ) :: R where {I1,I2,D,N,LS,LB, R<:AbstractReciprocalPoint{D}, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # fallback / fail due to missing implementation
    error("Implementation of reciprocal point " * string(I1) * " (instance " * string(I2) * ") not implemented yet " *
    "for concrete reciprocal point type " * string(R) * " and unitcells of dimension D=" * string(D) * " and N="*string(N))
end


# export the general interface function
export getReciprocalPoint





################################################################################
#
#	IMPLEMENTATIONS OF DEFAULT POINTS IN RECIPROCAL SPACE
#
################################################################################




#########################
#
#	GAMMA POINT
#
#########################

# Gamma point in arbitrary dimension (Fallback)
function getReciprocalPoint(
            :: Type{R},
            unitcell   :: U,
            identifier :: Val{:Gamma},
            instance   :: Val{I}
        ) :: R where {I,D,N,LS,LB, R<:AbstractReciprocalPoint{D}, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # fallback / fail due to missing implementation
    error("Implementation of Gamma point in " * string(D) * " spatial dimensions missing " *
    "for concrete reciprocal point type " * string(R))
end

# Gamma point in 2D
function getReciprocalPoint(
            :: Type{R},
            unitcell   :: U,
            identifier :: Val{:Gamma},
            instance   :: Val{I}
        ) :: R where {I,N,LS,LB, R<:AbstractReciprocalPoint{2}, S<:AbstractSite{LS,2}, B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # return the specific type
    return newReciprocalPoint(
        R,
        [0.0, 0.0],
        "Gamma",
        "\\Gamma"
    )
end

# Gamma point in 3D
function getReciprocalPoint(
            :: Type{R},
            unitcell   :: U,
            identifier :: Val{:Gamma},
            instance   :: Val{I}
        ) :: R where {I,N,LS,LB, R<:AbstractReciprocalPoint{3}, S<:AbstractSite{LS,3}, B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # return the specific type
    return newReciprocalPoint(
        R,
        [0.0, 0.0, 0.0],
        "Gamma",
        "\\Gamma"
    )
end
