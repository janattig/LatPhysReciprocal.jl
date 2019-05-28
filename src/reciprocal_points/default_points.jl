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

    # maybe instance is passed from another source in some notation as e.g. X_2
    if instance == 1 && occursin('_', string(identifier))
        return getReciprocalPoint(R, unitcell, Val(Symbol(split(string(identifier),'_')[1])), Val(Meta.eval(Meta.parse(split(string(identifier),'_')[2]))))
    end
    # return the specific val typed function
    return getReciprocalPoint(R, unitcell, Val(identifier), Val(instance))
end
function getReciprocalPoint(
            :: Type{R},
            unitcell   :: U,
            identifier :: Tuple{Symbol, Int64}
        ) :: R where {D,N,LS,LB, R<:AbstractReciprocalPoint{D}, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B}}

    # return the specific val typed function
    return getReciprocalPoint(R, unitcell, Val(identifier[1]), Val(identifier[2]))
end

# interface for own concrete reciprocal point type
function getReciprocalPoint(
            unitcell   :: U,
            identifier :: Symbol,
            instance   :: Int64 = 1
        ) :: ReciprocalPoint{D} where {D,N,LS,LB, S<:AbstractSite{LS,D}, B<:AbstractBond{LB,N}, U<:AbstractUnitcell{S,B}}

    # maybe instance is passed from another source in some notation as e.g. X_2
    if instance == 1 && occursin('_', string(identifier))
        return getReciprocalPoint(ReciprocalPoint{D}, unitcell, Val(Symbol(split(string(identifier),'_')[1])), Val(Meta.eval(Meta.parse(split(string(identifier),'_')[2]))))
    end
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







#########################
#
#	K POINT (and K')
#
#########################

# K point (Based on BZ object) in 2D
function getReciprocalPoint(
            :: Type{R},
            unitcell   :: U,
            identifier :: Val{:K},
            instance   :: Val{I}
        ) :: R where {D,I,N,LS,LB, R<:AbstractReciprocalPoint{2}, S<:AbstractSite{LS,2}, B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # get the bz
    bz = getBrillouinZone(getReciprocalUnitcell(unitcell))

    # return the specific type
    return newReciprocalPoint(
        R,
        corners(bz)[faces(bz)[1][mod(I,length(corners(bz))) + 1]],
        "K",
        "K"
    )
end

# K' point (Based on BZ object) in 2D
function getReciprocalPoint(
            :: Type{R},
            unitcell   :: U,
            identifier :: Val{:Kp},
            instance   :: Val{I}
        ) :: R where {D,I,N,LS,LB, R<:AbstractReciprocalPoint{2}, S<:AbstractSite{LS,2}, B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # get the bz
    bz = getBrillouinZone(getReciprocalUnitcell(unitcell))

    # return the specific type
    return newReciprocalPoint(
        R,
        corners(bz)[faces(bz)[1][mod(I+1,length(corners(bz))) + 1]],
        "K\'",
        "K\'"
    )
end



# K point (Based on BZ object) in 3D
function getReciprocalPoint(
            :: Type{R},
            unitcell   :: U,
            identifier :: Val{:K},
            instance   :: Val{I}
        ) :: R where {D,I,N,LS,LB, R<:AbstractReciprocalPoint{3}, S<:AbstractSite{LS,3}, B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # get the bz
    bz = getBrillouinZone(getReciprocalUnitcell(unitcell))

    # get the faceindex
    faceindex = mod(div(I, 10), length(faces(bz))) + 1
    face = faces(bz)[faceindex]
    # first and second index
    cornerindex = mod(I, 10)
    i = face[mod(cornerindex, length(face))+1]

    # return the specific type
    return newReciprocalPoint(
        R,
        corners(bz)[i],
        "K",
        "K"
    )
end

# K' point (Based on BZ object) in 3D
function getReciprocalPoint(
            :: Type{R},
            unitcell   :: U,
            identifier :: Val{:Kp},
            instance   :: Val{I}
        ) :: R where {D,I,N,LS,LB, R<:AbstractReciprocalPoint{3}, S<:AbstractSite{LS,3}, B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # get the bz
    bz = getBrillouinZone(getReciprocalUnitcell(unitcell))

    # get the faceindex
    faceindex = mod(div(I, 10), length(faces(bz))) + 1
    face = faces(bz)[faceindex]
    # first and second index
    cornerindex = mod(I, 10)
    i = face[mod(cornerindex+1, length(face))+1]

    # return the specific type
    return newReciprocalPoint(
        R,
        corners(bz)[i],
        "K\'",
        "K\'"
    )
end






#########################
#
#	M POINT (and M')
#
#########################

# M point (Based on BZ object) in 2D
function getReciprocalPoint(
            :: Type{R},
            unitcell   :: U,
            identifier :: Val{:M},
            instance   :: Val{I}
        ) :: R where {D,I,N,LS,LB, R<:AbstractReciprocalPoint{2}, S<:AbstractSite{LS,2}, B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # get the bz
    bz = getBrillouinZone(getReciprocalUnitcell(unitcell))

    # first and second index
    i1 = faces(bz)[1][mod(I-1,length(faces(bz)[1])) + 1]
    i2 = faces(bz)[1][mod(I,length(faces(bz)[1])) + 1]
    # find the first and next corner
    c1 = corners(bz)[i1]
    c2 = corners(bz)[i2]

    # return the specific type
    return newReciprocalPoint(
        R,
        (c1.+c2)./2,
        "M",
        "M"
    )
end

# M' point (Based on BZ object) in 2D
function getReciprocalPoint(
            :: Type{R},
            unitcell   :: U,
            identifier :: Val{:Mp},
            instance   :: Val{I}
        ) :: R where {D,I,N,LS,LB, R<:AbstractReciprocalPoint{2}, S<:AbstractSite{LS,2}, B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # get the bz
    bz = getBrillouinZone(getReciprocalUnitcell(unitcell))

    # first and second index
    i1 = faces(bz)[1][mod(I,length(faces(bz)[1])) + 1]
    i2 = faces(bz)[1][mod(I+1,length(faces(bz)[1])) + 1]
    # find the first and next corner
    c1 = corners(bz)[i1]
    c2 = corners(bz)[i2]

    # return the specific type
    return newReciprocalPoint(
        R,
        (c1.+c2)./2,
        "M\'",
        "M\'"
    )
end



# M point (Based on BZ object) in 3D
function getReciprocalPoint(
            :: Type{R},
            unitcell   :: U,
            identifier :: Val{:M},
            instance   :: Val{I}
        ) :: R where {D,I,N,LS,LB, R<:AbstractReciprocalPoint{3}, S<:AbstractSite{LS,3}, B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # get the bz
    bz = getBrillouinZone(getReciprocalUnitcell(unitcell))

    # get the faceindex
    faceindex = mod(div(I, 10), length(faces(bz))) + 1
    face = faces(bz)[faceindex]
    # first and second index
    cornerindex = mod(I, 10)
    i1 = face[mod(cornerindex-1, length(face)) + 1]
    i2 = face[mod(cornerindex,   length(face)) + 1]
    # find the first and next corner
    c1 = corners(bz)[i1]
    c2 = corners(bz)[i2]

    # return the specific type
    return newReciprocalPoint(
        R,
        (c1.+c2)./2,
        "M",
        "M"
    )
end

# M' point (Based on BZ object) in 3D
function getReciprocalPoint(
            :: Type{R},
            unitcell   :: U,
            identifier :: Val{:Mp},
            instance   :: Val{I}
        ) :: R where {D,I,N,LS,LB, R<:AbstractReciprocalPoint{3}, S<:AbstractSite{LS,3}, B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # get the bz
    bz = getBrillouinZone(getReciprocalUnitcell(unitcell))

    # get the faceindex
    faceindex = mod(div(I, 10), length(faces(bz))) + 1
    face = faces(bz)[faceindex]
    # first and second index
    cornerindex = mod(I, 10)
    i1 = face[mod(cornerindex, length(face)) + 1]
    i2 = face[mod(cornerindex+1,   length(face)) + 1]
    # find the first and next corner
    c1 = corners(bz)[i1]
    c2 = corners(bz)[i2]

    # return the specific type
    return newReciprocalPoint(
        R,
        (c1.+c2)./2,
        "M\'",
        "M\'"
    )
end







#########################
#
#	X POINT (center of face)
#
#########################

# X point (Based on BZ object) in 3D
function getReciprocalPoint(
            :: Type{R},
            unitcell   :: U,
            identifier :: Val{:X},
            instance   :: Val{I}
        ) :: R where {D,I,N,LS,LB, R<:AbstractReciprocalPoint{3}, S<:AbstractSite{LS,3}, B<:AbstractBond{LB,N},U<:AbstractUnitcell{S,B}}

    # get the bz
    bz = getBrillouinZone(getReciprocalUnitcell(unitcell))

    # get the faceindex
    faceindex = mod(I-1, length(faces(bz))) + 1
    face = faces(bz)[faceindex]
    # find the center
    center = [0.0, 0.0, 0.0]
    for c in face
        center .+= corners(bz)[c]
    end
    center ./= length(face)

    # return the specific type
    return newReciprocalPoint(
        R,
        center,
        "X",
        "X"
    )
end
