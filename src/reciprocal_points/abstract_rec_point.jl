################################################################################
#
#	ABSTRACT TYPE
#
#   ReciprocalPoint{D}
#   --> D is the dimension of embedding space
#
#   FILE CONTAINS
#       - abstract type definition
#       - interface definition
#       - interface testing
#
################################################################################






################################################################################
#
#   ABSTRACT TYPE DEFINITION
#
################################################################################
abstract type AbstractReciprocalPoint{D} <: AbstractSite{String, D} end

# export the type
export AbstractReciprocalPoint



################################################################################
#
#	INTERFACING / ACCESSING RECIPROCAL POINTS
#	(functions have to be overwritten by concrete types)
#
################################################################################

# default constructor interface
# used for creation of new reciprocal points
function newReciprocalPoint(
            :: Type{R},
            point   :: Vector{<:Real},
            label   :: AbstractString
        ) :: R where {R<:AbstractReciprocalPoint{D} where D}

    # print an error because implementation for concrete type is missing
    error(  "not implemented function 'newReciprocalPoint' for concrete reciprocal point type " *
            string(R) * " in reciprocal space dimension " * string(D) )
end

# export the newSite interface
export newReciprocalPoint



# interface of AbstractSite
# so it DOES NOT HAVE TO BE OVERWRITTEN by concrete recipocal point
function newSite(
            :: Type{R},
            point   :: Vector{<:Real},
            label   :: AbstractString
        ) :: R where {D, R <: AbstractReciprocalPoint{D}}

    # return the respective function
    return newReciprocalPoint(R, point, label)
end

# export the newSite interface
export newSite









# get label
function label(
            p :: AbstractReciprocalPoint{D}
        ) :: String where {D}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'label' for reciprocal point type " * string(typeof(p)))
end
# set label
function label!(
            p :: AbstractReciprocalPoint{D},
            l :: AbstractString
        ) where {D}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'label!' for reciprocal point type " * string(typeof(p)))
end

# export the label interface
export label, label!




# get point
function point(
            p :: AbstractReciprocalPoint{D}
        ) :: Vector{Float64} where {D}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'point' for reciprocal point type " * string(typeof(p)))
end
# set point
function point!(
            p :: AbstractReciprocalPoint{D},
            r :: Vector{<:Real}
        ) where {D}

    # print an error because implementation for concrete type is missing
    error("not implemented interface function 'point!' for reciprocal point type " * string(typeof(p)))
end

# export the point interface
export point, point!




# SIMILAR FUNCTION (can be overwritten but does not have to be overwritten)

# without new parameters
function similar(
            p :: R
        ) :: R where {D,R<:AbstractReciprocalPoint{D}}

    # return a new site object
    return newReciprocalPoint(R, deepcopy(point(p)), deepcopy(label(p)))
end
# with new parameters
function similar(
            p :: R,
            r :: Vector{<:Real},
            l :: AbstractString
        ) :: R where {D,R<:AbstractReciprocalPoint{D}}

    # create a new site object
    p_new = similar(p)
    # set parameters
    point!(p_new, r)
    label!(p_new, l)
    # return the new object
    return p_new
end

# export the similar interface
export similar




# LaTeX functionality

# get label
function labelLaTeX(
            p :: AbstractReciprocalPoint{D}
        ) :: String where {D}

    # return the usual label
    return label(p)
end
# set label
function labelLaTeX!(
            p :: AbstractReciprocalPoint{D},
            l :: AbstractString
        ) where {D}

    # set the usual label
    return label!(p,l)
end

# export the label interface
export labelLaTeX, labelLaTeX!



# new point with different LaTeX label
function newReciprocalPoint(
            :: Type{R},
            point      :: Vector{<:Real},
            label      :: AbstractString,
            labelLaTeX :: AbstractString
        ) :: R where {R<:AbstractReciprocalPoint{D} where D}

    # call the usual constructor by neglecting the LaTeX label
    return newReciprocalPoint(R, point, label)
end

# export the newSite interface
export newReciprocalPoint












################################################################################
#
#	OVERWRITING BASE.SHOW FOR RECIPROCAL POINTS
#
################################################################################
function show(io::IO, p::R) where {R<:AbstractReciprocalPoint{D} where D}
    print(io, R, " @", point(p), ": (\"", label(p), "\" / L\"", labelLaTeX(p), "\")")
end
