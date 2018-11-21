################################################################################
#
#	CONCRETE TYPE
#
#   ReciprocalPoint{D} <: AbstractReciprocalPoint{D}
#   --> D is the dimension of embedding space
#
#   FILE CONTAINS
#       - concrete struct definition
#       - interface implementation
#
################################################################################







################################################################################
#
#   CONCRETE STRUCT DEFINITION
#
################################################################################
mutable struct ReciprocalPoint{D} <: AbstractReciprocalPoint{D}

	# point
	point       :: Vector{Float64}

	# label
	label	    :: String
	# LaTeX label
	label_LaTeX :: String

end

# export the concrete type
export ReciprocalPoint





################################################################################
#
#	IMPLEMENTATION OF INTERFACE FOR CONCRETE RECIPROCAL POINT TYPE
#	(functions that had to be overwritten by concrete type)
#
################################################################################

# default constructor interface
# used for creation of new reciprocal points
function newReciprocalPoint(
            :: Type{ReciprocalPoint{D}},
            point   :: Vector{<:Real},
            label   :: AbstractString
        ) :: ReciprocalPoint{D} where {D}

    # check if correct dimension is given
    @assert length(point) == D

    # return a new object
    return ReciprocalPoint{D}(point, label, label)
end
function newReciprocalPoint(
            :: Type{ReciprocalPoint{D}},
            point       :: Vector{<:Real},
            label       :: AbstractString,
            label_LaTeX :: AbstractString
        ) :: ReciprocalPoint{D} where {D}

    # check if correct dimension is given
    @assert length(point) == D

    # return a new object
    return ReciprocalPoint{D}(point, label, label_LaTeX)
end




# get label
function label(
            p :: ReciprocalPoint{D}
        ) :: String where {D}

    # return the label
    return p.label
end
# set label
function label!(
            p :: ReciprocalPoint{D},
            l :: AbstractString
        ) where {D}

    # set the label
    p.label = l
end


# get label (LaTeX)
function labelLaTeX(
            p :: ReciprocalPoint{D}
        ) :: String where {D}

    # return the label
    return p.label_LaTeX
end
# set label (LaTeX)
function labelLaTeX!(
            p :: ReciprocalPoint{D},
            l :: AbstractString
        ) where {D}

    # set the label
    p.label_LaTeX = l
end


# get point
function point(
            p :: ReciprocalPoint{D}
        ) :: Vector{Float64} where {D}

    # return the point
    return p.point
end
# set point
function point!(
            p :: ReciprocalPoint{D},
            r :: Vector{<:Real}
        ) where {D}

    # set the point
    p.point = r
end
