################################################################################
#
#	CONCRETE TYPE
#
#   ReciprocalPath{P} <: AbstractReciprocalPath{P}
#   --> P is the reciprocal point type (<: AbstractReciprocalPoint{D})
#       --> D is the dimension of embedding space
#
#   FILE CONTAINS
#       - concrete struct definition
#       - interface implementation
#
################################################################################



################################################################################
#
#   STRUCT DEFINITION
#
################################################################################
mutable struct ReciprocalPath{P} <: AbstractReciprocalPath{P}

    # the only field is the array of points
    points :: Vector{P}

end

# export the type
export ReciprocalPath






################################################################################
#
#	INTERFACING / ACCESSING RECIPROCAL PATHS
#
################################################################################


# default constructor interface
# used for creation of new reciprocal paths
function newReciprocalPath(
            ::Type{ReciprocalPath{RPO}},
            points :: Vector{<:RPO}
        ) :: ReciprocalPath{RPO} where {D,RPO<:AbstractReciprocalPoint{D}}

    # return the newly created object
    return ReciprocalPath{RPO}(points)
end

# export the newUnitcell interface
export newReciprocalPath






# access the list of points
function points(
            p :: ReciprocalPath{RPO}
        ) :: Vector{RPO} where {D,RPO<:AbstractReciprocalPoint{D}}

    # return the respective field
    return p.points
end
# set the list of points
function points!(
            p :: ReciprocalPath{RPO},
            points :: Vector{<:RPO}
        ) where {D,RPO<:AbstractReciprocalPoint{D}}

    # set the respective field
    p.points = points
end


# export access to points
export points, points!









function saveReciprocalPath(
        rp :: RP,
        fn :: AbstractString,
        group :: AbstractString = "reciprocal_path"
        ;
        append :: Bool = false
    ) where {D,P<:ReciprocalPoint{D},RP<:ReciprocalPath{P}}

    # determine the mode based on if one wants to append stuff
    if append
        mode = "r+"
    else
        mode = "w"
    end

    # open the file in mode
    h5open(fn, mode) do file
        # create the group in which the bonds are saved
        group_path = g_create(file, group)
        # save the parameters
        attrs(group_path)["D"] = Int64(D)
        # save all Positions (D dimensions)
        if 0 < Int64(D)
            for n in 1:Int64(D)
                group_path["point_$(n)"] = Float64[point(p)[n] for p in rp]
            end
        end
        # save all labels
        group_path["label"]       = String[string(label(p)) for p in rp]
        group_path["label_LaTeX"] = String[string(labelLaTeX(p)) for p in rp]
    end

    # return nothing
    return nothing
end

function loadReciprocalPath(
        ::Type{RPI},
        fn :: AbstractString,
        group :: AbstractString = "reciprocal_path"
    ) where {DI,PI<:ReciprocalPoint{DI},RPI<:Union{ReciprocalPath{PI},ReciprocalPath}}

    # read attribute data
    attr_data = h5readattr(fn, group)
    # determine D based on this
    D  = attr_data["D"]

    # load all remaining data
    rp_label       = h5read(fn, group*"/label")
    rp_label_LaTeX = h5read(fn, group*"/label_LaTeX")

    if D == 0
        rp_point = [Float64[] for i in 1:length(rp_label)]
    else
        rp_point_parts = [
            h5read(fn, group*"/point_"*string(j)) for j in 1:D
        ]
        rp_point = Vector{Float64}[Float64[rp_point_parts[j][i] for j in 1:D] for i in 1:length(rp_label)]
    end

    # create list of points in the path
    rpp = ReciprocalPoint{D}[
        newReciprocalPoint(ReciprocalPoint{D}, rp_point[i], rp_label[i], rp_label_LaTeX[i])
        for i in 1:length(rp_label)
    ]

    # return the path
    return newReciprocalPath(ReciprocalPath{ReciprocalPoint{D}}, rpp)
end


# convinience function for standard type
function loadReciprocalPath(
        fn :: AbstractString,
        group :: AbstractString = "reciprocal_path"
    )

    return loadReciprocalPath(ReciprocalPath, fn, group)
end
