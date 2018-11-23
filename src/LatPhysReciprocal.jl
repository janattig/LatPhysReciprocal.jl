################################################################################
#
#   module LatPhysReciprocal
#   -> LatticePhysics_Base
#   -> LinearAlgebra
#
#   --> TYPE DEFINITIONS (Abstract / Concrete)
#           - Point
#           - Path
#           - BrillouinZone
#           - ReciprocalLattice ??
#
#   --> BUILDING FUNCTIONS FOR PATHS IN MOMENTUM SPACE
#           - modification of paths
#           - pre-implemented paths
#           - paths based on Unitcells
#
#   --> BRILLOUIN ZONES
#           - type definition
#           - pre-implemented BZs
#           - construction of default BZs (2D)
#           - construction of default BZs (3D)
#
################################################################################




# module start
module LatPhysReciprocal




# used libraries
using LatPhysBase

# imports (to overwrite functions)
import LatPhysBase.label
import LatPhysBase.label!
import LatPhysBase.point
import LatPhysBase.point!
import LatPhysBase.newSite
import LatPhysBase.similar
import Base.show




# RECIPROCAL POINTS

# abstract type definition
include("reciprocal_points/abstract_rec_point.jl")
# concrete type definition
include("reciprocal_points/concrete_rec_point.jl")

# Default definitions of various points
include("reciprocal_points/default_points.jl")




# RECIPROCAL PATHS

# abstract type definition
include("reciprocal_paths/abstract_rec_path.jl")
# concrete type definition
include("reciprocal_paths/concrete_rec_path.jl")




# module end
end
