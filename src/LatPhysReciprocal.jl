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



# RECIPROCAL POINTS

# abstract type definition
include("reciprocal_points/abstract_rec_point.jl")



# module end
end
