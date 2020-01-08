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
using LatPhysLatticeConstruction
using LinearAlgebra
using HDF5

# imports (to overwrite functions)
import LatPhysBase.label
import LatPhysBase.label!
import LatPhysBase.point
import LatPhysBase.point!
import LatPhysBase.newSite
import LatPhysBase.similar
import LatPhysBase.newUnitcell
import LatPhysBase.sites
import LatPhysBase.sites!
import LatPhysBase.bonds
import LatPhysBase.bonds!
import LatPhysBase.latticeVectors
import LatPhysBase.latticeVectors!
import Base.show




# RECIPROCAL POINTS

# abstract type definition
include("reciprocal_points/abstract_rec_point.jl")
# concrete type definition
include("reciprocal_points/concrete_rec_point.jl")

# Default definitions of various points
include("reciprocal_points/default_points.jl")

# Default definitions of various points
include("reciprocal_points/k_space_points.jl")



# RECIPROCAL PATHS

# abstract type definition
include("reciprocal_paths/abstract_rec_path.jl")
# concrete type definition
include("reciprocal_paths/concrete_rec_path.jl")





# RECIPROCAL LATTICE / UNITCELL / BRILLOUIN ZONE

# abstract type definition (unitcell)
include("reciprocal_lattice/abstract_reciprocal_unitcell.jl")
# concrete type definition (unitcell)
include("reciprocal_lattice/concrete_reciprocal_unitcell.jl")

# abstract type definition (Brillouin zone)
include("reciprocal_lattice/abstract_brillouin_zone.jl")
# concrete type definition (Brillouin zone)
include("reciprocal_lattice/concrete_brillouin_zone.jl")

# construction of reciprocal lattices from real space lattices
include("reciprocal_lattice/functions_reciprocal_lattice.jl")

# construction of Brillouin zones
include("reciprocal_lattice/functions_brillouin_zone.jl")


# module end
end
