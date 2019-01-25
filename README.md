# LatPhysReciprocal.jl

Reciprocal space tools and types for [`LatticePhysics.jl`](https://github.com/janattig/LatticePhysics.jl).



## Contents

Definition of concrete and abstract types for
1.  Reciprocal Point
2.  Reciprocal Path (composed of multiple reciprocal points)
3.  Reciprocal Unitcell
4.  (First) Brillouin Zone

Definition of functions to construct all objects.




## Installation

You can install the package via the package mode in Julia (Pkg). However, since the package
is not listed in the Julia package repositories, you have to first install the unregistered
dependencies manually with
```julia
(v1.0) pkg> add "https://github.com/janattig/LatPhysBase.jl"
(v1.0) pkg> add "https://github.com/janattig/LatPhysLatticeConstruction.jl"
```
to finally install the main package with
```julia
(v1.0) pkg> add "https://github.com/janattig/LatPhysReciprocal.jl"
```
