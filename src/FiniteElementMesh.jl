module FiniteElementMesh

using LinearAlgebra, SpecialPolynomials

import Base: show, read, write

include("types.jl")
include("relation.jl")
include("geometry.jl")
include("partition.jl")
include("interpolation.jl")

end # module FiniteElementMesh
