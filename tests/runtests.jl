using Test
using OpenKolmogorovFlow
import Random: seed!
import LinearAlgebra: dot, norm
# using Flows

# include("test_field.jl")      # OK
# include("test_operators.jl")  # OK
# include("test_ftfield.jl")    # OK
# include("test_broadcast.jl")  # OK
# include("test_indexing.jl")   # OK
# include("test_norms.jl")      # OK

# include("test_allocation.jl") # FIXME
# include("test_fft.jl")        # FIXME
# include("test_flow.jl")       # OK
# include("test_spectra.jl")
# include("test_distance.jl")
# include("test_shifts.jl")
include("test_tangent.jl")
# include("test_variational.jl")
# include("test_adjoint_identity.jl")
# include("test_adjoint_allocations.jl")
# # include("test_system.jl")