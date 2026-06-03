"""
    OpenKolmogorovFlow

Spectral solver utilities for the two-dimensional open Kolmogorov-flow
vorticity equation on a periodic square domain.

The package stores real-space fields on a uniform ``(2m+2) × (2m+2)`` grid and
stores vorticity in real-to-complex Fourier layout through [`FTField`](@ref).
The logical active cutoff `n` controls which Fourier modes are evolved, while
the storage cutoff `m` may be larger to provide dealiasing. Operators,
nonlinear terms, tangent equations, adjoint equations, forcings, symmetries,
and CFL-based timestep hooks are all built around this pair of resolutions.

Most mutating functions follow Julia's `f!(out, in)` convention and return the
modified output. Many solver objects carry preallocated work arrays so repeated
time stepping can avoid transient allocations.
"""
module OpenKolmogorovFlow

include("indexing.jl")
include("ftfield.jl")
include("field.jl")
include("operators.jl")
include("fft.jl")
include("norms.jl")
include("implicit.jl")
include("forcings.jl")
include("nonlinear.jl")
include("linearised.jl")
include("hooks.jl")
include("flow.jl")
include("spectra.jl")
include("shifts.jl")
# include("distance.jl")

end
