<p align="center">
  <img src="docs/src/assets/logo.png" alt="OpenKolmogorovFlow.jl logo" width="360">
</p>

# OpenKolmogorovFlow.jl

OpenKolmogorovFlow.jl provides spectral building blocks for the two-dimensional
open Kolmogorov-flow vorticity equation on a periodic square. It includes
Fourier and physical field containers, dealiased FFT transforms, spectral
derivatives, nonlinear terms, tangent and adjoint equations, forcings,
diagnostics, and symmetry operations.

The package is designed around mutating, preallocated operators so that
time-stepping workflows can keep allocations under control.

Documentation is available at:

https://Davide-Lasagna-s-Lab.github.io/OpenKolmogorovFlow.jl/dev/
