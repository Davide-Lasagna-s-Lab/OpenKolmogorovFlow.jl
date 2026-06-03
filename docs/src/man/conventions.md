```@meta
CurrentModule = OpenKolmogorovFlow
```

# Assumptions

OpenKolmogorovFlow.jl works with a two-dimensional vorticity formulation on the
periodic square ``[0,2\pi)^2``. Physical fields are sampled on a uniform
``(2m+2) \times (2m+2)`` grid, and Fourier fields use the packed layout
produced by real-to-complex FFTs.

## Resolution Parameters

Two integer cutoffs appear throughout the API:

- `n` is the active spectral cutoff. Solver loops evolve modes with
  ``j = 0,\ldots,n`` and ``k = -n,\ldots,n``.
- `m` is the storage cutoff. The physical grid has `2m+2` points per
  direction, and the stored Fourier array has size `(2m+2, m+2)`.

The storage cutoff may be larger than the active cutoff. This is how the
package represents dealiased nonlinear products: evolve the active modes, use a
larger physical grid for products, and mask inactive modes after transforming
back to Fourier space.

## Coordinates And Wavenumbers

Fourier modes are addressed with [`WaveNumber`](@ref). The package convention
is:

- `WaveNumber(k, j)` stores the signed full-FFT direction in `k`.
- `j` is the non-negative real-to-complex direction.
- Negative `j` modes are represented implicitly by Hermitian symmetry.

The differential operators follow this convention. [`ddx!`](@ref) multiplies
by ``i j`` and [`ddy!`](@ref) multiplies by ``i k``.

## Mutating Operators

Most expensive operations are written as `f!(out, in)` and return `out`. This
style is used for transforms, derivatives, elliptic solves, shifts, and
time-stepping terms. Solver objects such as [`ForwardEquation`](@ref) and
[`LinearisedEquation`](@ref) own their work arrays and FFT plans so repeated
calls can be allocation-conscious.

## Mean And Symmetry

The zero mode is treated specially in several operators. Inverse Laplacians set
the mean to zero, and FFT plans apply masking and Hermitian cleanup so inactive
or redundant modes do not leak into the active dynamics.
