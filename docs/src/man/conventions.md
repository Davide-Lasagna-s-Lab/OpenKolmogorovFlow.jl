```@meta
CurrentModule = OpenKolmogorovFlow
```

# Assumptions

OpenKolmogorovFlow.jl works with a two-dimensional vorticity formulation on the
periodic square ``[0,2\pi)^2``. Physical fields are sampled on a uniform
``(2m+2) \times (2m+2)`` grid, and Fourier fields use the packed layout
produced by real-to-complex FFTs.

The package is built for repeated time stepping. Constructors allocate storage
and FFTW plans; repeated calls to the operators mutate existing arrays.

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

For nonlinear simulations, the standard choice is:

```julia
n = 32
m = up_dealias_size(n)
```

For tests, linear operations, or deliberately undealiased experiments, `m = n`
is valid.

## Coordinates And Wavenumbers

Fourier modes are addressed with [`WaveNumber`](@ref). The package convention
is:

- `WaveNumber(k, j)` stores the signed full-FFT direction in `k`.
- `j` is the non-negative real-to-complex direction.
- Negative `j` modes are represented implicitly by Hermitian symmetry.

The differential operators follow this convention. [`ddx!`](@ref) multiplies
by ``i j`` and [`ddy!`](@ref) multiplies by ``i k``.

The concrete storage is an ordinary Julia matrix. The first storage dimension
contains signed `k` modes in FFT order and the second storage dimension
contains non-negative `j` modes:

```text
storage rows:     k = 0, 1, 2, ..., m, -(m+1), ..., -2, -1
storage columns:  j = 0, 1, 2, ..., m, m+1
```

Use [`WaveNumber`](@ref) for user-facing mode access:

```julia
Ω[WaveNumber(-2, 1)]
```

Raw integer indexing is available because [`FTField`](@ref) is an array, but
it exposes storage coordinates rather than physical wavenumber intent.

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

This convention is important for the streamfunction solve: the inverse
Laplacian is only defined after fixing the mean. If a calculation requires a
nonzero mean, keep that mode separate from the vorticity dynamics.

## Relationship To Flows.jl

OpenKolmogorovFlow.jl defines equation objects. It does not own the outer time
integration loop. Use [Flows.jl](https://github.com/Davide-Lasagna-s-Lab/Flows.jl)
for:

- `flow(...)`, the propagation object;
- IMEX methods such as `CB3R2R3e`;
- time-step policies such as `TimeStepConstant` and storage-driven policies;
- monitors and trajectory storage;
- tangent/adjoint propagation modes.

The Flows.jl manual is available at
[davide-lasagna-s-lab.github.io/Flows.jl](https://davide-lasagna-s-lab.github.io/Flows.jl/stable/).
