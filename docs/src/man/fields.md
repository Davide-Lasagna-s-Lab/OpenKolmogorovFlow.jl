```@meta
CurrentModule = OpenKolmogorovFlow
```

# Fields

The package has two primary field containers:

- [`Field`](@ref) stores a real physical-space scalar field.
- [`FTField`](@ref) stores a complex Fourier-space field in packed `rfft`
  layout.

Both behave like Julia arrays. Use `parent(field)` when direct access to the
underlying storage is needed.

## Physical Fields

[`Field(m)`](@ref) allocates a zero field on a uniform grid with `2m+2` points
in each direction. [`make_grid`](@ref) returns broadcastable coordinate arrays
for evaluating analytic functions on this grid.

```julia
u = Field(8) do x, y
    sin(x) * cos(y)
end
```

## Fourier Fields

[`FTField(n, m)`](@ref) allocates Fourier storage with active cutoff `n` and
storage cutoff `m`. Individual modes are indexed through [`WaveNumber`](@ref):

```julia
Ω = FTField(8, up_dealias_size(8))
Ω[WaveNumber(2, 1)] = 1im
```

The helper [`growto!`](@ref) copies active coefficients into larger Fourier
storage and clears inactive modes. This is useful when moving between active
and dealiased representations.

## Transforms

[`ForwardFFT!`](@ref) and [`InverseFFT!`](@ref) build reusable FFTW plans. The
allocating helpers [`FFT`](@ref) and [`IFFT`](@ref) are convenient for one-off
conversions or tests.

Forward transforms normalise by the number of physical grid points and apply
the active mask. Inverse transforms mask the input before calling FFTW, which
matches the solver convention that inactive storage modes are scratch space.

## Spectral Operators

The basic spectral operators are:

- [`ddx!`](@ref) and [`ddy!`](@ref) for first derivatives.
- [`laplacian!`](@ref) for ``\Delta``.
- [`invlaplacian!`](@ref) for the zero-mean inverse elliptic operators.

All of these are diagonal in Fourier space and require matching input/output
resolutions.
