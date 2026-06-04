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

The callable constructor is the most convenient way to set up exact test data:

```julia
m = 12
u = Field(m, (x, y) -> cos(2x + y))
```

The arrays returned by [`make_grid`](@ref) are shaped for broadcasting, so the
same expression can also be written explicitly:

```julia
x, y = make_grid(m)
u = Field(cos.(2x .+ y))
```

## Fourier Fields

[`FTField(n, m)`](@ref) allocates Fourier storage with active cutoff `n` and
storage cutoff `m`. Individual modes are indexed through [`WaveNumber`](@ref):

```julia
Ω = FTField(8, up_dealias_size(8))
Ω[WaveNumber(2, 1)] = 1im
```

For a real physical field, coefficients with negative `j` are not stored
explicitly. The missing half-plane is implied by Hermitian symmetry. The
following coefficients represent `sin(y)`:

```julia
Ω = FTField(4, 4)
Ω[WaveNumber( 1, 0)] = -0.5im
Ω[WaveNumber(-1, 0)] =  0.5im
u = IFFT(Ω)
```

The helper [`growto!`](@ref) copies active coefficients into larger Fourier
storage and clears inactive modes. This is useful when moving between active
and dealiased representations.

```julia
n = 16
small = FTField(n, n)
large = FTField(n, up_dealias_size(n))
growto!(large, small)
```

## Transforms

[`ForwardFFT!`](@ref) and [`InverseFFT!`](@ref) build reusable FFTW plans. The
allocating helpers [`FFT`](@ref) and [`IFFT`](@ref) are convenient for one-off
conversions or tests.

Forward transforms normalise by the number of physical grid points and apply
the active mask. Inverse transforms mask the input before calling FFTW, which
matches the solver convention that inactive storage modes are scratch space.

For repeated transforms, build plans once:

```julia
u = Field(m)
Ω = FTField(n, m)

fft! = ForwardFFT!(u)
ifft! = InverseFFT!(Ω)

fft!(Ω, u)
ifft!(u, Ω)
```

The allocating helpers copy their inputs before planning because FFTW planning
may destroy the arrays used to create plans.

## Spectral Operators

The basic spectral operators are:

- [`ddx!`](@ref) and [`ddy!`](@ref) for first derivatives.
- [`laplacian!`](@ref) for ``\Delta``.
- [`invlaplacian!`](@ref) for the zero-mean inverse elliptic operators.

All of these are diagonal in Fourier space and require matching input/output
resolutions.

```julia
Ω = FFT(Field(m, (x, y) -> sin(2x + y)), n)
dΩdx = similar(Ω)
dΩdy = similar(Ω)

ddx!(dΩdx, Ω)
ddy!(dΩdy, Ω)
```

Remember that the package convention is `ddx! -> im*j` and `ddy! -> im*k`.
This follows the naming used by the codebase even though the packed storage
uses `k` along the first matrix dimension.
