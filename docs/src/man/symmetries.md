```@meta
CurrentModule = OpenKolmogorovFlow
```

# Symmetries

The package includes in-place symmetry operations for Fourier fields.

## Shifts And Rotations

[`xshift!`](@ref) applies a continuous periodic shift in the `x` direction.
[`yshift!`](@ref) applies the discrete `y`-shift symmetry associated with the
default Kolmogorov forcing wavenumber. [`rotate!`](@ref) applies the
half-turn/conjugation symmetry.

The composed helpers [`shift!`](@ref) and [`shiftrotate!`](@ref) are useful
when comparing states modulo continuous and discrete symmetries.

```julia
Ω_shifted = copy(Ω)
shift!(Ω_shifted, π / 8, 1)

Ω_rotated = copy(Ω)
rotate!(Ω_rotated)
```

The `y`-shift helpers assume the default `kforcing = 4` symmetry class. If you
change the forcing wavenumber, check whether the same discrete symmetry still
represents the problem you are studying.

## Distances And Norms

[`normdiff`](@ref) and [`minnormdiff`](@ref) compare Fourier states. The
symmetry-aware minimum distance searches over streamwise shifts and the
package's discrete shift/rotation operations.

```julia
d = normdiff(Ω1, Ω2)
dmin, (sx, sy_index) = minnormdiff(Ω1, Ω2)
```

`normdiff` returns a squared distance. It is deliberately cheaper than
allocating `Ω1 - Ω2` and then taking a norm.

## Diagnostics

[`laminarflow`](@ref) constructs the analytic laminar vorticity field for a
chosen Reynolds number. [`dissrate`](@ref) and [`powinput`](@ref) evaluate
standard flow diagnostics. [`radial_mean`](@ref) groups Fourier coefficients
by radial wavenumber shells.

```julia
Ω_lam = laminarflow(n, m, Re, kforcing)
ϵ = dissrate(Ω_lam, Re)
P = powinput(Ω_lam, kforcing)
spectrum = radial_mean(abs2, Ω_lam)
```

These diagnostics use the same Fourier weighting as the package inner product,
so they are consistent with the norms used in tangent and adjoint checks.
