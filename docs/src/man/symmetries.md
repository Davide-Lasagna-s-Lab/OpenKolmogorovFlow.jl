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

## Distances And Norms

[`normdiff`](@ref) and [`minnormdiff`](@ref) compare Fourier states. The
symmetry-aware minimum distance searches over streamwise shifts and the
package's discrete shift/rotation operations.

## Diagnostics

[`laminarflow`](@ref) constructs the analytic laminar vorticity field for a
chosen Reynolds number. [`dissrate`](@ref) and [`powinput`](@ref) evaluate
standard flow diagnostics. [`radial_mean`](@ref) groups Fourier coefficients
by radial wavenumber shells.
