```@meta
CurrentModule = OpenKolmogorovFlow
```

# Equations

OpenKolmogorovFlow.jl represents equations as callable objects compatible with
the `Flows.jl` time-stepping interface.

## Forward Equation

[`ForwardEquation`](@ref) combines three pieces:

- an [`ImplicitTerm`](@ref) for viscous diffusion;
- a [`ForwardExplicitTerm`](@ref) for nonlinear advection and built-in
  Kolmogorov forcing;
- an optional [`AbstractForcing`](@ref) object for extra forcing.

The nonlinear term computes derivatives and velocities spectrally, transforms
to physical space for products, then transforms the result back to Fourier
space. The explicit term also updates the CFL estimate used by [`CFLHook`](@ref).

## Split Integration

[`splitexim`](@ref) returns explicit and implicit pieces for IMEX methods. The
explicit closure mutates the supplied right-hand side, while the implicit term
implements multiplication and shifted inversion for the linear viscous part.

## Tangent And Adjoint Equations

[`LinearisedEquation`](@ref) uses [`TangentMode`](@ref) or [`AdjointMode`](@ref)
to choose the linearised operator. The tangent equation evolves perturbations
about a forward state. The adjoint equation uses the package inner product and
supports adjoint forcing objects such as [`DissRateGradientForcing`](@ref).

## Forcing Objects

Forcings are small callable objects that mutate the destination right-hand side:

- [`DummyForcing`](@ref) leaves the right-hand side unchanged.
- [`SteadyForcing`](@ref) adds a stored Fourier field.
- [`WaveNumberForcing`](@ref) adds a mode-local forcing.
- [`ReForcing`](@ref) provides the tangent forcing associated with a Reynolds
  number perturbation.
- [`DissRateGradientForcing`](@ref) provides the adjoint forcing for
  dissipation-rate gradients.
