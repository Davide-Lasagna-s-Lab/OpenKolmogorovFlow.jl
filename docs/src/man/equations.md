```@meta
CurrentModule = OpenKolmogorovFlow
```

# Equations

OpenKolmogorovFlow.jl represents equations as callable objects compatible with
the [Flows.jl](https://github.com/Davide-Lasagna-s-Lab/Flows.jl)
time-stepping interface. Use this package for the PDE operators, and use
Flows.jl for the propagation loop.

## Forward Equation

[`ForwardEquation`](@ref) combines three pieces:

- an [`ImplicitTerm`](@ref) for viscous diffusion;
- a [`ForwardExplicitTerm`](@ref) for nonlinear advection and built-in
  Kolmogorov forcing;
- an optional [`AbstractForcing`](@ref) object for extra forcing.

The nonlinear term computes derivatives and velocities spectrally, transforms
to physical space for products, then transforms the result back to Fourier
space. The explicit term also updates the CFL estimate used by [`CFLHook`](@ref).

The smallest useful pattern is:

```julia
using Flows
using OpenKolmogorovFlow

n = 10
m = up_dealias_size(n)
Re = 1.25
kforcing = 4

Ω = FTField(n, m)
equation = ForwardEquation(n, m, Re, kforcing)
explicit, implicit = splitexim(equation)

integrator = flow(explicit,
                  implicit,
                  CB3R2R3e(Ω, Flows.NormalMode()),
                  TimeStepConstant(0.01))

integrator(Ω, (0.0, 1.0))
```

The state `Ω` is mutated in place. Allocate a copy in a monitor if you need to
store intermediate states.

## Split Integration

[`splitexim`](@ref) returns explicit and implicit pieces for IMEX methods. The
explicit closure mutates the supplied right-hand side, while the implicit term
implements multiplication and shifted inversion for the linear viscous part.

This is the bridge to Flows.jl. `splitexim(equation)` returns exactly the two
objects needed by IMEX methods:

```julia
explicit, implicit = splitexim(ForwardEquation(n, m, Re))
```

The explicit part computes advection and forcing. The implicit part represents
viscous diffusion and supports the shifted inverse operation required by the
IMEX schemes.

## Tangent And Adjoint Equations

[`LinearisedEquation`](@ref) uses [`TangentMode`](@ref) or [`AdjointMode`](@ref)
to choose the linearised operator. The tangent equation evolves perturbations
about a forward state. The adjoint equation uses the package inner product and
supports adjoint forcing objects such as [`DissRateGradientForcing`](@ref).

```julia
using FFTW

tangent = LinearisedEquation(n, m, Re, TangentMode())
adjoint = LinearisedEquation(n, m, Re, AdjointMode(),
                             FFTW.EXHAUSTIVE,
                             DissRateGradientForcing(n, Re))
```

Flows.jl supplies the linearised propagation modes and the storage/cache
machinery needed to replay a primal trajectory. See the
[Flows.jl linearised dynamics docs](https://davide-lasagna-s-lab.github.io/Flows.jl/stable/linearised/)
for the time-stepping side of that workflow.

## Forcing Objects

Forcings are small callable objects that mutate the destination right-hand side:

- [`DummyForcing`](@ref) leaves the right-hand side unchanged.
- [`SteadyForcing`](@ref) adds a stored Fourier field.
- [`WaveNumberForcing`](@ref) adds a mode-local forcing.
- [`ReForcing`](@ref) provides the tangent forcing associated with a Reynolds
  number perturbation.
- [`DissRateGradientForcing`](@ref) provides the adjoint forcing for
  dissipation-rate gradients.

For example, to add a single physical cosine forcing through its Fourier mode:

```julia
forcing = WaveNumberForcing(n, 1, 0, 1)
dΩdt = FTField(n, m)
forcing(0.0, Ω, Ω, dΩdt)
```

Forcings are additive by convention. They receive a destination right-hand-side
field and add their contribution to whatever is already stored there.
