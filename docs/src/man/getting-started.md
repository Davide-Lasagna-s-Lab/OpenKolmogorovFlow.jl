```@meta
CurrentModule = OpenKolmogorovFlow
```

# Getting Started

This page is the shortest path from a clean Julia environment to a useful
OpenKolmogorovFlow.jl run.

OpenKolmogorovFlow.jl supplies the PDE-specific pieces: field containers,
spectral operators, nonlinear terms, tangent/adjoint equations, forcing
objects, diagnostics, and symmetry maps. Time stepping is delegated to
[Flows.jl](https://github.com/Davide-Lasagna-s-Lab/Flows.jl), which provides
the `flow`, `CB3R2R3e`, `TimeStepConstant`, monitor, storage, and linearised
integration machinery.

## Install

The package is currently used as an unregistered package. Create or activate a
Julia project and develop both OpenKolmogorovFlow.jl and Flows.jl:

```julia
import Pkg

Pkg.activate("kolmogorov-run")
Pkg.develop(url="https://github.com/Davide-Lasagna-s-Lab/Flows.jl.git")
Pkg.develop(url="https://github.com/Davide-Lasagna-s-Lab/OpenKolmogorovFlow.jl.git")
Pkg.instantiate()
```

When working inside a local checkout of this repository, use the package
project directly and develop Flows.jl into it:

```julia
import Pkg

Pkg.activate(".")
Pkg.develop(url="https://github.com/Davide-Lasagna-s-Lab/Flows.jl.git")
Pkg.instantiate()
```

If Julia reports `Package Flows not found`, the Flows.jl dependency has not
been developed or resolved in the active project.

## Build A Field

Physical fields live on a uniform periodic grid with `2m+2` points in each
direction. Fourier fields store a real-to-complex half-plane and are indexed by
[`WaveNumber`](@ref).

```julia
using OpenKolmogorovFlow

n = 8
m = up_dealias_size(n)

u = Field(m, (x, y) -> sin(x) * cos(y))
Ω = FFT(u, n)
u_back = IFFT(Ω)
```

Use `n` for the active spectral cutoff and `m` for storage. For nonlinear
simulations, `m = up_dealias_size(n)` is the usual 3/2-rule choice.

## Run The Forward Equation

The example below starts from a small nonzero vorticity field and integrates at
low Reynolds number until the solution reaches the analytic laminar state. It
is intentionally close to the package test suite, so it exercises the same API
that is kept under regression tests.

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
method = CB3R2R3e(Ω, Flows.NormalMode())
timestep = TimeStepConstant(0.01)
integrator = flow(explicit, implicit, method, timestep)

Ω .= 0.01
Ω[WaveNumber(0, 0)] = 0

integrator(Ω, (0.0, 50.0))

target = laminarflow(n, m, Re, kforcing)
err = normdiff(Ω, target)
```

The important separation of responsibilities is:

- [`ForwardEquation`](@ref) stores the OpenKolmogorovFlow.jl operators and
  work arrays.
- [`splitexim`](@ref) exposes explicit and implicit pieces for an IMEX method.
- `CB3R2R3e` and `TimeStepConstant` come from Flows.jl.
- `flow(...)` comes from Flows.jl and performs the time integration.

## Record Samples

Flows.jl monitors can be used to record states without changing the equation:

```julia
monitor = Monitor(Ω, (t, Ω) -> copy(Ω))
integrator(Ω, (0.0, 10.0), monitor)

states = samples(monitor)
```

Use `copy(Ω)` in the monitor callback when you want snapshots. The state `Ω`
is mutated in place during integration, so storing `Ω` itself would only keep
references to the final object.

## Common Checks

- `m >= n` is required for [`FTField`](@ref) and equation constructors.
- `kforcing` must be representable in the active set, so choose
  `0 <= kforcing <= n`.
- The zero Fourier mode is set to zero by the solver conventions.
- Direct mode indexing should use [`WaveNumber(k, j)`](@ref), not raw storage
  indices, unless you are deliberately working with the packed FFT layout.
