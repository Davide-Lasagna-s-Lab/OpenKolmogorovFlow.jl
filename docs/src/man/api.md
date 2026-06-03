```@meta
CurrentModule = OpenKolmogorovFlow
```

# Reference

## Indexing

```@docs
WaveNumber
@loop_k
@loop_jk
```

## Fields

```@docs
AbstractField
Field
make_grid
AbstractFTField
FTField
growto!
```

## FFTs

```@docs
up_dealias_size
down_dealias_size
ForwardFFT!
InverseFFT!
FFT
IFFT
```

## Operators

```@docs
ddx!
ddy!
laplacian!
invlaplacian!
```

## Forward Equation

```@docs
ForwardExplicitTerm
ForwardEquation
splitexim
ImplicitTerm
```

## Linearised Equations

```@docs
AbstractLinearMode
TangentMode
AdjointMode
LinearisedExTerm
LinearisedEquation
```

## Forcings

```@docs
AbstractForcing
DummyForcing
SteadyForcing
WaveNumberForcing
DissRateGradientForcing
ReForcing
```

## Hooks

```@docs
CFLHook
```

## Diagnostics

```@docs
laminarflow
dissrate
powinput
radial_mean!
radial_mean
```

## Norms

```@docs
normdiff
minnormdiff
```

## Symmetries

```@docs
xshift!
yshift!
yshiftreflect!
rotate!
shift!
shiftrotate!
```
