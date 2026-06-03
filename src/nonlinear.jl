import Base.Threads: @spawn, threadid, nthreads
using FFTW

export ForwardEquation, splitexim

# ~~~ Explicit Term of Forward Equations ~~~
"""
    ForwardExplicitTerm(n, m, kforcing, flags, T)

Preallocated explicit nonlinear term for the forward vorticity equation.

This object owns Fourier and physical work arrays, plus reusable FFT plans. It
computes advection, updates the CFL estimate `β`, and adds the built-in
Kolmogorov forcing modes.
"""
struct ForwardExplicitTerm{n, m, FT, F, ITT, FTT}
     FTCache::Vector{FT} # storage
      FCache::Vector{F}
       ifft!::ITT        # transforms
        fft!::FTT
    kforcing::Int        # forcing wave number
           β::Vector{Float64}    # used for time step selection based on CFL number
end

# Outer constructor
function ForwardExplicitTerm(n::Int,
                             m::Int,
                             kforcing::Int,
                             flags::UInt32,
                             ::Type{S}) where {S}
    # stupidity check
    m ≥ n || throw(ArgumentError("`m` must be bigger than `n`, got `n, m= $n, $m`. Are you sure?"))

    # real fields have larger size `m`
    FTCache = [FTField(n, m, S) for i = 1:4]
    FCache  = [Field(m, S)      for i = 1:4]

    # transforms
     fft! = ForwardFFT!( FCache[1], flags)
    ifft! = InverseFFT!(FTCache[1], flags)

    # construct object. Initialise β to a small value, so we do a small first time step
    ForwardExplicitTerm{n, m, eltype(FTCache), eltype(FCache),
                 typeof(ifft!), typeof(fft!)}(FTCache, FCache, ifft!, fft!, kforcing, [1e-6])
end

# Adhere to callable interface for IMEXRKCB
function (Eq::ForwardExplicitTerm{n, m, FT})(t::Real,
                                             Ω::FT,
                                             dΩdt::FT,
                                             add::Bool=false) where {n, m, FT<:AbstractFTField{n, m}}
    # extract aliases
    dΩdx, dΩdy, U, V = Eq.FTCache
    dωdx, dωdy, u, v = Eq.FCache

    # Work entirely spectrally until physical multiplication is unavoidable.
    # These two temporaries become ∂ω/∂x and ∂ω/∂y after inverse transforms.
    ddx!(dΩdx, Ω)
    ddy!(dΩdy, Ω)

    # Recover velocity from vorticity using the streamfunction relation.
    # U = -∂ψ/∂y and V = ∂ψ/∂x when Ω = Δψ.
    invlaplacian!(U, dΩdy); U .*= -1
    invlaplacian!(V, dΩdx)

    # inverse transform to physical space into temporaries
    @sync for i = 1:4
        @spawn begin
            Eq.ifft!(Eq.FCache[i], Eq.FTCache[i])
        end
    end

    # calculate β = Δx/u_max. Use absolute velocity so large negative
    # velocities constrain the CFL step just like positive ones.
    u_max = max(maximum(abs, u), maximum(abs, v))
    Eq.β[1] = iszero(u_max) ? Inf : (π/(m+1))/u_max

    # Multiply in physical space. Reuse `u` as the nonlinear RHS buffer:
    # -(u ∂xω + v ∂yω).
    u  .= .- u.*dωdx .- v.*dωdy

    # forward transform to Fourier space into destination
    add == true ? (Eq.fft!(U, u); dΩdt .+= U) : Eq.fft!(dΩdt, u)

    # ~~~ FORCING TERM ~~~
    dΩdt[WaveNumber( Eq.kforcing, 0)] -= Eq.kforcing/2
    dΩdt[WaveNumber(-Eq.kforcing, 0)] -= Eq.kforcing/2

    return nothing
end


# ~~~ SOLVER OBJECT FOR THE GOVERNING EQUATIONS ~~~
"""
    ForwardEquation(n, m, Re[, kforcing=4, flags=FFTW.EXHAUSTIVE,
                    forcing=DummyForcing(n), T=Float64])

Forward vorticity equation for open Kolmogorov flow.

The equation combines the implicit viscous term, the explicit nonlinear
advection and built-in Kolmogorov forcing, and an optional user-supplied
forcing object. It follows the callable interface expected by `Flows.jl`.
"""
struct ForwardEquation{n, m, 
                       IT<:ImplicitTerm, 
                       ET<:ForwardExplicitTerm,
                       GT<:AbstractForcing{n}}
    imTerm::IT
    exTerm::ET
    forcing::GT  # extra forcing

end

# outer constructor: main entry point
function ForwardEquation(n::Int,
                         m::Int,
                         Re::Real,
                         kforcing::Int=4,
                         flags::UInt32=FFTW.EXHAUSTIVE,
                         forcing::AbstractForcing=DummyForcing(n),
                         ::Type{S}=Float64) where {S<:Real}
    imTerm = ImplicitTerm(Re)
    exTerm = ForwardExplicitTerm(n, m, kforcing, flags, S)
    ForwardEquation{n, m, 
                    typeof(imTerm), 
                    typeof(exTerm), 
                    typeof(forcing)}(imTerm, exTerm, forcing)
end

# evaluate right hand side of governing equations
(eq::ForwardEquation{n, m})(t::Real,
                            Ω::FTField{n, m},
                            dΩdt::FTField{n, m}) where {n, m} =
    (mul!(dΩdt, eq.imTerm, Ω);
      eq.exTerm(t, Ω, dΩdt, true);
      eq.forcing(t, Ω, dΩdt)) # note forcing always adds to dVdt
  

# /// SPLIT EXPLICIT AND IMPLICIT PARTS ///
"""
    splitexim(eq::ForwardEquation) -> (explicit, implicit)

Return the explicit and implicit pieces of a [`ForwardEquation`](@ref).

The returned `explicit` closure has the signature expected by IMEX integrators
in `Flows.jl`; the returned implicit term is the equation's [`ImplicitTerm`](@ref).
"""
function splitexim(eq::ForwardEquation{n, m}) where {n, m}
    function wrapper(t::Real,
                     Ω::FTField{n, m},
                  dΩdt::FTField{n, m},
                   add::Bool=false)
        eq.exTerm(t, Ω, dΩdt, add)
        eq.forcing(t, Ω, dΩdt) # note forcing always adds to dΛdt
        return dΩdt
    end

    return wrapper, eq.imTerm
end
