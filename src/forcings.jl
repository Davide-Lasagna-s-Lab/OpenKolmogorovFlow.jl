# ////// FORCINGS FOR THE DIRECT, TANGENT OR ADJOINT EQUATIONS//////

export DummyForcing,
       DissRateGradientForcing,
       WaveNumberForcing,
       ReForcing,
       SteadyForcing

"""
    AbstractForcing{n}

Supertype for additive forcing objects acting on fields with active cutoff `n`.

Forcing objects are callable and mutate the right-hand-side field passed as the
last argument. Implementations may support the forward, tangent, or adjoint
equation signatures depending on where the forcing is meaningful.
"""
abstract type AbstractForcing{n} end


# ////// DUMMY FORCING - DOES NOTHING //////
"""
    DummyForcing(n)

No-op forcing with active cutoff `n`.

Use this as the default forcing object when the equation should only include
its built-in terms.
"""
struct DummyForcing{n} <: AbstractForcing{n} end

DummyForcing(n::Int) = DummyForcing{n}()

(::DummyForcing{n})(t::Real,
                    Ω::FT,
                 dΩdt::FT) where {n, FT<:AbstractFTField{n}} = dΩdt

(::DummyForcing{n})(t::Real,
                    Ω::FT,
                    Λ::FT,
                 dΛdt::FT) where {n, FT<:AbstractFTField{n}} = dΛdt

(::DummyForcing{n})(t::Real,
                    Ω::FT,
                 dΩdt::FT,
                    Λ::FT,
                 dΛdt::FT) where {n, FT<:AbstractFTField{n}} = dΛdt

# ////// GENERIC FORCING   //////
"""
    SteadyForcing(f)

Time-independent additive forcing stored as an [`FTField`](@ref).

Calling the forcing adds `f` to the destination right-hand side.
"""
 struct SteadyForcing{n, FT} <: AbstractForcing{n} 
    f::FT
    SteadyForcing(h::FT) where {n, FT<:AbstractFTField{n}} = new{n, FT}(h)
end

# provide multiple methods for the same function
function (f::SteadyForcing{n})(t::Real,
                               Ω::FT,
                            dΩdt::FT) where {n, FT<:AbstractFTField{n}} 
    dΩdt .+= f.f
    return dΩdt
end

# ////// FORCING AT ONE SPECIFIC WAVENUMBER  //////
"""
    WaveNumberForcing(n, k, j, value)

Add a forcing at one Fourier mode.

For `j == 0`, the conjugate `-k` mode is also updated so the physical-space
forcing remains real. The forcing object can be used with direct, tangent, and
adjoint equation call signatures.
"""
struct WaveNumberForcing{n, T} <: AbstractForcing{n} 
      k::Int
      j::Int
    val::T
end

WaveNumberForcing(n::Int, k::Int, j::Int, f::T) where {T<:Number} = 
    WaveNumberForcing{n, T}(k, j, f)

# provide multiple methods for the same function
(f::WaveNumberForcing{n})(t::Real,
                          Ω::FT,
                       dΩdt::FT) where {n, FT<:AbstractFTField{n}} = f(t, Ω, Ω, dΩdt)

(f::WaveNumberForcing{n})(t::Real,
                          Ω::FT,
                          Λ::FT,
                       dΛdt::FT) where {n, FT<:AbstractFTField{n}} = f(t, Ω, Ω, Λ, dΛdt)
                          
function (f::WaveNumberForcing{n})(t::Real,
                                   Ω::FT,
                                dΩdt::FT,
                                   Λ::FT,
                                dΛdt::FT) where {n, FT<:AbstractFTField{n}}
                 dΛdt[WaveNumber( f.k, f.j)] +=      f.val/2
    f.j == 0 && (dΛdt[WaveNumber(-f.k, f.j)] += conj(f.val)/2)
    return dΛdt
end

# ////// FORCING FOR ADJOINT EQUATION FOR DISSIPATION RATE AS THE FUNCTIONAL //////
"""
    DissRateGradientForcing(n, Re)

Adjoint forcing corresponding to the gradient of dissipation rate.

When called on an adjoint right-hand side, this subtracts `2Ω/Re`.
"""
struct DissRateGradientForcing{n} <: AbstractForcing{n} 
    Re::Float64
end

DissRateGradientForcing(n::Int, Re::Real) = DissRateGradientForcing{n}(Re)

(f::DissRateGradientForcing{n})(t::Real,
                                Ω::FT,
                                Λ::FT,
                             dΛdt::FT) where {n, FT<:AbstractFTField{n}} = 
    (dΛdt .-= 2.0.*Ω./f.Re; dΛdt)

# ////// FORCING FOR TANGENT EQUATION FOR PERTURBATION WRT TO RE //////
"""
    ReForcing(n, m, Re)
    ReForcing(Re, TMP)

Tangent-equation forcing for a perturbation with respect to Reynolds number.

The object owns a temporary Fourier field `TMP` used to apply the Laplacian of
the base vorticity before adding the parameter derivative contribution.
"""
struct ReForcing{n, FT} <: AbstractForcing{n} 
     Re::Float64
    TMP::FT
end

ReForcing(n::Int, m::Int, Re::Real) = ReForcing(Re, FTField(n, m))
ReForcing(Re::Real, TMP::FTField{n}) where {n} = ReForcing{n, typeof(TMP)}(Re, TMP)

(f::ReForcing{n})(t::Real,
                  Ω::FT,
                  Λ::FT,
               dΛdt::FT) where {n, FT<:AbstractFTField{n}} = 
    (dΛdt .-= laplacian!(f.TMP, Ω)./(f.Re.^2); dΛdt)

(f::ReForcing{n})(t::Real,
                  Ω::FT,
               dΩdt::FT,
                  Λ::FT,
               dΛdt::FT) where {n, FT<:AbstractFTField{n}} = f(t, Ω, Λ, dΛdt)
