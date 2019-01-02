using FFTW

export ForwardEquation, splitexim

# ~~~ Explicit Term of Forward Equations ~~~
struct ForwardExplicitTerm{n, m, FT<:AbstractFTField,   F<:AbstractField,
                                ITT<:InverseFFT!,     FTT<:ForwardFFT!}
     FTCache::Vector{FT} # storage
      FCache::Vector{F}
       ifft!::ITT        # transforms
        fft!::FTT
    kforcing::Int        # forcing wave number
           β::Vector{Float64}    # used for time step selection based on CFL number
end

# Outer constructor: TODO: merge with above
function ForwardExplicitTerm(n::Int, m::Int, kforcing::Int, ::Type{S}) where {S}
    # stupidity check
    m ≥ n || throw(ArgumentError("`m` must be bigger than `n`, got `n, m= $n, $m`. Are you sure?"))

    # real fields have larger size `m`
    FTCache = [FTField(n, m, S) for i = 1:4]
    FCache  = [Field(m, S)      for i = 1:4]

    # transforms
    ifft! = InverseFFT!(FTCache[1])
     fft! = ForwardFFT!(FCache[1])

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

    # obtain vorticity derivatives
    ddx!(dΩdx, Ω)
    ddy!(dΩdy, Ω)

    # obtain velocity components
    invlaplacian!(U, dΩdy); U .*= -1
    invlaplacian!(V, dΩdx)

    # inverse transform to physical space into temporaries
    Eq.ifft!(u,    U)
    Eq.ifft!(v,    V)
    Eq.ifft!(dωdx, dΩdx)
    Eq.ifft!(dωdy, dΩdy)

    # calculate β = Δx/u_max
    Eq.β[1] = (π/(m+1))/max(maximum(u), maximum(v))

    # multiply in physical space. Overwrite u
    u  .= .- u.*dωdx .- v.*dωdy

    # forward transform to Fourier space into destination
    add == true ? (Eq.fft!(U, u); dΩdt .+= U) : Eq.fft!(dΩdt, u)

    # ~~~ FORCING TERM ~~~
    dΩdt[WaveNumber( Eq.kforcing, 0)] -= Eq.kforcing/2
    dΩdt[WaveNumber(-Eq.kforcing, 0)] -= Eq.kforcing/2

    return nothing
end


# ~~~ SOLVER OBJECT FOR THE GOVERNING EQUATIONS ~~~
struct ForwardEquation{n, m, IT<:ImplicitTerm, ET<:ForwardExplicitTerm}
    imTerm::IT
    exTerm::ET
end

# outer constructor: main entry point
function ForwardEquation(n::Int,
                         m::Int,
                         Re::Real,
                         kforcing::Int=4,
                         ::Type{S}=Float64) where {S<:Real}
    imTerm = ImplicitTerm(Re)
    exTerm = ForwardExplicitTerm(n, m, kforcing, S)
    ForwardEquation{n, m, typeof(imTerm), typeof(exTerm)}(imTerm, exTerm)
end

# evaluate right hand side of governing equations
(eq::ForwardEquation{n, m})(t::Real,
                            Ω::FTField{n, m},
                            dΩdt::FTField{n, m}) where {n, m} =
    (A_mul_B!(dΩdt, eq.imTerm, Ω); eq.exTerm(t, Ω, dΩdt, true))

# obtain two components
splitexim(eq::ForwardEquation) = (eq.exTerm, eq.imTerm)