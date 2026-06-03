export laminarflow, dissrate, powinput

"""
    laminarflow(n, m, Re, [kforcing=4]) -> FTField

Return the steady laminar vorticity field for the sinusoidal Kolmogorov
forcing.

The returned field has active cutoff `n`, storage cutoff `m`, and nonzero
coefficients only at `WaveNumber(±kforcing, 0)`. The forcing wavenumber must be
representable in the active mode set.
"""
function laminarflow(n::Int, m::Int, Re::Real, kforcing::Int=4)
    0 ≤ kforcing ≤ n  ||
        throw(ArgumentError("forcing wave number must be in [0, n]"))
    Ω = FTField(n, m)
    Ω[WaveNumber(-kforcing, 0)] = -Re/kforcing/2
    Ω[WaveNumber( kforcing, 0)] = -Re/kforcing/2
    return Ω
end

"""
    dissrate(Ω, Re) -> Real

Energy-dissipation-rate density associated with vorticity field `Ω`.

The value is computed as `norm(Ω)^2 / Re` using the package's Fourier-space
inner product normalization.
"""
dissrate(Ω::FTField{n, m}, Re::Real) where {n, m} = norm(Ω)^2/Re

"""
    powinput(Ω, [kf=4]) -> Real

Return the power input associated with the active forcing mode `kf`.

This diagnostic reads the vorticity coefficient at `WaveNumber(kf, 0)`.
"""
powinput(Ω::FTField, kf::Int=4) = -imag(im * Ω[WaveNumber(kf, 0)]/kf)
