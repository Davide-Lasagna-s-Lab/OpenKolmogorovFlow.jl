export shift!, shiftrotate!, xshift!, yshift!, yshiftreflect!, rotate!

"""
    xshift!(U, s) -> U

Apply a continuous periodic shift by distance `s` in the `x` direction.

The operation is diagonal in Fourier space and mutates `U` in place.
"""
function xshift!(U::FTField{n}, s::Real) where {n}
    (s == 0) && return U
    @inbounds for j = 0:n
        # precompute this, since it's expensive
        val = cis(j*s)
        @simd for k = -n:n
            U[WaveNumber(k, j)] *= val
        end
    end
    return U
end

"""
    yshift!(U, m) -> U

Apply the discrete `y`-shift symmetry used by the default Kolmogorov forcing.

`m` is an integer quarter-period index; the implementation assumes the forcing
wavenumber `kf = 4`.
"""
function yshift!(U::FTField{n}, m::Int) where {n}
    (m == 0) && return U
    @inbounds for k = -n:n
        # precompute this, since it's expensive. Note also this
        # assumes that kf = 4
        val = cis(k*m*π/2)
        @simd for j = 0:n
            U[WaveNumber(k, j)] *= val
        end
    end
    return U
end

"""
    yshiftreflect!(U) -> U

Apply the package's combined vertical shift/reflection symmetry in place.

This helper assumes the default `kf = 4` symmetry class.
"""
function yshiftreflect!(U::FTField{n}) where {n}
    @inbounds for k = -n:n
        # precompute this, since it's expensive. Note also this
        # assumes that kf = 4
        val = cis(k*π/4)
        @simd for j = 0:n
            U[WaveNumber(k, j)] *= -val
        end
    end
    return U
end

"""
    rotate!(U) -> U

Apply the half-turn/conjugation symmetry in Fourier space.
"""
rotate!(U::FTField) = map!(conj, U, U)

"""
    shift!(U, s, m) -> U

Apply `xshift!(U, s)` followed by `yshift!(U, m)`.
"""
shift!(U::FTField, s::Real, m::Int) = yshift!(xshift!(U, s), m)

"""
    shiftrotate!(U, s, m) -> U

Apply [`shift!`](@ref) followed by [`rotate!`](@ref).
"""
shiftrotate!(U::FTField, s::Real, m::Int) = rotate!(yshift!(xshift!(U, s), m))
