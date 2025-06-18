export shift!, shiftrotate!, xshift!, yshift!, yshiftreflect!, rotate!

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

rotate!(U::FTField) = map!(conj, U, U)

shift!(U::FTField, s::Real, m::Int) = yshift!(xshift!(U, s), m)
shiftrotate!(U::FTField, s::Real, m::Int) = rotate!(yshift!(xshift!(U, s), m))
