export ddx!, ddy!, invlaplacian!, laplacian!

"""
    ddx!(OUT, U) -> OUT

Differentiate the Fourier-space field `U` along the `x` direction.

The derivative is applied spectrally by multiplying each active mode by
`im*j`, following the package convention that `j` is the non-negative
real-to-complex direction. `OUT` and `U` must have matching active and storage
resolutions.
"""
function ddx!(OUT::FTField{n,m}, U::FTField{n,m}) where {n, m}
    @loop_jk n m OUT[_k, _j] = im * j * U[_k, _j]
    return OUT
end

"""
    ddy!(OUT, U) -> OUT

Differentiate the Fourier-space field `U` along the `y` direction.

The derivative is applied spectrally by multiplying each active mode by
`im*k`, where `k` is the signed full-FFT wavenumber.
"""
function ddy!(OUT::FTField{n,m}, U::FTField{n,m}) where {n, m}
    @loop_jk n m OUT[_k, _j] = im * k * U[_k, _j]
    return OUT
end

"""
    invlaplacian!(OUT, U) -> OUT
    invlaplacian!(OUT, U, c) -> OUT

Apply a diagonal inverse elliptic operator in Fourier space.

Without `c`, this computes ``-Δ^{-1} U`` mode-by-mode and sets the mean mode to
zero. With `c`, this computes ``(I - cΔ)^{-1} U`` and also zeros the mean mode,
matching the implicit time-stepping operator used by `Flows.ImcA!`.
"""
function invlaplacian!(OUT::FTField{n,m}, U::FTField{n,m}) where {n, m}
    @loop_jk n m OUT[_k, _j] = - U[_k, _j] / (j^2 + k^2)
    @inbounds OUT[WaveNumber(0, 0)] = 0
    return OUT
end

function invlaplacian!(OUT::FTField{n,m}, U::FTField{n,m}, c::Real) where {n, m}
    @loop_jk n m OUT[_k, _j] = U[_k, _j] / (1 + c * (j^2 + k^2))
    @inbounds OUT[WaveNumber(0, 0)] = 0
    return OUT
end

"""
    laplacian!(OUT, U) -> OUT

Apply the spectral Laplacian to `U`.

Each active Fourier coefficient is multiplied by `-(j^2 + k^2)`. The zero
mode is therefore mapped to zero.
"""
function laplacian!(OUT::FTField{n,m}, U::FTField{n,m}) where {n, m}
    @loop_jk n m OUT[_k, _j] = -U[_k, _j] * (j^2 + k^2)
    return OUT
end
