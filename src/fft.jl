import FFTW: unsafe_execute!, plan_rfft, plan_brfft

export FFT, IFFT, ForwardFFT!, InverseFFT!, up_dealias_size, down_dealias_size

# ~~~ UTILS ~~~
# Set coefficients outside the active truncation to zero. This is the mask that
# keeps dealiased storage modes from feeding back into the evolved active set.
function _apply_mask(U::AbstractFTField{n, m, T}) where {n, m, T}
    @inbounds begin
        # middle block
        for jj = 1:n+1
            @simd for kk = (n+2):2m+2-n
                U.data[kk, jj] = zero(Complex{T})
            end
        end
        # vertical block
        for jj = n+2:m+2
            @simd for kk = 1:2m+2
                U.data[kk, jj] = zero(Complex{T})
            end
        end
        # zero mean
        U.data[1, 1] = 0
    end
    return U
end

# Enforce the Hermitian symmetry required on the rfft zero column. The rest of
# the Hermitian half-plane is represented implicitly by FFTW's real transform.
function _apply_symmetry(U::AbstractFTField{n, m, T}) where {n, m, T}
    @inbounds @simd for k = 1:m
        pos = U.data[k+1, 1]
        neg = U.data[2*(m+1)-k+1, 1]
        _re = 0.5 * (real(pos) + real(neg))
        _im = 0.5 * (imag(pos) - imag(neg))
        U.data[k+1, 1] = _re + im * _im 
        U.data[2*(m+1)-k+1, 1] = _re - im * _im 
    end
    return U
end

"""
    up_dealias_size(n::Int) -> Int

Return the storage cutoff `m` used for the 3/2-rule dealiased grid associated
with active cutoff `n`.
"""
up_dealias_size(n::Int) = n + n>>1

"""
    down_dealias_size(m::Int) -> Int

Return the largest active cutoff `n` whose 3/2-rule storage cutoff fits inside
an existing physical-space storage cutoff `m`.
"""
down_dealias_size(m::Int) = findlast(n->(up_dealias_size(n) ≤ m), 1:m)

# ~~~ NON ALLOCATING VERSION ~~~

"""
    ForwardFFT!(u::AbstractField; flags=FFTW.EXHAUSTIVE)

Reusable in-place forward FFT plan from [`Field`](@ref) to [`FTField`](@ref).

Calling the plan writes Fourier coefficients into the supplied `FTField`,
normalises by the number of physical grid points, applies the active-mode mask,
and restores the stored Hermitian symmetry on the zero `j` column.
"""
struct ForwardFFT!{m, P}
    plan::P
    function ForwardFFT!(u::AbstractField{m}, flags=FFTW.EXHAUSTIVE) where {m}
        plan = plan_rfft(parent(u), [2, 1], flags=flags)
        new{m, typeof(plan)}(plan)
    end
end

# callable interface
(f::ForwardFFT!{m})(U::FTField{n, m}, u::Field{m}) where {n, m} =
    (unsafe_execute!(f.plan, parent(u), parent(U)); 
        U .*= 1/(2m+2)^2; _apply_symmetry(_apply_mask(U)))


"""
    InverseFFT!(U::AbstractFTField; flags=FFTW.EXHAUSTIVE)

Reusable in-place inverse FFT plan from [`FTField`](@ref) to [`Field`](@ref).

The input Fourier field is masked before FFTW executes. This is appropriate for
solver caches where inactive modes are scratch storage; callers that need to
preserve every stored coefficient should pass a copy.
"""
struct InverseFFT!{m, P}
    plan::P
    function InverseFFT!(U::AbstractFTField{n, m}, flags=FFTW.EXHAUSTIVE) where {n, m}
        plan = plan_brfft(parent(U), 2m+2, [2, 1], flags=flags)
        new{m, typeof(plan)}(plan)
    end
end

# callable interface
(i::InverseFFT!{m})(u::Field{m}, U::FTField{n, m}) where {n, m} =
    (unsafe_execute!(i.plan, parent(_apply_mask(U)), parent(u)); u)


"""
    FFT(u::AbstractField, n::Int) -> FTField

Allocate a Fourier field with active cutoff `n` and transform the physical
field `u` into it.

The storage cutoff is inherited from `u`, so `u::Field{m}` produces
`FTField{n,m}`.
"""
function FFT end

# We need copies because the plan destroys the input
function FFT(u::AbstractField{m, T}, n::Int) where {m, T}
    v = copy(u)
    fun = ForwardFFT!(v, FFTW.ESTIMATE); v .= u
    return fun(FTField(n, m, T), v)
 end

"""
    IFFT(U::AbstractFTField) -> Field

Allocate a physical-space field and inverse transform `U` into it.

The returned field has the storage cutoff `m` of `U`.
"""
function IFFT(U::AbstractFTField{n, m, T}) where {n, m, T}
    V = copy(U)
    fun = InverseFFT!(V, FFTW.ESTIMATE); V .= U
    return fun(Field(m, T), V)
end
