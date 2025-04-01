export DistanceCache, distance!

# Storage type to assist the computation of distance between two vorticity field
struct DistanceCache{nc, T, F<:Field{nc, T}, FT<:FTField{nc, nc, T}}
         r::F                    # the cross correlation field in physical space
    r_line::Vector{T}            # cross correlation along Δy = nc*π/4
         R::FT                   # the cross correlation in Fourier space
    R_line::Vector{Complex{T}}   # Fourier series or r_line
       fft                       # 2d fft
      plan                       # 1d fft plan
    function DistanceCache(nc::Int, ::Type{T}=Float64, flags=FFTW.EXHAUSTIVE) where {T}
        # number of points in the physical grid
        M = 2*nc + 2
        
        # we want to index exact location for vertical integer
        # shifts, so grid size must be divisible by 4
        # ! this condition contradicts the condition that nc should be even
        rem(M, 4) == 0 || error("cache size must be divisible by four")
        
        # store cross correlation of the two velocity fields
        # in physical and Fourier space
        r, R = Field(nc, T), FTField(nc, nc, T)

        # cross correlation along a line
        r_line, R_line = zeros(T, M), zeros(Complex{T}, nc+2)

        # define plans and transforms
         fft = ForwardFFT!(r, flags)
        plan = plan_rfft(r_line, flags=flags) # forward transform
        new{nc, T, typeof(r), typeof(R)}(r, r_line, R, R_line, fft, plan)
    end
end

"""
    distance!(Ω₁, Ω₂, cache)

Calculate the minimum distance (L2 norm) between vorticity fields `Ω₁` and `Ω₂`,
among all possible continuous shift `s ∈ [0, 2π)` along the `x` direction and
for discrete vertical shifts `m*π/4, m ∈ {0, 2, 4, 6}`. It returns the minimum
distance `d` and the shifts `s` and `m`, such that

julia> dotdiff(shift!(Ω₁, (s, m)), Ω₂) ≈ d

where the equality holds approximately because different algorithms are used.
This function exploits the correlation theorem to compute the optimal shift
more efficiently.

We currently use the algorithm that is fastest and allocates less memory, but
can be significantly less precise to to cancellation errors. This is not a big
issue, because the main purpose of this function is to detect near recurrences
and we do not care if the distance is not calculated to 16 decimal digits, but
only, say, 11.

The object `cache` is created as

julia> cache = DistanceCache(nc)

where `nc` is the size of the grid where the correlation function is evaluated
on. This can be the same size as the fields `Ω₁` and `Ω₂` above, or lower, for
faster, but maybe less accurate calculations. The number `nc/4` must be even.
"""
function distance!(Ω₁::FTField{n, m, T}, Ω₂::FTField{n, m, T}, cache::DistanceCache{nc, T}) where {n, m, nc, T}
    nc ≤ n || throw(ArgumentError("cache size must be lower or equal that field size"))
    
    # evaluate product in fourier space
    _evalconjprod!(Ω₁, Ω₂, cache.R)

    # transform cross correlation to physical space
    cache.fft(cache.R, cache.r)

    # the find correlation peak, where distance is minimum
    s_opt, m_opt, r_opt = _peakdetect(cache)

    # return minimum distance and the shifts to obtain it
    return norm(Ω₁)^2 + norm(Ω₂)^2 - 8*π^2*r_opt, (s_opt, m_opt)
end


# ~~~~ Below is the private interface, mainly helper functions ~~~

# Evaluate product of two fields in Fourier space, 
# obtaining the cross correlation R in Fourier space
function _evalconjprod!(Ω₁::FTField,
                        Ω₂::FTField,
                         R::FTField{nc}) where {nc}
    R .= 0
    @inbounds begin
        for j = 0:nc, k = -nc:nc
            R[WaveNumber(k, j)] = conj(Ω₁[WaveNumber(k, j)])*Ω₂[WaveNumber(k, j)]
        end
    end
    R[WaveNumber(0, 0)] = 0
    
    return nothing
end

# Find peak of correlation function for y shifts m = 0, 1, 2 and 3.
function _peakdetect(cache::DistanceCache{nc}) where {nc}
    # scan correlation data to find maximum for discrete vertical shifts
    M  = div(2*nc+2, 4)
    (r_opt_0, s_opt_int_0) = findmax(view(cache.r, 0*M + 1, :))
    (r_opt_2, s_opt_int_2) = findmax(view(cache.r, 1*M + 1, :))
    (r_opt_4, s_opt_int_4) = findmax(view(cache.r, 2*M + 1, :))
    (r_opt_6, s_opt_int_6) = findmax(view(cache.r, 3*M + 1, :))

    # this calculates the maximum by r_opt
    ((r_opt, m_opt), s_opt_int) = findmax(((r_opt_0, 0, s_opt_int_0),
                                           (r_opt_2, 1, s_opt_int_2),
                                           (r_opt_4, 2, s_opt_int_4),
                                           (r_opt_6, 3, s_opt_int_6)))

    # Copy value of correlation function at `Δy = m_opt*π/4` into a temporary
    # buffer and compute its Fourier transform, so that we can use it for
    # interpolation and for finding the correlation peak.
    cache.r_line .= view(cache.r, m_opt*M + 1, :)
    unsafe_execute!(cache.plan, cache.r_line, cache.R_line)
    cache.R_line .*= 1/length(cache.r_line)

    # use newton method to find correlation peak accurately
    s_opt, r_max = _peroptim(cache.R_line, 2π/(2*nc+2)*s_opt_int)

    # return x shift, y shift and correlation maximum
    return s_opt, m_opt, r_max
end

# Use Newton method to find an extremum of the periodic function defined
# by the Fourier Series `R` near `x`. We assume that this is going to
# converge in less than 10 iteration, otherwise we raise an error.
function _peroptim(R::Vector{<:Complex}, x::Real, stol::Real=1e-12)
    for i = 1:10 # this must end
        val, grad, hess = _perinterp(R, x)
        Δ = grad/hess
        x = x - Δ
        abs(Δ) < stol && return (x, _perinterp(R, x)[1])
    end
    throw(error("iterations did not converge"))
end

# Interpolate periodic function at `x` from its Fourier series `R`. Returns
# value and first and second derivatives at `x`. Only works for even sized
# input arrays (cover the use cases in this package). We assume the grid
# is over the domain [0, 2π].
function _perinterp(R::Vector{Complex{T}}, x::Real) where {T}
    lX = length(R)
    lX % 2 == 0 || throw(ArgumentError("array length must be even"))
    @inbounds begin
        val, grad, hess  = zero(T), zero(T), zero(T)
        for kk = 2:lX-1
            coskx, sinkx = cos((kk-1)*x), sin((kk-1)*x)
            rexk, imxk = reim(R[kk])
            val  +=  (rexk*coskx - imxk*sinkx)
            grad -=  (rexk*sinkx + imxk*coskx)*(kk-1)
            hess -=  (rexk*coskx - imxk*sinkx)*(kk-1)^2
        end
        val  = 2*val  + real(R[lX])*cos((lX-1)*x) + real(R[1])
        grad = 2*grad - real(R[lX])*sin((lX-1)*x)*(lX-1)
        hess = 2*hess - real(R[lX])*cos((lX-1)*x)*(lX-1)^2
    end
    val, grad, hess
end