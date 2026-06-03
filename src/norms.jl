import LinearAlgebra: dot, norm

# Definitions of inner products, norms and distances
export normdiff, minnormdiff

"""
    dot(Ω₁, Ω₂) -> Real

Fourier-space inner product between two vorticity fields.

The stored `rfft` half-plane is weighted so modes with implicit conjugate
partners contribute twice, while the `j == 0` column contributes once. The mean
mode correction avoids double counting `(0, 0)`.
"""
# Inner product between two vorticity fields (Ω₁, Ω₂) with possible non-zero mean
function dot(Ω₁::AbstractFTField{n, m, T}, Ω₂::AbstractFTField{n, m, T}) where {n, m, T}
    # we split the loop because different part carry different weight
    s = zero(T)
    @inbounds for _j = 1:n
        @loop_k n m (s += 2 * real(Ω₁[_k, _j] * conj(Ω₂[_k, _j])))
    end
    @loop_k n m (s += real(Ω₁[_k, 0] * conj(Ω₂[_k, 0])))
    # do not count the mean flow twice
    s -= 0.5*real(Ω₁[0, 0] * conj(Ω₂[0, 0]))
    return s
end

"""
    normdiff(Ω₁, Ω₂) -> Real

Return the squared norm of `Ω₁ - Ω₂` using the same Fourier weights as
the Fourier-space `dot` method.
"""
function normdiff end

# Inner product of the difference of two vorticity fields (Ω₁-Ω₂, Ω₁-Ω₂)
function normdiff(Ω₁::AbstractFTField{n, m, T}, Ω₂::AbstractFTField{n, m, T}) where {n, m, T}
    s = zero(T)
    @inbounds for _j = 1:n
        @loop_k n m (s += 2*abs2(Ω₁[_k, _j] - Ω₂[_k, _j]))
    end
    @loop_k n m (s += abs2(Ω₁[_k, 0] - Ω₂[_k, 0]))
    s -= 0.5*abs2(Ω₁[0, 0] - Ω₂[0, 0])
    return s
end

"""
    norm(Ω, [p=2]) -> Real

Two-norm of a Fourier-space vorticity field.

Only `p == 2` is supported.
"""
# The two norm of a field ||Ω|| = sqrt(inner(Ω, Ω))
norm(Ω::AbstractFTField, n::Int=2) =
    (n==2 || throw(ArgumentError("only the 2-norm is defined"));
    sqrt(dot(Ω, Ω)))


"""
    minnormdiff(Ω₁, Ω₂[, TMP, N]) -> (distance, (xshift, yshift_index))

Scan over package symmetries and return the smallest squared distance between
shifted copies of `Ω₁` and `Ω₂`.

`N` controls the number of trial shifts in the continuous `x` direction. The
integer `yshift_index` is passed to [`yshift!`](@ref), corresponding to quarter
period shifts of the default `kforcing = 4` symmetry.
"""
function minnormdiff end

# return minimum distance across shifts
function minnormdiff(Ω₁::FTField{n, m},
                     Ω₂::FTField{n, m},
                     TMP::FTField{n, m} = similar(Ω₁),
                     N::Int = 15) where {n, m}

    # mimimum values 
    dmin = normdiff(Ω₁, Ω₂)
    smin = 0.0
    mmin = 0

    # for every y shift
    for _m = (0, 1, 2, 3)
        TMP .= Ω₁
        yshift!(TMP, _m)

        # scan through x shifts
        for k = 1:N-1
            xshift!(TMP, 2*π/N)
            d = normdiff(TMP, Ω₂)
            if d < dmin
                dmin = d
                smin = 2*π/N * k
                mmin = _m
            end
        end
    end
   
    return dmin, (smin, mmin)
end
