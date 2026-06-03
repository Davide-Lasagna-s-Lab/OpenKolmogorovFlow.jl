export AbstractFTField,
       FTField,
       growto!

# ~~~ ABSTRACT TYPE FOR STRUCTURED, ARRAY-LIKE OBJECTS FROM FFT ~~~
# The parameter n denotes the largest wave number that we allow
# indexing on a grid with number of points equal to 2n+2. We do
# not consider the Nyquist component, but always set it to zero
# since it is just a hassle to take it into account properly. The
# parameter m is the actual largest wave number that can be
# represented by the underlying data.
abstract type AbstractFTField{n, m, T<:Real} <: AbstractMatrix{Complex{T}} end

Base.axes(::AbstractFTField{n, m}) where {n, m} = (0:2*m+1, 0:m+1)
Base.size(U::AbstractFTField{n, m}) where {n, m} = (2*m+2, m+2)
Base.IndexStyle(::Type{<:AbstractFTField}) = IndexLinear()

# ~~~ MAIN TYPE ~~~
struct FTField{n, m, T<:Real, M<:AbstractMatrix{Complex{T}}} <: AbstractFTField{n, m, T}
    data::M
    # pass field sizes and data type
    function FTField{n, m, T}() where {n, m, T}
        m ≥ n || throw(ArgumentError("m can't be smaller than n"))
        data = zeros(Complex{T}, 2m+2, m+2)
        new{n, m, T, typeof(data)}(data)
    end
    # pass the data over the larger grid, and the actual largest wave number
    function FTField(input::AbstractMatrix{Complex{T}}, n::Int) where {T<:Real}
        _checksize(input)
        # get largest wave number
        m = size(input, 1) >> 1 - 1
        m ≥ n || throw(ArgumentError("n can't be larger than m"))
        # instantiate object
        new{n, m, T, typeof(input)}(input)
    end
end

# OUTER CONSTRUCTORS
FTField(n::Int, m::Int, ::Type{T}=Float64) where {T} = FTField{n, m, T}()

# size checker
function _checksize(data::AbstractMatrix{<:Complex})
    M, N = size(data)
    iseven(M)   || throw(ArgumentError("input row size must be even"))
    M>>1 == N-1 || throw(ArgumentError("inconsistent input size"))
    return nothing
end

# this is for JFOp in NKSearch
Base.hash(U::FTField) = hash(parent(U))

# ~~~ array interface ~~~
Base.@propagate_inbounds function Base.getindex(U::FTField{n, m},
                                               _k::Int,
                                               _j::Int) where {n, m}
    @boundscheck checkbounds(U, _k, _j)
    @inbounds ret = U.data[_k+1, _j+1]
    return ret
end

Base.@propagate_inbounds function Base.setindex!(U::FTField{n, m},
                                               val::Number,
                                                _k::Int,
                                                _j::Int) where {n, m}
    @boundscheck checkbounds(U, _k, _j)
    @inbounds U.data[_k+1, _j+1] = val
    return val
end

# Indexing with wavenumber
Base.@propagate_inbounds function Base.getindex(U::FTField{n, m},
                                                w::WaveNumber) where {n, m}
    _k, _j = _reindex(w.k, w.j, m)
    @boundscheck checkbounds(U, _k, _j)
    @inbounds ret = U.data[_k+1, _j+1]
    return ret
end

Base.@propagate_inbounds function Base.setindex!(U::FTField{n, m},
                                               val::Number,
                                                 w::WaveNumber) where {n, m}
    _k, _j = _reindex(w.k, w.j, m)
    @boundscheck checkbounds(U, _k, _j)
    @inbounds U.data[_k+1, _j+1] = val
    return val
end

# Linear indexing
Base.@propagate_inbounds Base.getindex(U::FTField, i::Int) =
    (@boundscheck checkbounds(U.data, i);
     @inbounds ret = U.data[i]; ret)

 Base.@propagate_inbounds Base.setindex!(U::FTField, val::Number, i::Int) =
    (@boundscheck checkbounds(U.data, i);
     @inbounds U.data[i] = val; val)

# accessors functions
Base.parent(U::FTField) = U.data

# allow constructing similar fields
Base.similar(U::FTField{n, m, T}) where {n, m, T} = FTField(n, m, T)
Base.zero(U::FTField) = similar(U)
Base.copy(U::FTField) = (V = similar(U); V .= U; V)

# ~~~ interpolate field to larger grid ~~~

# same size is just a copy
@inline growto!(OUT::FTField{n, m}, U::FTField{n, m}) where {n, m} = OUT .= U

# different size
@inline function growto!(OUT::AbstractFTField{n, m}, U::AbstractFTField{p, q}) where {n, m, p, q}
    m >= q ||
        throw(ArgumentError("output storage size must be larger than or equal to input storage size"))
    n >= p ||
        throw(ArgumentError("output active wavenumber cutoff must be larger than or equal to input active cutoff"))
    @inbounds begin
        OUT .= 0
        for j = 0:p, k = -p:p
            OUT[WaveNumber(k, j)] = U[WaveNumber(k, j)]
        end
    end
    OUT
end
