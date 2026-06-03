export Field, make_grid

"""
    AbstractField{m,T} <: AbstractMatrix{T}

Abstract matrix-like representation of a real physical-space field.

The type parameter `m` determines the uniform grid size `2m+2` in both
periodic directions. Concrete fields are stored as ordinary matrices in the
same orientation used by FFTW plans in this package.
"""
abstract type AbstractField{m, T} <: AbstractMatrix{T} end

Base.size(f::AbstractField{m}) where {m} = (2m+2, 2m+2)
Base.IndexStyle(::Type{<:AbstractField}) = Base.IndexLinear()

"""
    Field(data::AbstractMatrix{<:Real})
    Field(m::Int, [T=Float64])
    Field(m::Int, fun)

Physical-space scalar field on the periodic square.

`Field(data)` wraps a square, even-sized real matrix. The grid parameter is
inferred from `size(data,1) == 2m+2`. `Field(m, T)` allocates a zero field, and
`Field(m, fun)` evaluates `fun.(x, y)` on [`make_grid`](@ref).

The object behaves as an `AbstractMatrix`; use `parent(u)` to access the
underlying storage.
"""
struct Field{m, T<:Real, M<:AbstractMatrix{T}} <: AbstractField{m, T}
    data::M
    function Field(data::M) where {T<:Real, M<:AbstractMatrix{T}}
        _checksize(data)
        new{size(data, 1)>>1 - 1, T, M}(data)
    end
end


# OUTER CONSTRUCTORS
Field(m::Int, ::Type{T}=Float64) where {T} = Field(zeros(T, 2m+2, 2m+2))

Field(m::Int, fun::Base.Callable) = Field(fun.(make_grid(m)...))

function _checksize(data::AbstractMatrix{<:Real})
    M, N = size(data)
    N == M    || throw(ArgumentError("input matrix must be square: got $M×$N"))
    iseven(M) || throw(ArgumentError("size must be even"))
    return nothing
end

# ~~~ array interface ~~~
@inline function Base.getindex(f::Field{n}, i::Int, j::Int) where {n}
    @boundscheck checkbounds(f, i, j)
    @inbounds ret = f.data[i, j]
    return ret
end

@inline function Base.setindex!(f::Field{n}, val::Number, i::Int, j::Int) where {n}
    @boundscheck checkbounds(f, i, j)
    @inbounds f.data[i, j] = val
    return val
end

# Linear indexing
@inline function Base.getindex(f::Field, i::Int)
    @boundscheck checkbounds(f, i)
    @inbounds ret = f.data[i]
    return ret
end

@inline function Base.setindex!(f::Field, val::Number, i::Int)
    @boundscheck checkbounds(f, i)
    @inbounds f.data[i] = val
    return val
end

# accessors functions
Base.parent(U::Field) = U.data

Base.similar(u::Field{m, T}) where {m, T} = Field(m, T)
Base.copy(u::Field{m, T}) where {m, T} = (v = similar(u); v .= u; v)

"""
    make_grid(m::Int) -> (x, y)
    make_grid(u::Field) -> (x, y)

Return broadcastable coordinate arrays for the periodic square `[0, 2π)^2`.

Both directions have `2m+2` points. The endpoint `2π` is omitted so periodic
functions are sampled without duplicating the first point.
"""
function make_grid(m::Int)
    x = range(0, stop=2π, length=2m+3)[1:(2m+2)]
    return reshape(x, 1, 2m+2), reshape(x, 2m+2, 1)
end

make_grid(u::Field{m}) where {m} = make_grid(m)
