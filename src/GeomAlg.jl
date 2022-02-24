module GeomAlg

import Cassette
using Reexport, StaticArrays
using Base.Cartesian: @nexprs, @ncall

include("util.jl")
using .Utils

include("QuadraticForms.jl")
@reexport using .QuadraticForms

include("interface.jl")
export
    AbstractMultivector, AbstractBlade, vectorspacedim, scalarfieldtype,
    grade, rightbasisdual, leftbasisdual, basiselem,
    similarmv, similarmvtype

abstract type AbstractMultivector{K, N} end
const AM = AbstractMultivector
abstract type AbstractBlade{K, N} <: AM{K, N} end

scalarfieldtype(T::Type) =
    error("scalarfieldtype is not implemented for type $T")
scalarfieldtype(x) = scalarfieldtype(typeof(x))
scalarfieldtype(::Type{<:AbstractQuadraticForm{K}}) where K = K
scalarfieldtype(MV::Type{<:AM{K}}) where K = K

vectorspacedim(T::Type) = error("vectorspacedim is not implemented for type $T")
vectorspacedim(x) = vectorspacedim(typeof(x))
vectorspacedim(::Type{<:AbstractQuadraticForm{<:Any, N}}) where N = N
vectorspacedim(MV::Type{<:AM{<:Any, N}}) where N = N

Base.getindex(x::AM, g::Int) = grade(x, g)

struct MultivectorBasis{MV<:AM} end
MultivectorBasis(MV::Type{<:AM}) = MultivectorBasis{MV}()
MultivectorBasis(x::AM) = MultivectorBasis(typeof(x))
multivectortype(::Type{MultivectorBasis{MV}}) where MV<:AM = MV
multivectortype(x::MultivectorBasis) = multivectortype(typeof(x))
scalarfieldtype(MB::Type{<:MultivectorBasis}) =
    scalarfieldtype(multivectortype(MV))
vectorspacedim(MB::Type{<:MultivectorBasis}) =
    vectorspacedim(multivectortype(MV))

let CACHE=Dict{Tuple{DataType, Vector{Int}}, AM}()
function Base.getindex(x::MultivectorBasis, idxs::NTuple{<:Any, Int})
    MV = multivectortype(x)
    key = (MV, collect(idxs))
    if haskey(CACHE, key)
        CACHE[key]
    else
        CACHE[key] = basiselem(MV, idxs...)
    end
end
end
Base.getindex(x::MultivectorBasis, idxs::Int...) = x[idxs]
Base.getindex(x::MultivectorBasis, idxs...) = x[to_indices(x, idxs)]

function Base.show(io::IO, mime::MIME"text/plain", x::MultivectorBasis)
    print(io, "MultivectorBasis for ")
    show(io, mime, multivectortype(x))

    nothing
end

quadraticform(a::AM, b::AM) =
    quadraticform(promote_type(typeof(a), typeof(b)))

Cassette.@context QFormCtx
function Cassette.overdub(
    ctx::QFormCtx{QF}, ::typeof(quadraticform),
    x::AM{K,N}, y::AM{K,N}
) where {K, N, QF<:AM{K,N}}
    if hasmethod(quadraticform, typeof(x))
        quadraticform(xs...)
    elseif hasmethod(quadraticform, Tuple{Type{K}, typeof(N)})
        quadraticform(K, N)
    else
        ctx.metadata
    end
end
function Cassette.overdub(
    ctx::QFormCtx{QF}, ::typeof(quadraticform), x::AM{K,N}
) where {K, N, QF<:AbstractQuadraticForm{K,N}}
    if hasmethod(quadraticform, Tuple{typeof(x)})
        quadraticform(x)
    elseif hasmethod(quadraticform, Tuple{Type{K}, typeof(N)})
        quadraticform(K, N)
    else
        ctx.metadata
    end
end
function Cassette.overdub(
    ctx::QFormCtx, ::typeof(quadraticform), xs...
)
    if hasmethod(quadraticform, typeof(x))
        quadraticform(xs...)
    else
        error("No quadraticform associated with $(typeof.(xs))")
    end
end

macro with_quadratic_form(form, body)
    :(Cassette.recurse(QFormCtx(; metadata=$(esc(form))), () -> $(esc(body))))
end

include("interface.jl")
include("operators.jl")
include("conversion.jl")
include("basisindexiter.jl")
#include("basisindexiter_reference.jl")

basiselem(MV::Type{<:AM}) = one(MV)

basiscoeff(x::AM, I::CartesianIndex) = basiscoeff(x, Tuple(I))
basiscoeff(x::AM, idxs::Int...) = basiscoeff(x, idxs)

hasgrade(x::AM) = !iszero(x[0])

gradeinv(x::AM) =
    sum(iseven(i) ? x[i] : -x[i] for i = 0:vectorspacedim(x))

notrevgrade(i) = iszero(i & oftype(i, 2))
rev(x::AM) =
    sum(revgrade(i) ? x[i] : -x[i] for i = 0:vectorspacedim(x))

### Requires a notion of sign over the base field.
### Assumes a geometric product basis.
metricinv(x::AM) = metricinv(quadraticform(x), x)
function metricinv(q::AbstractQuadraticForm{K,N}, x::AM{K,N}) where {K,N}
    sum(eachbasisindex(N)) do I
        revsq = prod(q(i) for i in Tuple(I) if !iszero(q(i)))
        x / sign(notrevgrade(length(I)) ? revsq : -revsq)
    end
end

graderev(x) = rev(gradeinv(x))
metricrev(x) = rev(metricinv(x))
metricgraderev(x) = rev(metricinv(gradeinv(x)))

include("TreeMV.jl")

function Base.promote_rule(
    MV1::Type{<:AM}, MV2::Type{<:AM}
)
    K1 = scalarfieldtype(MV1)
    K2 = scalarfieldtype(MV2)
    N1 = vectorspacedim(MV1)
    N2 = vectorspacedim(MV2)
    TreeMV.TreeMultivector{promote_type(K1, K2), max(N1, N2), UInt}
end
Base.promote_rule(MV::Type{<:AM}, T::Type) = similarmvtype(MV, T)

end # module GeomAlg
