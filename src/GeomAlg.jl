module GeomAlg

import Cassette
using Reexport, StaticArrays
using Base.Cartesian: @nexprs, @ncall

include("util.jl")
export randtmv, randtmv2

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
    ctx::QFormCTX{QF}, ::typeof(quadraticform),
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

∧(args...) = outerprod(args...)
a::AM  +  b::AM = add(a, b)
a::AM  +  b     = add(a, b)
a      +  b::AM = add(a, b)
a::AM  -  b::AM = add(a, b)
a::AM  -  b     = sub(a, b)
a      -  b::AM = sub(a, b)
a::AM  *  b::AM = add(a, b)
a::AM  *  b     = mul(a, b)
a      *  b::AM = mul(a, b)
a::AM  /  b::AM = add(a, b)
a::AM  /  b     = div(a, b)
a      /  b::AM = div(a, b)
a::AM  \  b::AM = add(a, b)
a::AM  \  b     = div(b, a)
a      \  b::AM = div(a, b)
a  ⊣  b = leftcontract(a, b)
a  ⊢  b = rightcontract(a, b)
a  ×₂ b =  a*b - b*a
a  ×  b = (a*b - b*a)/2
a  ∥₂ b =  a*b + b*a
a  ∥  b = (a*b + b*a)/2
mul(a, b) = mul(quadraticform(a, b), a, b)
inv(a) = inv(quadraticform(a), a)
div(a, b) = div(quadraticform(a, b), a, b)
div(q, a, b) = mul(q, a, inv(q, b))
leftcontract(a, b) = leftcontract(quadraticform(a, b), a, b)
rightcontract(a, b) = rightcontract(quadraticform(a, b), a, b)

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

similarmv(x::AM, args...) = convert(similartype(args...), x)
Base.convert(::Type{AM{K}}, x::AM) where K =
    convert(AM{K, vectorspacedim(x)}, x)
function Base.convert(
    ::Type{AM{<:Any, N}}, x::AM
) where N
    convert(AM{scalarfieldtype(x), N}, x)
end

include("basisindexiter.jl")
#include("basisindexiter_reference.jl")

function Base.:(==)(
    x::AM{<:Any,N}, y::AM{<:Any,N}
) where N
    all(basiscoeff(x, I) == basiscoeff(x, I) for I in eachbasisindex(N))
end

Base.:-(x::MV, y::MV) where MV<:AM = x + -y
Base.:*(a::K, x::AM{K}) where K = x*a
Base.:/(x::Union{MV, K}, y::MV) where {K, MV<:AM{K}} = x*inv(y)
for op in (:+, :-, :*, :/) @eval begin
    Base.$op(x::AM, y::AM) =
        $op(promote(x,y)...)
end end
for op in (:+, :-) @eval begin
    Base.$op(x::AM, y) = $op(promote(x,y)...)
    Base.$op(x, y::AM) = $op(promote(x,y)...)
end end
for op in (:*, :/) @eval begin
    Base.$op(x::AM, y) = $op(x, convert(scalarfieldtype(x), y))
    Base.$op(x, y::AM) = $op(convert(scalarfieldtype(y), x), y)
end end

basiselem(MV::Type{<:AM}) = one(MV)

basiscoeff(x::AM, I::CartesianIndex) = basiscoeff(x, Tuple(I))
basiscoeff(x::AM, idxs::Int...) = basiscoeff(x, idxs)

hasgrade(x::AM) = !iszero(x[0])

gradeinv(x::AM) =
    sum(iseven(i) ? x[i] : -x[i] for i = 0:vectorspacedim(x))
function metricinv(x::AbstractMultivector)
    s = x^2
    @assert hasgrade(x, 0)
    abs(basiscoeff(s))^2 / x
end
rev(x::AbstractMultivector) =
    sum(i % 4 in (0, 1) ? x[i] : -x[i] for i = 0:vectorspacedim(x))
graderev(x) = rev(gradeinv(x))
metricrev(x) = rev(metricinv(x))
metricgraderev(x) = rev(metricinv(gradeinv(x)))

#include("TreeMV.jl")
include("TreeMV2.jl")

function Base.convert(
    TMV::Type{<:TreeMV2.TreeMultivector{K,N}}, x::AM{K,N}
) where {K,N}
    T = TreeMV2.codetype(TMV)
    codes = T[]
    coeffs = K[]

    for I in eachindex(MultivectorBasis(x))
        val = basiscoeff(x, I)
        if !iszero(val)
            c = reduce(|, T(1) << (n-1) for n in Tuple(I))
            push!(codes, c)
            push!(coeffs, val)
        end
    end

    @unsafe TMV(codes, val)
end

function Base.promote_rule(
    MV1::Type{<:AM}, MV2::Type{<:AM}
)
    K1 = scalarfieldtype(MV1)
    K2 = scalarfieldtype(MV2)
    N1 = vectorspacedim(MV1)
    N2 = vectorspacedim(MV2)
    TreeMV2.TreeMultivector{promote_type(K1, K2), max(N1, N2), UInt}
end
Base.promote_rule(MV::Type{<:AM}, T::Type) = similartype(MV, T)

end # module GeomAlg
