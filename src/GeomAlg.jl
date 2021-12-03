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
abstract type AbstractBlade{K, N} <: AbstractMultivector{K, N} end

struct MultivectorBasis{MV<:AbstractMultivector} end
MultivectorBasis(MV::Type{<:AbstractMultivector}) = MultivectorBasis{MV}()
MultivectorBasis(x::AbstractMultivector) = MultivectorBasis(typeof(x))
multivectortype(::Type{MultivectorBasis{MV}}) where MV<:AbstractMultivector = MV
multivectortype(x::MultivectorBasis) = multivectortype(typeof(x))
scalarfieldtype(MB::Type{<:MultivectorBasis}) =
    scalarfieldtype(multivectortype(MV))
vectorspacedim(MB::Type{<:MultivectorBasis}) =
    vectorspacedim(multivectortype(MV))

let CACHE=Dict{Tuple{DataType, Vector{Int}}, AbstractMultivector}()
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

Cassette.@context QFormCtx
function Cassette.overdub(
    ctx::QFormCtx{QF}, ::typeof(quadraticform), x::AbstractMultivector{K,N}
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
    ctx::QFormCtx, ::typeof(quadraticform), x
)
    if hasmethod(quadraticform, Tuple{typeof(x)})
        quadraticform(x)
    else
        error("No quadraticform associated with $(typeof(x))")
    end
end

macro with_quadratic_form(form, body)
    :(Cassette.recurse(QFormCtx(; metadata=$(esc(form))), () -> $(esc(body))))
end

include("interface.jl")

∧(args...) = outerprod(args...)
a  ⊣  b = leftcontract(a, b)
a  ⊢  b = rightcontract(a, b)
a  ×₂ b =  a*b - b*a
a  ×  b = (a*b - b*a)/2
a  ∥₂ b =  a*b + b*a
a  ∥  b = (a*b + b*a)/2

scalarfieldtype(T::Type) =
    error("scalarfieldtype is not implemented for type $T")
scalarfieldtype(x) = scalarfieldtype(typeof(x))
scalarfieldtype(::Type{<:AbstractQuadraticForm{K}}) where K = K
scalarfieldtype(MV::Type{<:AbstractMultivector{K}}) where K = K

vectorspacedim(T::Type) = error("vectorspacedim is not implemented for type $T")
vectorspacedim(x) = vectorspacedim(typeof(x))
vectorspacedim(::Type{<:AbstractQuadraticForm{<:Any, N}}) where N = N
vectorspacedim(MV::Type{<:AbstractMultivector{<:Any, N}}) where N = N

Base.getindex(x::AbstractMultivector, g::Int) = grade(x, g)

similarmv(x::AbstractMultivector, args...) = convert(similartype(args...), x)
Base.convert(::Type{AbstractMultivector{K}}, x::AbstractMultivector) where K =
    convert(AbstractMultivector{K, vectorspacedim(x)}, x)
function Base.convert(
    ::Type{AbstractMultivector{<:Any, N}}, x::AbstractMultivector
) where N
    convert(AbstractMultivector{scalarfieldtype(x), N}, x)
end

struct BasisIndexIter{N}
    GeomAlg.eachbasisindex(N) = new{N}()
end
struct MVIndexView{N,A<:AbstractArray{Int,1}}
    parent::A
    last::Int
end
MVIndexView{N}(a::AbstractArray{Int,1}, last::Int) where N =
    MVIndexView{N,typeof(a)}(a, last)
Base.length(::BasisIndexIter{N}) where N = 2^N
Base.eltype(::BasisIndexIter{N}) where N =
    if N < 15
        MVIndexView{N,MVector{N,Int}}
    else
        MVIndexView{N,Vector{Int}}
    end

include("basisindexiter.jl")
#include("basisindexiter_reference.jl")

function Base.:(==)(
    x::AbstractMultivector{<:Any,N}, y::AbstractMultivector{<:Any,N}
) where N
    all(basiscoeff(x, I) == basiscoeff(x, I) for I in eachbasisindex(N))
end

Base.:-(x::MV, y::MV) where MV<:AbstractMultivector = x + -y
Base.:*(a::K, x::AbstractMultivector{K}) where K = x*a
Base.:/(x::Union{MV, K}, y::MV) where {K, MV<:AbstractMultivector{K}} = x*inv(y)
for op in (:+, :-, :*, :/) @eval begin
    Base.$op(x::AbstractMultivector, y::AbstractMultivector) =
        $op(promote(x,y)...)
end end
for op in (:+, :-) @eval begin
    Base.$op(x::AbstractMultivector, y) = $op(promote(x,y)...)
    Base.$op(x, y::AbstractMultivector) = $op(promote(x,y)...)
end end
for op in (:*, :/) @eval begin
    Base.$op(x::AbstractMultivector, y) = $op(x, convert(scalarfieldtype(x), y))
    Base.$op(x, y::AbstractMultivector) = $op(convert(scalarfieldtype(y), x), y)
end end

basiselem(MV::Type{<:AbstractMultivector}) = one(MV)

basiscoeff(x::AbstractMultivector, I::CartesianIndex) = basiscoeff(x, Tuple(I))
basiscoeff(x::AbstractMultivector, idxs::Int...) = basiscoeff(x, idxs)

hasgrade(x::AbstractMultivector) = !iszero(x[0])

gradeinv(x::AbstractMultivector) =
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
    TMV::Type{<:TreeMV2.TreeMultivector{K,N}}, x::AbstractMultivector{K,N}
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
    MV1::Type{<:AbstractMultivector}, MV2::Type{<:AbstractMultivector}
)
    K1 = scalarfieldtype(MV1)
    K2 = scalarfieldtype(MV2)
    N1 = vectorspacedim(MV1)
    N2 = vectorspacedim(MV2)
    TreeMV2.TreeMultivector{promote_type(K1, K2), max(N1, N2), UInt}
end
Base.promote_rule(MV::Type{<:AbstractMultivector}, T::Type) = similartype(MV, T)

end # module GeomAlg
