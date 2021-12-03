### I'm realizing that, currently,
###     the way these tree algorithms are implemented
###     I'm essentially using the type system for dispatch,
###     which is a no-no.
### I don't think there's really a better way to structure the types,
###     but if this ends up being exceedingly slow
###     it might be faster to switch the dispatch to isa conditionals.
###
### Ideas:
### - Have a `similar` type function to get a multivector of a similar type.
### - Have *intermediate representations*.
###   For example, adding together TreeMultivectors pairwise is costly,
###       since at each step we have to walk the trees
###       and construct an entirely new one.
###   We can delay this computation by just accruing a list of summands
###       and then doing the actual computation when necessary.
###   We can achieve this for incompatible operations automatically
###       by leveraging the promotion system.
### ^ Instead of the above,
###       create a macro (say `@fusetree`)
###       that looks at each arithmetic expression
###       and turns it into a fused expression
###       that only goes over the tree once.
### - Can we make the tree representation dimension agnostic?
###   Is that even worthwhile?
module TreeMV
using ..GeomAlg
using ..GeomAlg: unsafe, @unsafe

struct TreeMultivectorNode end

struct TreeMVZero end
const TMVZERO = TreeMVZero()
struct TreeMVPreleaf{K}
    left::K
    right::K
end

struct TreeMVTrunk{K}
    left::Union{TreeMVTrunk{K}, TreeMVPreleaf{K}, TreeMVZero}
    right::Union{TreeMVTrunk{K}, TreeMVPreleaf{K}, TreeMVZero}

    @unsafe TreeMVTrunk{K}(left, right) where K = new{K}(left, right)
end
function TreeMVTrunk{K}(left, right) where K
    TreeMV = TreeMVTrunk{K}
    Preleaf = TreeMVPreleaf{K}
    MaybeTrunk = Union{TreeMVTrunk, TreeMVZero}
    MaybePreleaf = Union{TreeMVPreleaf, TreeMVZero}
    @assert left isa MaybeTrunk && right isa MaybeTrunk ||
            left isa MaybePreleaf && right isa MaybePreleaf

    @unsafe TreeMVTrunk{K}(left, right)
end

const InnerTreeMV{K} = Union{TreeMVTrunk{K}, TreeMVPreleaf{K}, TreeMVZero}
const NZInnerTreeMV{K} = Union{TreeMVTrunk{K}, TreeMVPreleaf{K}}

struct TreeMultivector{K, N} <: AbstractMultivector{K, N}
    data::Union{TreeMVTrunk{K}, TreeMVPreleaf{K}, TreeMVZero}
    # Inner constructors necessary to check for correct tree depth, etc.
end
TreeMultivector{K, N}(x::TreeMultivector{K, N}) where {K, N} = x
function TreeMultivector{K, N}(scalar::K) where {K, N}
    if iszero(scalar) return TreeMultivector{K, N}(TMVZERO) end
    ret = TreeMVPreleaf{K}(zero(K), scalar)
    for _ = 2:N
        ret = TreeMVTrunk{K}(TMVZERO, ret)
    end

    TreeMultivector{K, N}(ret)
end

GeomAlg.scalarfieldtype(::Type{<:InnerTreeMV{K}}) where K = K

Base.iszero(::NZInnerTreeMV) = false
Base.iszero(::TreeMVZero) = true

Base.one(TMV::Type{<:TreeMultivector}) = TMV(oneunit(scalarfieldtype(TMV)))
Base.zero(TMV::Type{<:TreeMultivector}) = TMV(TMZERO)

treeorzero(construct, left, right) =
    if iszero(left) && iszero(right)
        TMVZERO
    else
        construct(left, right)
    end

#function GeomAlg.basiscoeff(x::TreeMultivector, basiscode::UInt)
#    K = scalarfieldtype(x)
#    # Technically shouldn't need this, but it serves as another sanity check.
#    d = vectorspacedim(x)
#
#    @assert iszero(basiscode >> d)
#
#    curnode::InnerTreeMV{K} = x.data
#    if iszero(curnode) return zero(K) end
#
#    mask = UInt(1)
#    for i = 0:(d - 1)
#        if iszero(basiscode & mask)
#            curnode = curnode.right
#        else
#            curnode = curnode.left
#        end
#        if iszero(curnode) return zero(K) end
#        mask << 1
#    end
#    @assert curnode isa K
#
#    curnode
#end

function GeomAlg.grade(x::TreeMultivector, g::Int)
    ret = g < 0 || g > vectorspacedim(x) ? TMVZERO : _grade(x.data, g)
    typeof(x)(ret)
end
_grade(::TreeMVZero, ::Int) = TMVZERO
_grade(x, g::Int) = g == 0 ? x : zero(x)
_grade(x::NZInnerTreeMV, g::Int) =
    treeorzero(unsafe(typeof(x)), _grade(x.left, g-1), _grade(x.right, g))

function GeomAlg.basiselem(TMV::Type{<:TreeMultivector}, is::Int...)
    K = scalarfieldtype(TMV)
    N = vectorspacedim(TMV)

    isvec = collect(is)
    @assert all(i -> 0 < i <= N, isvec) && allunique(isvec)
    perm = sortperm(isvec)
    
    TMV(_basiselem(K, N, 1, isevenperm(perm), isvec[perm]))
end
function _basiselem(K, N, i, shouldnotnegate, is)
    val = shouldnotnegate ? oneunit(K) : -oneunit(K)
    if i == N
        if !isempty(is) && i == is[1]
            TreeMVPreleaf{K}(val, zero(K))
        else
            TreeMVPreleaf{K}(zero(K), val)
        end
    else
        if !isempty(is) && i == is[1]
            nextnode = _basiselem(K, N, i+1, shouldnotnegate, @view is[2:end])
            TreeMVTrunk{K}(nextnode, TMVZERO)
        else
            nextnode = _basiselem(K, N, i+1, shouldnotnegate, is)
            TreeMVTrunk{K}(TMVZERO, nextnode)
        end
    end
end

### This is not strictly necessary given the above general algorithm,
###     but I came up with this first and it ought to be more efficient anyway.
function GeomAlg.basiselem(TMV::Type{<:TreeMultivector}, i::Int)
    K = scalarfieldtype(TMV)
    N = vectorspacedim(TMV)

    @assert 0 < i <= N

    ret = if i == N
        TreeMVPreleaf{K}(oneunit(K), zero(K))
    else
        TreeMVPreleaf{K}(zero(K), oneunit(K))
    end

    for _ = (i+2):N
        ret = @unsafe TreeMVTrunk{K}(TMVZERO, ret)
    end

    if i != N
        ret = @unsafe TreeMVTrunk{K}(ret, TMVZERO)
    end

    for _ = 1:(i-1)
        ret = @unsafe TreeMVTrunk{K}(TMVZERO, ret)
    end

    TMV(ret)
end

Base.:-(x::TreeMultivector) = typeof(x)(_neg(x.data))
_neg(x) = -x
_neg(::TreeMVZero) = TMVZero
_neg(x::NZInnerTreeMV) = @unsafe typeof(x)(_neg(x.left), _neg(x.right))

Base.:-(x::TMV, y::TMV) where TMV<:TreeMultivector =
    TMV(_sub(x.data, y.data))
_sub(x, y) = _sub(x, y)
_sub(::TreeMVZero, ::TreeMVZero) = TMVZERO
_sub(x::NZInnerTreeMV, ::TreeMVZero) = x
_sub(::TreeMVZero, y::NZInnerTreeMV) = _neg(y)
_sub(x::ITMV, y::ITMV) where ITMV<:NZInnerTreeMV =
    treeorzero(unsafe(ITMV), _sub(x.left, y.left), _sub(x.right, y.right))

Base.promote_rule(TMV::Type{<:TreeMultivector{K}}, ::Type{K}) where K = TMV

Base.:+(x::TMV, y::TMV, zs::TMV...) where TMV<:TreeMultivector =
    TMV(_add(scalarfieldtype(TMV), x.data, y.data, map(z -> z.data, zs)))
_add(K::Type, ::TreeMVZero, ::TreeMVZero, ::NTuple{<:Any, TreeMVZero}) = TMVZERO
function _add(K::Type, x, y, zs)
    treezero(a) = a isa TreeMVZero ? zero(K) : a
    treezero(x) + treezero(y) + (isempty(zs) ? zero(K) : sum(treezero, zs))
end
function _add(
    ::Type{K}, x::InnerTreeMV{K}, y::InnerTreeMV{K},
    zs::NTuple{<:Any, InnerTreeMV{K}}
) where K
    left(a) = a isa TreeMVZero ? TMVZERO : a.left
    right(a) = a isa TreeMVZero ? TMVZERO : a.right

    ctor = if x isa TreeMVTrunk || y isa TreeMVTrunk ||
              any(a -> a isa TreeMVTrunk, zs)
        unsafe(TreeMVTrunk{K})
    else
        TreeMVPreleaf{K}
    end

    treeorzero(ctor,
        _add(K, left(x), left(y), map(left, zs)),
        _add(K, right(x), right(y), map(right, zs))
    )
end

Base.:*(a::K, x::TreeMultivector{K}) where K = x*a
Base.:*(x::TreeMultivector{K}, a::K) where K =
    iszero(a) ? zero(x) : typeof(x)(_mul(x.data, a))
Base.:/(x::TreeMultivector{K}, a::K) where K =
    typeof(x)(_div(x.data, a))
for (op, f) in ((:*, :_mul), (:/, :_div)) @eval begin
    $f(x, a) = $op(x, a)
    $f(::TreeMVZero, _) = TMVZERO
    $f(x::NZInnerTreeMV{K}, a::K) where K =
        @unsafe typeof(x)($f(x.left, a), $f(x.right, a))
end end

### Rather than doing an actual computation here,
###     we could allow types to provide
###     their own implementation-specific dualization.
### In this case, for example,
###     this would allow us to choose the "easy" operation of swapping the tree,
###     and then we could implement the dual as just a wrapper
###     that aliases left -> right and right -> left.

GeomAlg.rightbasisdual(x::TreeMultivector) =
    typeof(x)(_rightbasisdual(x.data, 0, 0))
_rightbasisdual(::TreeMVZero, ::Int, ::Int) = TMVZERO
function _rightbasisdual(x::TreeMVPreleaf, i::Int, j::Int)
    # IDK if multiplying by -1 is better than (explicit) branching.
    # At least I know this will work for arbitrary fields.
    left = iseven(i) ? x.right : -x.right
    right = iseven(i+j) ? x.left : -x.left
    typeof(x)(left, right)
end
function _rightbasisdual(x::TreeMVTrunk, i::Int, j::Int)
    left = _rightbasisdual(x.right, i, j+1)
    right = _rightbasisdual(x.left, i+j, j)
    @unsafe typeof(x)(left, right)
end

GeomAlg.leftbasisdual(x::TreeMultivector) =
    typeof(x)(_leftbasisdual(x.data, 0, 0))
_leftbasisdual(::TreeMVZero, ::Int, ::Int) = TMVZERO
function _leftbasisdual(x::TreeMVPreleaf, i::Int, j::Int)
    # IDK if multiplying by -1 is better than (explicit) branching.
    # At least I know this will work for arbitrary fields.
    left = iseven(i+j) ? x.right : -x.right
    right = iseven(i) ? x.left : -x.left
    typeof(x)(left, right)
end
function _leftbasisdual(x::TreeMVTrunk, i::Int, j::Int)
    left = _leftbasisdual(x.right, i+j, j)
    right = _leftbasisdual(x.left, i, j+1)
    @unsafe typeof(x)(left, right)
end

### Perhaps refactor to pass each branches string down
function Base.show(io::IO, mime::MIME"text/plain", x::TreeMultivector)
    println(io, "TreeMultivector{$(repr(mime, scalarfieldtype(x)))}")
    _showplain(io, mime, x.data, 1, 0, 1, true, -1)

    nothing
end

function _showplain(io, mime, ::TreeMVZero, i, j, k, isindented, indent)
    print(io, "∅")

    nothing
end

function _showplain(io, mime, x::TreeMVPreleaf, i, j, k, isindented, indent)
    leftindent = isindented ? "  "^(indent+1) : ""
    leftval = iszero(x.left) ? "∅" : " "*repr(mime, x.left)
    println(io, leftindent, i, ":", leftval)

    rightval = iszero(x.right) ? "∅" : " "*repr(mime, x.right)
    print(io, "  "^(max(indent, 0)), j, ":", rightval)

    nothing
end

function _showplain(io, mime, x::TreeMVTrunk, i, j, k, isindented, indent)
    leftindent = isindented ? "  "^(indent+1) : ""
    print(io, leftindent, i, ":")

    _showplain(io, mime, x.left, i+1, j+k, 1, false, indent+1)
    println(io)
    _showplain(io, mime, x.right, i+1, j, k+1, true, indent)

    nothing
end

end # module TreeMV
