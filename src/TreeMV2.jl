module TreeMV2
using TupleTools
using ..GeomAlg
using ..GeomAlg: @unsafe, isevenperm

### Should really use `SVector`s, probably.
### Description:
### - 'codes' stores `N` bits,
###   with the i^th bit corresponding to whether or not
###   the i^th basis vector is included.
### - `coeffs` stores the coefficient
###   for the corresponding basis blade in `codes`.
### Invariants:
### - `codes[i]` should correspond to `coeffs[i]`.
### - `codes` should be sorted in ascending order by `Unsigned` value;
###   this puts them in (reversed) tree order.
struct TreeMultivector{K, N, T<:Unsigned} <: AbstractMultivector{K, N}
    codes::Vector{T}
    coeffs::Vector{K}

    @unsafe TreeMultivector{K, N, T}(codes, coeffs) where {K, N, T<:Unsigned} =
        new{K, N, T}(codes, coeffs)
end
function TreeMultivector{K, N, T}(codes, coeffs) where {K, N, T<:Unsigned}
    len = length(codes)
    @assert length(coeffs) == len

    codesvec = collect(codes)
    perm = sortperm(codesvec)
    @unsafe TreeMultivector{K, N, T}(codesvec[perm], collect(coeffs)[perm])
end

codetype(x::TreeMultivector) = codetype(typeof(x))
codetype(::Type{<:TreeMultivector{<:Any, <:Any, T}}) where T<:Unsigned = T

GeomAlg.similarmvtype(TMV::Type{<:TreeMultivector}, K::Type, N::Int) =
    TreeMultivector{K, N, codetype(TMV)}
GeomAlg.similarmvtype(TMV::Type{<:TreeMultivector}, K::Type) =
    TreeMultivector{K, vectorspacedim(TMV), codetype(TMV)}
GeomAlg.similarmvtype(TMV::Type{<:TreeMultivector}, N::Int) =
    TreeMultivector{scalarfieldtype(TMV), N, codetype(TMV)}

function Base.convert(
    ::Type{AbstractMultivector{K, N}}, x::TreeMultivector
) where {K, N}
    @assert N >= vectorspacedim(x)
    @unsafe TreeMultivector{K, N}(
        convert(AbstractArray{codetype(x)}, x.codes),
        convert(AbstractArray{K}, x.coeffs)
    )
end

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


Base.:(==)(x::TreeMultivector{<:Any,N}, y::TreeMultivector{<:Any,N}) where N =
    x.codes == y.codes && x.coeffs == y.coeffs

Base.zero(TMV::Type{<:TreeMultivector}) =
    @unsafe TMV(codetype(TMV)[], scalarfieldtype(TMV)[])

function Base.one(TMV::Type{<:TreeMultivector})
    K = scalarfieldtype(TMV)
    T = codetype(TMV)
    i = one(K)
    @unsafe TreeMultivector{T, typeof(i)}([zero(T)], [i])
end

function Base.oneunit(TMV::Type{<:TreeMultivector})
    K = scalarfieldtype(TMV)
    T = codetype(TMV)
    @unsafe TMV([zero(T)], [oneunit(K)])
end

for f in (:zero, :one, :oneunit) @eval begin
    Base.$f(x::TreeMultivector) = $f(typeof(x))
end end

Base.promote_rule(TMV::Type{<:TreeMultivector{K}}, ::Type{K}) where K = TMV
Base.convert(TMV::Type{<:TreeMultivector{K}}, a::K) where K =
    @unsafe TMV([zero(K)], [a])

Base.:(==)(x::TreeMultivector, y::TreeMultivector) =
    x.codes == y.codes && x.coeffs == y.coeffs

function GeomAlg.add(x::TMV, y::TMV) where TMV<:TreeMultivector
    K = scalarfieldtype(TMV)
    T = codetype(TMV)

    xnblades, ynblades = nbasisblades(x), nbasisblades(y)

    if iszero(xnblades)
        return y
    elseif iszero(ynblades)
        return x
    end

    maxnblades = xnblades + ynblades
    codes  = Vector{T}(undef, maxnblades)
    coeffs = Vector{K}(undef, maxnblades)
    idx = xi = yi = 1
    
    @inbounds while xi <= xnblades && yi <= ynblades
        xcode, ycode = x.codes[xi], y.codes[yi]
        curcode = min(xcode, ycode)

        if xcode == curcode
            if ycode == curcode
                xysum = x.coeffs[xi] + y.coeffs[yi]
                if !iszero(xysum)
                    codes[idx] = curcode
                    coeffs[idx] = xysum
                    idx += 1; xi += 1; yi += 1
                end
            else
                codes[idx] = curcode
                coeffs[idx] = x.coeffs[xi]
                idx += 1; xi += 1
            end
        else # ycode == curcode
            codes[idx] = curcode
            coeffs[idx] = y.coeffs[yi]
            idx += 1; yi += 1
        end
    end

    if xi > xnblades
        @inbounds for yi = yi:ynblades
            codes[idx] = y.codes[yi]
            coeffs[idx] = y.coeffs[yi]
            idx += 1
        end
    else # yi > xnblades
        @inbounds for xi = xi:xnblades
            codes[idx] = x.codes[xi]
            coeffs[idx] = x.coeffs[xi]
            idx += 1
        end
    end

    resize!(codes,  idx-1)
    resize!(coeffs, idx-1)
    @unsafe TMV(codes, coeffs)
end

GeomAlg.sub(x::TreeMultivector) = @unsafe typeof(x)(x.codes, -x.coeffs)

GeomAlg.mul(_, x::TreeMultivector{K}, a::K) where K =
    @unsafe typeof(x)(x.codes, x.coeffs .* a)

function GeomAlg.mul(
    q::AbstractQuadraticForm{K,N}, x::TreeMultivector{K,N}
) where {K,N}
end

GeomAlg.div(_, x::TreeMultivector{K}, a::K) where K =
    @unsafe typeof(x)(x.codes, x.coeffs ./ a)

GeomAlg.hasgrade(x::TreeMultivector, g::Int) =
    any(==(g) ∘ count_ones, x.codes)

function GeomAlg.grade(x::TreeMultivector, g::Int)
    K = scalarfieldtype(x)
    N = vectorspacedim(x)
    T = codetype(x)

    if !(0 <= g <= N)
        zero(x)
    end

    nblades = nbasisblades(x)

    codes = Vector{T}(undef, nblades)
    coeffs = Vector{K}(undef, nblades)
    i = 1
    @inbounds for j in 1:nblades
        xcode = x.codes[j]
        if count_ones(xcode) == g
            codes[i] = xcode
            coeffs[i] = x.coeffs[j]
            i += 1
        end
    end
    @show i

    resize!(codes,  i-1)
    resize!(coeffs, i-1)

    @unsafe typeof(x)(codes, coeffs)
end

GeomAlg.basiselem(TMV::Type{<:TreeMultivector}) = oneunit(TMV)
function GeomAlg.basiselem(TMV::Type{<:TreeMultivector}, idxs::Int...)
    K = scalarfieldtype(TMV)
    N = vectorspacedim(TMV)
    T = codetype(TMV)

    if any(i -> !(1 <= i <= N), idxs) || !allunique(idxs)
        return zero(TMV)
    end

    code = reduce(|, (oneunit(T) << (n-1) for n in idxs))
    coeff = isevenperm(TupleTools.sortperm(idxs)) ? oneunit(K) : -oneunit(K)

    @unsafe TMV([code], [coeff])
end

function GeomAlg.basiscoeff(x::TreeMultivector, vecnums::NTuple{<:Any, Int})
    K = scalarfieldtype(x)
    T = codetype(x)
    code = reduce(|, (T(1) << (i-1) for i in vecnums))
    for (i, c) in enumerate(x.codes)
        if c > code
            break
        elseif c == code
            return x.coeffs[i]
        end
    end

    zero(K)
end

function _dualsign(N, T, code)
    ret = T(0)
    ip = T(0) # intervening parity
    for k = 0:(N-1)
        hasvec = (code >> k) & T(1)
        ret = ret ⊻ (ip & hasvec)
        ip = T(1) ⊻ ip ⊻ hasvec
    end

    ret
end

function GeomAlg.rightbasisdual(x::TreeMultivector)
    N = vectorspacedim(x)
    T = codetype(x)

    flipbits(c) = T(2^N-1) ⊻ c
    codes = reverse!(flipbits.(x.codes))

    coeffs = reverse!([
        iszero(_dualsign(N, T, c)) ? v : -v
        for (c, v) in zip(x.codes, x.coeffs)
    ])

    @unsafe typeof(x)(codes, coeffs)
end

function GeomAlg.leftbasisdual(x::TreeMultivector)
    N = vectorspacedim(x)
    T = codetype(x)

    flipbits(c) = T(2^N-1) ⊻ c
    codes = reverse!(flipbits.(x.codes))

    coeffs = reverse!([
        iszero(_dualsign(N, T, c)) ? -v : v
        for (c, v) in zip(x.codes, x.coeffs)
    ])

    @unsafe typeof(x)(codes, coeffs)
end

function gradeinv(x::TreeMultivector)
    coeffs = copy(x.coeffs)
    @inbounds for (i, c) in enumerate(x.codes)
        if isodd(count_ones(c))
            coeffs[i] = -coeffs[i]
        end
    end

    @unsafe typeof(x)(x.codes, coeffs)
end

### Requires a notion of sign over the base field.
### Assumes a geometric product basis.
function metricinv(q, x::TreeMultivector)
    K = scalarfieldtype(x)
    N = vectorspacedim(x)
    T = codetype(x)
    I = one(K)

    coeffs = copy(x.coeffs)
    @inbounds for (i, c) in enumerate(x.codes)
        product = I
        for j = 0:(N-1)
            sq = q(j)
            if !iszero((T(1) << j) & c) && !iszero(sq)
                product *= sq
            end
        end

        notrev = GeomAlg.notrevgrade(count_ones(c))
        coeffs[i] /= sign(notrev ? product : -product)
    end

    @unsafe typeof(x)(x.codes, coeffs)
end

function Base.show(io::IO, ::MIME"text/plain", TMV::Type{<:TreeMultivector})
    K = scalarfieldtype(TMV)
    N = vectorspacedim(TMV)
    T = codetype(TMV)
    print(io, N, "D ", K, "-TreeMultivector with ", T, " codes")
end

function Base.show(io::IO, mime::MIME"text/plain", x::TreeMultivector)
    N = vectorspacedim(x)
    K = scalarfieldtype(x)
    T = codetype(x)

    width = ndigits(N)

    function printbasisblade(code, coeff)
        if iszero(code)
            print(io, " "^((width+1)*(N-1)), lpad("", width), "| ", coeff)
        else
            paddedvectornums = Iterators.map(1:N) do j
                vectorbit = UInt(1) << (j-1)
                iszero(vectorbit & code) ? " "^width : lpad(j, width)
            end
            print(io, join(paddedvectornums, " "), "| ", coeff)
        end
    end

    show(io, mime, typeof(x))
    println()
    if iszero(x)
        print(io, "∅")
    else
        for (code, coeff) in zip((@view x.codes[1:end-1]), (@view x.coeffs[1:end-1]))
            printbasisblade(code, coeff)
            println(io)
        end
        printbasisblade(x.codes[end], x.coeffs[end])
    end

    nothing
end

end # module TreeMV2
