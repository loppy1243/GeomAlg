import Base.Threads
using StaticArrays

struct BasisIndex{N}
    last::Int
    idx::SVector{N,Int}

    @unsafe function BasisIndex{N}(last, idx::SVector{N,Int}) where N
        new{N}(last, idx)
    end
end
@unsafe function BasisIndex{N}(last, idx::NTuple{<:Any,Int}) where N
    idxmut = MVector{N,Int}(undef)

    for i = 1:last
        @inbounds idxmut[i] = idx[i]
    end

    @unsafe BasisIndex{N}(last, SVector{N,Int}(idxmut))
end
function BasisIndex{N}(last, idx::Vararg{Int}) where N
    @assert last <= N && last <= length(idx)
    @unsafe BasisIndex{N}(last, idx)
end
BasisIndex{N}(idx::Int...) where N =
    @unsafe BasisIndex{N}(min(N, length(idx)), idx)
BasisIndex(idx::Int...) = @unsafe BasisIndex{length(idx)}(length(idx), idx)

Base.:(==)(x::BasisIndex, y::BasisIndex) =
    x.last == y.last && all(@inbounds x.idx[i] == y.idx[i] for i = 1:x.last)

function Base.show(io::IO, ::MIME"text/plain", I::BasisIndex{N}) where N
    print(io, "BasisIndex{$N}(")
    if !iszero(I.last)
        print(io, join(I.idx[1:I.last], ", "))
    end
    print(io, ')')
end

struct BasisIndexIter{N}
    GeomAlg.eachbasisindex(N) = new{N}()
end
GeomAlg.eachbasisindex(T::Type) = eachbasisindex(vectorspacedim(T))
GeomAlg.eachbasisindex(x::AbstractMultivector) =
    eachbasisindex(vectorspacedim(x))

Base.eltype(::BasisIndexIter{N}) where N = BasisIndex{N}
Base.length(::BasisIndexIter{N}) where N = 2^N

let STATES =
    SVector{Threads.nthreads(), Dict{Int, MVector{<:Any, Int}}}(
        Dict{Int, MVector{<:Any, Int}}()
        for _ = 1:Threads.nthreads()
    )
function Base.iterate(::BasisIndexIter{N}) where N
    state::MVector{N,Int} = get!(STATES[Threads.threadid()], N) do
        MVector{N,Int}(undef)
    end
    I = @unsafe BasisIndex{N}(0, SVector(state))

    if N != 0
        @inbounds state[1] = 1
    end
    (I, (1, state))
end
end

function Base.iterate(iter::BasisIndexIter{N}, (grade, state)) where N
    if N == 1 && grade == 1
        I = @unsafe BasisIndex{N}(grade, SVector(state))

        (I, (2, state))
    elseif grade >= N
        nothing
    elseif @inbounds(state[grade]) > N
        if @inbounds(state[1]) > N - grade
            _rollover_next_grade(iter, grade, state)
        else
            _rollover(iter, grade, state)
        end
    else
        I = @unsafe BasisIndex{N}(grade, SVector(state))
        @inbounds state[grade] += 1

        (I, (grade, state))
    end
end

@inline function _rollover_next_grade(::BasisIndexIter{N}, grade, state) where N
    for j = 1:(grade+1)
        @inbounds state[j] = j
    end

    I = @unsafe BasisIndex{N}(grade+1, SVector(state))
    @inbounds state[grade+1] += 1

    (I, (grade+1, state))
end

@inline function _rollover(::BasisIndexIter{N}, grade, state) where N
    i = 0
    for j in 1:(grade-1) if @inbounds(state[j+1]) > N - grade + j
        i = j
        break
    end end

    @inbounds state[i] += 1
    for j = (i+1):grade
        @inbounds state[j] = state[j-1] + 1
    end

    I = @unsafe BasisIndex{N}(grade, SVector(state))
    @inbounds state[grade] += 1

    (I, (grade, state))
end
