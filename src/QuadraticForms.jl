### To allow arbitrary forms,
###     we have to carry it around with each multivector.
### For the special case of an orthonormalizable field,
###     (more generally, an isbitstype field)
###     we can carry the form through the type system as an NTuple{Dim, Int}.
### Is is possible to generically program for both of these cases
###     in an ergonomic way?
module QuadraticForms
    import LinearAlgebra
    using ..GeomAlg
    using ..GeomAlg: @unsafe
    export AbstractQuadraticForm

    abstract type AbstractQuadraticForm{K, N} end
    Base.@propagate_inbounds (q::AbstractQuadraticForm)(x::Int) = q(x, x)

    Cassette.@context QFormCtx
    function Base.@propagate_inbounds Cassette.overdub(
        ctx::QFormCtx{QF}, ::typeof(quadraticform)
    ) where {K,N,QF<:AbstractQuadraticForm{K,N}}
        ctx.metadata
    end

    struct Nonexistent <: AbstractQuadraticForm{Nothing, 0} end

    (::Nonexistent)(x::Int) = error("Default quadratic form not set")
    (::Nonexistent)(x::Int, y::Int) = error("Default quadratic form not set")

    function show(io::IO, ::MIME"text/plain", ::Nonexistent)
        print(io, "Nonexistent 0D Nothing-form")
    end

    function Base.Matrix(q::AbstractQuadraticForm)
        N = vectorspacedim(q)
        [@inbounds q(i, j) for i=1:N, j=1:N]
    end
    Base.Matrix{K}(q::AbstractQuadraticForm{K}) where K = Matrix(q)

    const DEFAULT = Ref{AbstractQuadraticForm}(Nonexistent())

    struct Diagonal{K, N} <: AbstractQuadraticForm{K, N}
        sig::NTuple{N, K}
    end

    Base.@propagate_inbounds (q::Diagonal)(x::Int) = q.sig[x]
    Base.@propagate_inbounds (q::Diagonal)(x::Int, y::Int) =
        x == y ? q.sig[x] : zero(scalarfieldtype(q))

    Base.Matrix(q::Diagonal) = LinearAlgebra.Diagonal(collect(q.sig))

    struct BasisChange{
        K, N, QT<:AbstractQuadraticForm{K, N}
    } <: AbstractQuadraticForm{K, N}
        inner::QT
        matrix::Matrix{K} # (new <- standard) change-of-basis matrix

        @unsafe function BasisChange{K, N, QT}(inner, matrix) where {
            K, N, QT<:AbstractQuadraticForm
        }
            new{K, N, QT}(inner, matrix)
        end
    end
    function BasisChange{K, N, QT}(inner::QT, matrix::Matrix{K}) where
             {K, N, QT<:AbstractQuadraticForm{K, N}}
        @assert size(matrix) == (N, N)
        # Technically should be checking for invertibility too...
        @unsafe BasisChange{K, N, QT}(inner, matrix)
    end

    Base.@propagate_inbounds function (q::BasisChange)(x::Int, y::Int)
        N = vectorspacedim(q)

        @boundscheck if !(1 <= x <= N && 1 <= y <= N)
            throw(BoundError(q, (x, y)))
        end

        sum(@inbounds q.matrix[i,x]*q.inner(i,j)*q.matrix[j,y] for i=1:N, j=1:N)
    end

    Base.Matrix(q::BasisChange) = transpose(q.matrix)*Matrix(q.inner)*q.matrix

    function Base.show(io::IO, ::MIME"text/plain", q::AbstractQuadraticForm)
        K = scalarfieldtype(q)
        N = vectorspacedim(q)
        print(io, nameof(typeof(q)), " ", N, "D ", nameof(K), "-form ")
        terms = []
        for j in 1:N
            qj = q(j)
            if !iszero(qj) push!(terms, "$(qj)X$(j)²") end
            for i = (j+1):N
                qij = q(i, j)
                if !iszero(qij)
                    push!(push!(terms, "$(qij)X$(i)*X$(j)"))
                end
            end
        end
        print(io, join(terms, " + "))

        nothing
    end
end # module QuadraticForms
