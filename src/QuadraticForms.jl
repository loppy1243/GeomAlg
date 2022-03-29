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

function Base.Matrix(q::AbstractQuadraticForm)
    N = vectorspacedim(q)
    [@inbounds q(i, j) for i=1:N, j=1:N]
end
Base.Matrix{K}(q::AbstractQuadraticForm) where K =
    convert(Matrix{K}, Matrix(q))

struct Unsafe{K, N} <: AbstractQuadraticForm{K, N}
    mat::Matrix{K}
end
Base.@propagate_inbounds (q::Unsafe)(x::Int) = q.mat[x,x]
Base.@propagate_inbounds (q::Unsafe)(x::Int, y::Int) = q.mat[x,y]
Base.Matrix(x::Unsafe) = x.mat

struct Diagonal{K, N} <: AbstractQuadraticForm{K, N}
    sig::NTuple{N, K}
end
Diagonal{K,N}(sig::Vararg{K,N}) where {K,N} = Diagonal{K,N}(sig)
Diagonal{K}(sig::NTuple{N,K}) where {K,N} = Diagonal{K,N}(sig)
Diagonal{K}(sig::K...) where K = Diagonal{K,length(sig)}(sig)
Diagonal(sig::K...) where K = Diagonal{K,length(sig)}(sig)

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
        throw(BoundsError(q, (x, y)))
    end

    sum(@inbounds q.matrix[i,x]*q.inner(i,j)*q.matrix[j,y] for i=1:N, j=1:N)
end

Base.Matrix(q::BasisChange) = transpose(q.matrix)*Matrix(q.inner)*q.matrix

VARNAMES = ['a'+i for i = 0:25]
NTERMS = 6

### This is more so a shorthand since our biliear forms are doubled
###     and our quadratic forms are not.
function Base.show(io::IO, ::MIME"text/plain", q::AbstractQuadraticForm)
    K = scalarfieldtype(q)
    N = vectorspacedim(q)
    maxvars = N <= length(VARNAMES) ? N : length(VARNAMES)

    print(io,
        nameof(typeof(q)), " ",
        N, "D ",
        nameof(K),
        if K isa DataType && !isempty(K.parameters)
            "{"*join(K.parameters, ", ")*"}"
        else
            ""
        end,
        "-form ")

    terms = []
    for j in 1:maxvars
        qj = q(j)
        if !iszero(qj)
            if isone(qj)
                push!(terms, "$(VARNAMES[j])²")
            else
                push!(terms, "$(qj)*$(VARNAMES[j])²")
            end

            if length(terms) >= NTERMS break end
        end
        for i = (j+1):maxvars
            qij = q(i, j)
            vars = sort!(VARNAMES[[i,j]])
            if !iszero(qij)
                if isone(qij)
                    push!(terms, "$(vars...)")
                else
                    push!(terms, "$(qij)*$(vars...)")
                end
            end

            if length(terms) >= NTERMS break end
        end
    end

    print(io, join(terms, " + "))
    if length(terms) >= NTERMS || N > length(VARNAMES)
        print(io, " + ...")
    end

    nothing
end

end # module QuadraticForms
