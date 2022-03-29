struct QuadraticForm{K, N, T}
    data::T
    QuadraticForm{K,N}(data) = new{K,N,typeof(data)}(data)
end

_qform_datatype(x::QuadraticForm) = _qform_datatype(typeof(x))
_qform_datatype(::Type{<:QuadraticForm{<:Any, <:Any, T}}) where T = T

module QForms
    abstract type QFormRepr end
    struct FunctionLike<:QFormRepr end
    struct ArrayLike<:QFormRepr end
end
using .QForms: QFormRepr

Base.@propagate_inbounds (q::QF)(i::Int) =
    _apply_qform(QFormRepr(_qform_datatype(T)), q, i)
Base.@propagate_inbounds (q::QF)(i::Int, j::Int)
    _apply_qform(QFormRepr(_qform_datatype(T)), q, i, j)

Base.@propagate_inbounds _apply_qform(::QForms.FunctionLike, q, i) = q(i)
Base.@propagate_inbounds _apply_qform(::QForms.FunctionLike, q, i, j) = q(i,j)
Base.@propagate_inbounds _apply_qform(::QForms.MatrixLike, q, i) = q[i,i]
Base.@propagate_inbounds _apply_qform(::QForms.MatrixLike, q, i, j) = q[i,j]

quadraticform(a::AM, b::AM) =
    quadraticform(promote_type(typeof(a), typeof(b)))

Cassette.@context QFormCtx
function Cassette.overdub(
    ctx::QFormCtx{<:QuadraticForm{K,N}}, ::typeof(quadraticform), x::AM{K,N}
) where {K,N}
    ctx.metadata
end

macro with_quadratic_form(form, body)
    form, body = esc.((form, body))
    :(Cassette.recurse(QFormCtx(; metadata=$form), () -> $body))
end
macro with_quadratic_form(K, N, form, body)
    K, N, form, body = esc.((K, N, form, body))
    quote
        Cassette.recurse(
            QFormCtx(; metadata=QuadraticForm{$K,$N}($form)),
            () -> $body
        )
    end
end
