function Base.:(==)(
    x::AM{<:Any,N}, y::AM{<:Any,N}
) where N
    all(basiscoeff(x, I) == basiscoeff(x, I) for I in eachbasisindex(N))
end

∧(args...) = outerprod(args...)
Base.:+(a::AM)  = a
Base.:-(a::AM)  = sub(a)
for (op, f) in ((:+, :add), (:-, :sub), (:/, :div), (:*, :mul)) @eval begin
    Base.$op(a::AM, b::AM) = $f(a, b)
    Base.$op(a::AM, b)     = $f(a, b)
    Base.$op(a,     b::AM) = $f(a, b)
end end
a  ⊣  b = leftcontract(a, b)
a  ⊢  b = rightcontract(a, b)
a  ×₂ b = antisymmprod2(a, b)
a  ×  b = antisymmprod2(a, b)/2
a  ∥₂ b = symmprod2(a, b)
a  ∥  b = symmprod2(a, b)/2
for f in (
    :mul, :div, :leftcontract, :rightcontract, :symmprod2, :antisymmprod2
) @eval begin
    $f(a::AM, b::AM) = $f(quadraticform(a, b), a, b)
    $f(q, a, b) = $f(q, promote(a, b)...)
end end

mul(q, a::T, b::T) where T = error("multiplication unimplemented")
mul(a::AM{K}, b::K) where K = error("right scalar multiplication unimplemented")
function mul(a::AM, b)
    K = promote_type(scalarfieldtype(a), typeof(b))
    a′ = convert(AM{K}, a)
    mul(a′, b)
end
mul(a::K, b::AM{K}) where K = error("left scalar multiplication unimplemented")
function mul(a, b::AM)
    K = promote_type(scalarfieldtype(b), typeof(a))
    b′ = convert(AM{K}, b)
    mul(a, b′)
end

div(a::AM{K}, b::K) where K = mul(a, Base.inv(b))
function div(a::AM, b)
    K = promote_type(scalarfieldtype(a), typeof(b))
    a′ = convert(AM{K}, a)
    div(a′, b)
end
div(a, b::AM) = mul(a, inv(quadraticform(b), b))
div(q, a::T, b::T) where T = mul(q, a, inv(q, b))

inv(a) = Base.inv(a)
inv(a::AM) = inv(quadraticform(a), a)

leftcontract(a::AM, b) = zero(a)
leftcontract(a, b::AM) = mul(a, b)
rightcontract(a::AM, b) = mul(a, b)
rightcontract(a, b::AM) = zero(b)
symmprod2(a::AM, b) = mul(a, b)
symmprod2(a, b::AM) = mul(a, b)
antisymmprod2(a::AM, b) = zero(a)
antisymmprod2(a, b::AM) = zero(b)

sub(x, y) = add(x, -y)
