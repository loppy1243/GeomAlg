Base.promote_rule(MV::Type{<:AM{K}}, ::Type{K}) where K = MV

similarmv(x::AM, args...) = convert(similarmvtype(args...), x)
Base.convert(::Type{AM{K}}, x::AM) where K =
    convert(AM{K, vectorspacedim(x)}, x)
function Base.convert(
    ::Type{AM{<:Any, N}}, x::AM
) where N
    convert(AM{scalarfieldtype(x), N}, x)
end
