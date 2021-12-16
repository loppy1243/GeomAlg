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

Base.iterate(x::BasisIndexIter{N}) where N =
    if N == 0
        idxs = MVector{0,Int}()
        (MVIndexView{N}(idxs, 0), (val=Val(1), idxs=idxs))
    elseif N < 15
        idxs = MVector{N,Int}(undef)
        idxs[1] = 0
        (MVIndexView{N}(idxs, 0), (val=Val(1), idxs=idxs))
    else
        idxs = Vector{Int}(undef, N)
        idxs[1] = 0
        (MVIndexView{N}(idxs, 0), (val=Val(1), idxs=idxs))
    end
@generated function Base.iterate(
    x::BasisIndexIter{N},
    st::NamedTuple{(:val, :idxs), Tuple{Val{G}, V}}
) where {N,G,V<:Union{MVector{N,Int},Vector{Int}}}
    if N == 1 && G == 1
        return :(MVIndexView{$N,$V}(st.idxs, 1), (val=$(Val(2)), idxs=st.idxs))
    elseif G >= N 
        return nothing
    end

    firstif = if G >= 2 quote
        if st.idxs[$G] >= N
            st.idxs[$(G-1)] += 1
            st.idxs[$G] = st.idxs[$(G-1)] + 1
        else
            st.idxs[$G] += 1
            return (MVIndexView{$N,$V}(st.idxs, $G), st)
        end end
    else
        :()
    end

    ifs = Base.Generator((G-1):-1:2) do i quote
        if st.idxs[$i] > $(N-G+i)
            st.idxs[$i-1] += 1
            @nexprs $(G-i+1) j -> st.idxs[$i+j-1] = st.idxs[$i-1] + j
            st.idxs[$i] = st.idxs[$i-1] + 1
        else
            return (MVIndexView{$N,$V}(st.idxs, $G), st)
        end
    end end
        
    quote
        @inbounds begin
        $firstif
        $(ifs...)
        if $(G == 1 ? :(>=) : :>)(st.idxs[1], $(N-G+1))
            @nexprs $(G+1) i -> st.idxs[i] = i
            (MVIndexView{$N,$V}(st.idxs, $(G+1)), (val=$(Val(G+1)), idxs=st.idxs))
        else
            $(G == 1 ? :(st.idxs[1] += 1) : :())
            (MVIndexView{$N,$V}(st.idxs, $G), st)
        end
        end # inbounds
    end
end
