Base.iterate(x::BasisIndexIter{N}) where N =
    if N == 0
        (CartesianIndex(), (nidxs=1, idxs=MVector{0,Int}()))
    else
        (CartesianIndex(), (nidxs=1, idxs=MVector(0:(N-1)...)))
    end
function Base.iterate(x::BasisIndexIter{N}, st) where N
    if N == 1 && st.nidxs == 1
        return (CartesianIndex(1), (nidxs=2, idxs=st.idxs))
    elseif st.nidxs >= N
        return nothing
    end

    if st.nidxs >= 2
        if st.idxs[st.nidxs] >= N
            st.idxs[st.nidxs-1] += 1
            st.idxs[st.nidxs] = st.idxs[st.nidxs-1] + 1
        else
            st.idxs[st.nidxs] += 1
            return (CartesianIndex(ntuple(i -> st.idxs[i], st.nidxs)), st)
        end
    end

    for i = (st.nidxs-1):-1:2
        if st.idxs[i] > N - st.nidxs + i
            st.idxs[i-1] += 1
            for j = i:st.nidxs
                st.idxs[j] = st.idxs[i-1] + (j-i+1)
            end
            st.idxs[i] = st.idxs[i-1] + 1
        else
            return (CartesianIndex(ntuple(i -> st.idxs[i], st.nidxs)), st)
        end
    end

    if st.nidxs == 1 ? st.idxs[1] >= N-st.nidxs+1 : st.idxs[1] > N-st.nidxs+1
        st.idxs[1:(st.nidxs+1)] .= 1:(st.nidxs+1)
        (CartesianIndex(ntuple(identity, st.nidxs+1)),
         (nidxs=st.nidxs+1, idxs=st.idxs))
    else
        if st.nidxs == 1; st.idxs[1] += 1 end
        (CartesianIndex(ntuple(i -> st.idxs[i], st.nidxs)), st)
    end
end
