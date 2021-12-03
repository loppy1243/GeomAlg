using BenchmarkTools, GLMakie
using GeomAlg.TreeMV2: _dualsign, _dualsign2

function bench(Nmax, T)
    Ns = 1:Nmax
    ds1 = Float64[]
    ds2 = Float64[]
    code = rand(T)
    function _bench(::Val{N}, ::Val{S}) where {N, S}
        f1(c) = _dualsign(N, S, c)
        f2(c) = _dualsign2(N, S, c)

        f1(code); f2(code)

        print("Benchmarking _dualsign($N, ...) ...")
        bm = @benchmark $f1($code)
        push!(ds1, minimum(bm.times))

        print("Done\nBenchmarking _dualsign2($N, ...)...")
        bm = @benchmark $f2($code)
        push!(ds2, minimum(bm.times))

        println("Done\n")
    end
    for N in Ns
        _bench(Val(N), Val(T))
    end

    colors = Makie.wong_colors()
    groups = repeat([1,2]; inner=length(Ns))
    figaxplt = barplot(
        repeat(Ns; outer=2),
        [ds1; ds2];
        dodge=groups,
        color=colors[groups],
        axis=(xticks=(Ns, string.(Ns)),)
    )
    elements = [PolyElement(;polycolor=colors[i]) for i in (1,2)]
    labels = ["_dualsign", "_dualsign2"]
    Legend(figaxplt.figure[1,2], elements, labels)

    figaxplt
end
