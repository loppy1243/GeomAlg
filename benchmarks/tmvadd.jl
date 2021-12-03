using GeomAlg, BenchmarkTools, GLMakie

function bench(N, s)
    ns = 1:s:2^N
    ntrials = length(ns)
    newvec() = Vector{Float64}(undef, ntrials)
#    tmv1 = (mins=newvec(), meds=newvec())
    tmv2 = (mins=newvec(), meds=newvec())

    rand(2^N) + rand(2^N)
    print("Benchmarking control...")
    bc = @benchmark x+y+z setup=((x,y,z)=(rand($(2^N)),rand($(2^N)),rand($(2^N))))
    println("Done\n")
    control = (min=minimum(bc.times), med=median(bc.times))

    for (i, n) in enumerate(ns)
#        randtmv(N, n) + randtmv(N, n)
#        print("Benchmarking tmv1 with ", n, " elements...")
#        bc = @benchmark x+y setup=((x,y)=(randtmv($N,$n),randtmv($N,$n))) samples=100
#        println("Done")
#        tmv1.mins[i] = minimum(bc.times)
#        tmv1.meds[i] = median(bc.times)
    
        randtmv2(N, n) + randtmv2(N, n)
        print("Benchmarking tmv2 with ", n, " elements...")
        bc = @benchmark x+y+z setup=((x,y,z)=(randtmv2($N,$n),randtmv2($N,$n),randtmv2($N,$n))) samples=100
        println("Done")
        tmv2.mins[i] = minimum(bc.times)
        tmv2.meds[i] = median(bc.times)
        println()
    end

    print("Plotting...")
    colors = Makie.wong_colors()
    groups = repeat([1,2]; inner=ntrials)
    figaxplt = barplot(
        repeat(ns; outer=2),
        [tmv2.mins; tmv2.meds];
        dodge=groups,
        color=groups,
        axis=(xticks=(ns, string.(ns)),)
    )
    fig, ax, _ = figaxplt

    cminline = hlines!(ax, [control.min]; color=colors[3])
    cmedline = hlines!(ax, [control.med]; color=colors[4])

    labels = ["tmv2 min", "tmv2 med", "control min", "control med"]
    elements = [
        [PolyElement(;polycolor=colors[i]) for i in (1,2)];
        [cminline, cmedline]
    ]

    Legend(fig[1,2], elements, labels)
    println("Done")

    figaxplt
end
