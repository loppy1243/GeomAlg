using MacroTools

function randtmv2(N, n)
    maxcode = UInt(2^N - 1)
    codes = collect(UInt(0):maxcode)
    for _ = 1:(2^N - n)
        popat!(codes, rand(eachindex(codes)))
    end
    GeomAlg.TreeMV2.TreeMultivector{Float64, N, UInt}(
        codes, rand(Float64, length(codes))
    )
end

unsafe(T::Type) = (args...; kwargs...) -> T(args...; kwargs...)

macro unsafe(ex::Expr)
    doerror() = error("Expected function definition or call expression")

    isfuncdef =
       ex.head === :function ||
       ex.head === :(=) &&
           isexpr(ex.args[1]) &&
           ex.args[1].head in (:where, :call)

    if isfuncdef
        def = splitdef(ex)
        haskey(def, :name) || doerror()

        for (k, v) in def
            def[k] = v isa Array || v isa Tuple ? map(esc, v) : esc(v)
        end

        name = def[:name]
        params = get(def, :params, nothing)
        whereparams = get(def, :whereparams, ())
        delete!(def, :name)
        delete!(def, :whereparams)
        closure = combinedef(def)
        if isnothing(params) quote
            function GeomAlg.unsafe(::Type{$name}) where {$(whereparams...)}
                $closure
            end end
        else quote
            function GeomAlg.unsafe(::Type{$name{$(params...)}}) where {$(whereparams...)}
                $closure
            end end
        end
    elseif @capture(ex, T_(xs__))
        quote
            GeomAlg.unsafe($(esc(T)))($(map(esc, xs)...))
        end
    else
        doerror()
    end
end

### Stolen and modified from `Combinatorics.levicivita`
function isevenperm(p)
    n = length(p)

    todo = trues(n)
    first = 1
    cycles = flips = 0

    while cycles + flips < n
        first = coalesce(findnext(todo, first), 0)
        (todo[first] = !todo[first]) && return 0
        j = p[first]
        (0 < j <= n) || return 0
        cycles += 1
        while j ≠ first
            (todo[j] = !todo[j]) && return 0
            j = p[j]
            (0 < j <= n) || return 0
            flips += 1
        end
    end

    iseven(flips)
end
