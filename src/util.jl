module Utils

using MacroTools
export @unsafe, isevenperm

struct Unsafe{T} end

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
            def[k] = v isa Tuple || v isa Array ? map(esc, v) : esc(v)
        end
        name = def[:name]
        def[:name] = :Unsafe
        def[:params] = 
            if haskey(def, :params)
                (:($name{$(def[:params]...)}),)
            else
                (name,)
            end

        combinedef(def)
    elseif @capture(ex, T_(xs__))
        T = esc(T)
        xs = map(esc, xs)
        quote
            Unsafe{$T}($(xs...))
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

end # module Utils
