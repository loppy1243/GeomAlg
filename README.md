# Performant Geometric Algebra over an Arbitrary Field

The goals of this package are to:
* implement primitives and useful procedures
    for arbitrary geometric/Clifford algebras
    (including those of degenerate signatures).
* allow an arbitrary base field __including those of characteristic 2__.
* do all the above performantly.

The basic multivector type, `TreeMultivector`,
  is inspired by [[1](#Breuils2017)].
I also would like to try implementing [[2](#DeKeninck2020)].

## How to Deal with Quadratic Forms
This has been an interesting problem so far,
  and I think I've hit upon a solution that makes me happy.

The basic issue is, essentially,
  "how do we _ergonomically_ pass around a quadratic form?"

Your first reaction is probably to attach it to the multivector type.
If we were just working with real numbers (as floats),
  this would be great since floats are `isbitstype`
  and so we can pass the quadratic form around through the type system;
  but if we want arbitrary fields,
  then we can't assume our field is `isbitstype`.
So attaching to the multivector type would have to mean
  carrying the quadratic form around as part of the multivector data...
  so now _every_ multivector has a copy of the quadratic form attached to it.
Perhaps this is ok at the end of the day,
  but it's sort of yucky given how we can (and usually do)
  think of a Clifford algebra as the exterior algebra
  together with the geometric product;
  the multivector data should be separate from product prescription.

We could just explicitly pass a quadratic form around everywhere,
  but this is anti-ergonomic since then, for example,
  it wouldn't even be possible to use multiplication syntax.

The solution I've settled on
  (barring future discovery of untenability)
  is the following:
There is a function `quadraticform(x)` where `x` is a multivector.
This function will return a quadratic form for use with `x`;
  the magic, however, comes from the
  [Cassette.jl](https://github.com/JuliaLabs/Cassette.jl)
  package.
With Cassette,
  we can define a macro `@quadraticform q expr`
  where `q` evaluates to a quadratic form
  and wherein _every compatible occurrence_
  of `quadraticform(x)` called within `expr`
  is modified to return the value of `q`.
(Doesn't even have to be a top-level call!)

An example would be nice:
```julia
# In reality, at this level of (lack of) abstraction,
#   we would check that `a` and `b` are the same type.
⊣(a, b) = leftcontract(quadraticform(a), a, b)
function _leftcontract(q, a, b)
    # Compute the left contraction defined by `q`
end

inv(a) = inv(quadraticform(a), a)
function inv(q, a)
    # Compute the multivector inverse defined by `q`
end

# Project the multivector `a` onto the multivector `b`.
project(a, b) = (a ⊣ b) ⊣ inv(b)

# `x` and `y` multivectors
function mygeomalgroutine(x, y)
    I = one(scalarfieldtype(x))
    q = QuadraticForms.Diagonal(I, I, I)
    
    @with_quadratic_form q begin
        # do some stuff...
        z = project(x, y)
        # do some more stuff...
    end
end
```
The left contractions and multivector inverse are metric operations,
  and require a quadratic form to be specified.
In `mygeomalgroutine`,
  __that quadratic form will be the `q`__ which is defined there.
Notice how we don't even need to worry about the quadratic form
  when we're writing generic algorithms like `project`,
  and if we want to dispatch on it
  then we just pass a call to `quadraticform` through a function barrier
  like in `leftcontract` and `inv`.
Most importantly,
  this setup will (_cough_ should _cough_) incur __no overhead__
  except for possibly slightly longer compile times;
  in Julia, this could very well mean that computation of values of `q`
  __will be completely inlined__.

We get other niceties as well.
```julia
quadraticform(_) = QuadraticForms.Diagonal(1.0, 1.0)
```
This will define a quadratic form to be used everywhere,
  so that `@with_quadratic_form` isn't necessary.
(A caution though: this does have precedence over `@with_quadratic_form`.)
For convenience, we can also do
```julia
quadraticform(K, N) = QuadraticForms.Diagonal(ones(K, N)...)
```
And in the case where you _do_
  want to attach a quadratic form to a given multivector type,
  this system also allows
```julia
quadraticform(x::MultivectorType1) = # ...
quadraticform(x::MultivectorType2) = # ...
quadraticform(x::MultivectorType3) = # ...
```

# References
<a id="Breuils2017">
[1] S. Breuils, V. Nozick, L. Fuchs,
  A Geometric Algebra Implementation using Binary Tree,
  Adv. Appl. Clifford Algebras. 27 (2017) 2133–2151.
</a>
<br/><br/>
<a id="DeKeninck2020">
[2] S. De Keninck, L. Dorst, Hyperwedge,
  in: N. Magnenat-Thalmann, C. Stephanidis, E. Wu, D. Thalmann, B. Sheng,
  J. Kim, G. Papagiannakis, M. Gavrilova (Eds.),
  Advances in Computer Graphics,
  Springer International Publishing, Cham, 2020: pp. 549–554.
</a>
