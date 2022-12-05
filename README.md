# SimplicialCubature

[![Build Status](https://github.com/stla/SimplicialCubature.jl/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/stla/SimplicialCubature.jl/actions/workflows/test.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/stla/SimplicialCubature.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/stla/SimplicialCubature.jl)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://stla.github.io/SimplicialCubature.jl/dev)

A simplex is a triangle in dimension 2, a tetrahedron in dimension 3. 
This package provides two main functions: `integrateOnSimplex`, to integrate 
an arbitrary function on a simplex, and `integratePolynomialOnSimplex`, to 
get the exact value of the integral of a multivariate polynomial on a 
simplex.

A `n`-dimensional simplex must be given by `n+1` vectors of length `n`, 
which represent the simplex vertices. For example, the 3-dimensional 
unit simplex is encoded as follows:

```julia
S = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
```

Or you can get it by running `CanonicalSimplex(3)`.

Suppose you want to integrate the function 
$$
f(x, y ,z) = x + yz
$$
on the unit simplex. To use `integrateOnSimplex`, you have to define $f$ 
as a function of a 3-dimensional vector:

```julia
function f(x)
  return x[1] + x[2]*x[3]
end

using SimplicialCubature
I = integrateOnSimplex(f, S)
```

Then the value of the integral is given in `I.integral`.

Since the function $f$ of this example is polynomial, you can use 
`integratePolynomialOnSimplex`:

```julia
using SimplicialCubature
using TypedPolynomials

@polyvar x y z
P = x + y*z
integratePolynomialOnSimplex(P, S)
```



