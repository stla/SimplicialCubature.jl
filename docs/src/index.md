# SimplicialCubature.jl documentation

This package is a port of the R package **SimplicialCubature**, 
written by John P. Nolan, and which contains R translations of 
some Matlab and Fortran code written by Alan Genz.

___

A simplex is a triangle in dimension 2, a tetrahedron in dimension 3. 
This package provides two main functions: `integrateOnSimplex`, to integrate 
an arbitrary function on a simplex, and `integratePolynomialOnSimplex`, to 
get the exact value of the integral of a multivariate polynomial on a 
simplex.

A ``n``-dimensional simplex must be given by ``n+1`` vectors of length ``n``, 
which represent the simplex vertices. For example, the ``3``-dimensional 
unit simplex is encoded as follows:

```julia
S = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
```

Or you can get it by running `CanonicalSimplex(3)`.

Suppose you want to integrate the function 
```math
f(x, y ,z) = x + yz
```
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

Since the function ``f`` of this example is polynomial, you can use 
`integratePolynomialOnSimplex`:

```julia
using SimplicialCubature
using TypedPolynomials

@polyvar x y z
P = x + y*z
integratePolynomialOnSimplex(P, S)
```

Be careful if your polynomial does not involve one of the variables. 
For example if ``P(x, y, z) = x + y``, you have to encode it as a polynomial 
depending on ``z``: type `P = x + y + 0*z`.

In addition, on this example where the vertex coordinates of ``S`` and the 
coefficients of ``P`` are integer numbers, there is a more clever way to 
proceed: while `integratePolynomialOnSimplex` implements an exact proedure, 
it is not free of (small) numerical errors, but the returned value in this 
situation will be really exact if you use a polynomial with *rational* 
coefficients:

```julia
@polyvar x y z
P = 1//1*x + y*z
integratePolynomialOnSimplex(P, S)
```


## Member functions

```@autodocs
Modules = [SimplicialCubature]
Order   = [:type, :function]
```


## References

- A. Genz and R. Cools. 
*An adaptive numerical cubature algorithm for simplices.* 
ACM Trans. Math. Software 29, 297-308 (2003).

- Jean B. Lasserre.
*Simple formula for the integration of polynomials on a simplex.* 
BIT Numerical Mathematics 61, 523-533 (2021).