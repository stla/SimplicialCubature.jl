var documenterSearchIndex = {"docs":
[{"location":"#SimplicialCubature.jl-documentation","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.jl documentation","text":"","category":"section"},{"location":"","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.jl documentation","text":"This package is a port of the R package SimplicalCubature,  written by John P. Nolan, and which contains R translations of  some Matlab and Fortran code written by Alan Genz.","category":"page"},{"location":"","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.jl documentation","text":"___","category":"page"},{"location":"","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.jl documentation","text":"A simplex is a triangle in dimension 2, a tetrahedron in dimension 3.  This package provides two main functions: integrateOnSimplex, to integrate  an arbitrary function on a simplex, and integratePolynomialOnSimplex, to  get the exact value of the integral of a multivariate polynomial on a  simplex.","category":"page"},{"location":"","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.jl documentation","text":"A n-dimensional simplex must be given by n+1 vectors of length n,  which represent the simplex vertices. For example, the 3-dimensional  unit simplex is encoded as follows:","category":"page"},{"location":"","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.jl documentation","text":"S = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]","category":"page"},{"location":"","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.jl documentation","text":"Or you can get it by running CanonicalSimplex(3).","category":"page"},{"location":"","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.jl documentation","text":"Suppose you want to integrate the function ","category":"page"},{"location":"","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.jl documentation","text":"f(x y z) = x + yz","category":"page"},{"location":"","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.jl documentation","text":"on the unit simplex. To use integrateOnSimplex, you have to define f  as a function of a 3-dimensional vector:","category":"page"},{"location":"","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.jl documentation","text":"function f(x)\n  return x[1] + x[2]*x[3]\nend\n\nusing SimplicialCubature\nI = integrateOnSimplex(f, S)","category":"page"},{"location":"","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.jl documentation","text":"Then the value of the integral is given in I.integral.","category":"page"},{"location":"","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.jl documentation","text":"Since the function f of this example is polynomial, you can use  integratePolynomialOnSimplex:","category":"page"},{"location":"","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.jl documentation","text":"using SimplicialCubature\nusing TypedPolynomials\n\n@polyvar x y z\nP = x + y*z\nintegratePolynomialOnSimplex(P, S)","category":"page"},{"location":"","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.jl documentation","text":"Be careful if your polynomial does not involve one of the variables.  For example if P(x y z) = x + y, you have to encode it as a polynomial  depending on z: type P = x + y + 0*z.","category":"page"},{"location":"#Member-functions","page":"SimplicialCubature.jl documentation","title":"Member functions","text":"","category":"section"},{"location":"","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.jl documentation","text":"Modules = [SimplicialCubature]\nOrder   = [:type, :function]","category":"page"},{"location":"#SimplicialCubature.CanonicalSimplex-Tuple{Any}","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.CanonicalSimplex","text":"CanonicalSimplex(n)\n\nCanonical n-dimensional simplex.\n\nArgument\n\nn: positive integer\n\n\n\n\n\n","category":"method"},{"location":"#SimplicialCubature.integrateOnSimplex-Union{Tuple{T}, Tuple{Function, Union{Array{Vector{T}, 1}, Array{Array{Vector{T}, 1}, 1}}}} where T<:Real","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.integrateOnSimplex","text":"integrateOnSimplex(f, S; dim, maxEvals, absError, tol, rule, info, fkwargs...)\n\nIntegration of a function over one or more simplices.\n\nArguments\n\nf: function to be integrated; must return a real scalar value or a real vector\nS: simplex or vector of simplices; a simplex is given by n+1 vectors of dimension n\ndim: number of components of f\nmaxEvals: maximum number of calls to f\nabsError: requested absolute error\ntol: requested relative error\nrule: integration rule, an integer between 1 and 4; a 2*rule+1 degree integration rule is used\ninfo: Boolean, whether to print more info\nfkwargs: keywords arguments of f\n\n\n\n\n\n","category":"method"},{"location":"#SimplicialCubature.integratePolynomialOnSimplex-Tuple{Any, Any}","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.integratePolynomialOnSimplex","text":"integratePolynomialOnSimplex(P, S)\n\nExact integral of a polynomial over a simplex.\n\nArgument\n\nP: polynomial\nS: simplex, given by a vector of n+1 vectors of dimension n, the simplex vertices \n\n\n\n\n\n","category":"method"},{"location":"#References","page":"SimplicialCubature.jl documentation","title":"References","text":"","category":"section"},{"location":"","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.jl documentation","text":"A. Genz and R. Cools. ","category":"page"},{"location":"","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.jl documentation","text":"An adaptive numerical cubature algorithm for simplices.  ACM Trans. Math. Software 29, 297-308 (2003).","category":"page"},{"location":"","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.jl documentation","text":"Jean B. Lasserre.","category":"page"},{"location":"","page":"SimplicialCubature.jl documentation","title":"SimplicialCubature.jl documentation","text":"Simple formula for the integration of polynomials on a simplex.  BIT Numerical Mathematics 61, 523-533 (2021).","category":"page"}]
}