using SimplicialCubature
using Test
using TypedPolynomials

@testset "A polynomial." begin
    function f(x)
        return x[1]^3 + x[1]*x[2]*x[3]
    end
    @polyvar x y z
    P = x^3 + x*y*z
    S = CanonicalSimplex(3)
    I_f = integrateOnSimplex(f, S)
    I_P = integratePolynomialOnSimplex(P, S)
    @test isapprox(
        I_f.integral, I_P
    )
end

@testset "A polynomial with rational coefficients." begin
    @polyvar x y z
    P = 1//1*x + y*z
    S = CanonicalSimplex(3)
    I_P = integratePolynomialOnSimplex(P, S)
    @test I_P == 1//20
end

@testset "A function." begin
    function f(x)
        return exp(x[1] + x[2] + x[3])
    end
    S = [[0, 0, 0], [1, 1, 1], [0, 1, 1], [0, 0, 1]]
    I_f = integrateOnSimplex(f, S)
    @test isapprox(
        I_f.integral, (exp(1) - 1)^3 / 6
    )
end
