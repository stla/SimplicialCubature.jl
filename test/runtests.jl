using SimplicialCubature
using Test
using TypedPolynomials

@testset "SimplicialCubature.jl" begin
    function f(x)
        return x[1] + x[2]*x[3]
    end
    @polyvar x y z
    P = x + y*z
    S = CanonicalSimplex(3)
    I_f = integrateOnSimplex(f, S)
    I_P = integratePolynomialOnSimplex(P, S)
    @test isapprox(
        I_f.integral, I_P
    )
end
