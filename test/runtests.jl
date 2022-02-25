using Test
using GeomAlg

@testset "Quadratic Forms" begin
    VARS = QuadraticForms.VARNAMES
    NTERMS = QuadraticForms.NTERMS

    @test repr("text/plain", QuadraticForms.Diagonal(1.0, 2.0, 3.0)) == """
        Diagonal 3D Float64-form $(VARS[1])² + 2.0*$(VARS[2])² + 3.0*$(VARS[3])²"""
    @test repr("text/plain", QuadraticForms.Unsafe{Int64,5}(
        Int64[0 2 0 0 0; 2 0 0 0 0; 0 0 0 1 0; 0 0 1 0 0; 0 0 0 0 5]
    )) == """
        Unsafe 5D Int64-form 2*$(VARS[1])$(VARS[2]) + $(VARS[3])$(VARS[4]) + 5*$(VARS[5])²"""
    @test repr("text/plain",
        QuadraticForms.Diagonal(ntuple(_ -> 1//1, NTERMS + 1))
    ) == """
        Diagonal $(NTERMS+1)D Rational{$Int}-form $(join((x*"²" for x in VARS[1:NTERMS]), " + ")) + ..."""
end

@testset "Basis Index Iteration" begin
    collectindices(N) = collect(GeomAlg.eachbasisindex(N))

    B = GeomAlg.BasisIndex{0}
    @test collectindices(0) == [B()]

    B = GeomAlg.BasisIndex{1}
    @test collectindices(1) == [B(), B(1)] 

    B = GeomAlg.BasisIndex{2}
    @test collectindices(2) == [B(), B(1), B(2), B(1,2)]

    B = GeomAlg.BasisIndex{3}
    @test collectindices(3) == [
        B(),
        B(1),   B(2),   B(3),
        B(1,2), B(1,3), B(2,3),
        B(1,2,3)
    ]

    B = GeomAlg.BasisIndex{4}
    @test collectindices(4) == [
        B(),
        B(1),     B(2),     B(3),     B(4),
        B(1,2),   B(1,3),   B(1,4),   B(2,3), B(2,4), B(3,4),
        B(1,2,3), B(1,2,4), B(1,3,4), B(2,3,4),
        B(1,2,3,4)
    ]

    B = GeomAlg.BasisIndex{5}
    @test collectindices(5) == [
        B(),
        B(1),       B(2),       B(3),       B(4),       B(5),
        B(1,2),     B(1,3),     B(1,4),     B(1,5),     B(2,3),   B(2,4),   B(2,5),   B(3,4),   B(3,5),   B(4,5),
        B(1,2,3),   B(1,2,4),   B(1,2,5),   B(1,3,4),   B(1,3,5), B(1,4,5), B(2,3,4), B(2,3,5), B(2,4,5), B(3,4,5),
        B(1,2,3,4), B(1,2,3,5), B(1,2,4,5), B(1,3,4,5), B(2,3,4,5),
        B(1,2,3,4,5)
    ]
end
