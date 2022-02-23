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
