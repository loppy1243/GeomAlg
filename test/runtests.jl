using Test
using GeomAlg

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
