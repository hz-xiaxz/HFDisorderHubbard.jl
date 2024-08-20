using HFDisorderHubbard
using Test
using Counters

function compare_tuples(A::Tuple, B::Tuple)
    return counter(collect(A)) == counter(collect(B))
end

@testset "cubic_lattice" begin
    C2 = CubicLattice(2)
    @test C2.N == 8
    @test C2.neigh[1] isa Tuple
    @test compare_tuples(C2.neigh[1], (2, 2, 3, 3, 5, 5))
    @test compare_tuples(C2.neigh[2], (1, 1, 4, 4, 6, 6))
    @test compare_tuples(C2.neigh[8], (7, 7, 4, 4, 6, 6))

    C3 = CubicLattice(3)
    @test C3.N == 27
    @test C3.neigh[1] isa Tuple
    @test compare_tuples(C3.neigh[1], (2, 3, 4, 7, 10, 19))
    @test compare_tuples(C3.neigh[2], (1, 3, 11, 20, 5, 8))
    @test compare_tuples(C3.neigh[27], (26, 25, 24, 21, 18, 9))
end
