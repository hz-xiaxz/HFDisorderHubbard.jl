using HFDisorderHubbard
using Test
using Random

@testset "SCF_init" begin
    rng = Random.Xoshiro(42)
    lat = CubicLattice(2)
    init_guess = HFDisorderHubbard.init(rng, lat)
    @test length(init_guess.n_up) == 8
    @test length(init_guess.n_down) == 8
    @test length(init_guess.S_up) == 8
    @test allequal(init_guess.n_up .<= 1.0)
    @test allequal(init_guess.n_up .>= 0.0)
    @test allequal(init_guess.n_down .<= 1.0)
    @test allequal(init_guess.n_down .>= 0.0)
    @test allequal(abs.(init_guess.S_up) .<= 0.5)
    @test allequal(real.(init_guess.S_up) .<= 0.5)
    @test allequal(real.(init_guess.S_up) .>= -0.5)
    @test allequal(imag.(init_guess.S_up) .<= 0.5)
    @test allequal(imag.(init_guess.S_up) .>= -0.5)
end
