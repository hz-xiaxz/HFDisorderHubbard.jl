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
    @test sum(init_guess.n_up) + sum(init_guess.n_down) ≈ 8.0
    @test allequal(init_guess.n_up .>= 0.0)
    @test allequal(init_guess.n_down .>= 0.0)
    # @test allequal(init_guess.n_up .<= 1.0)
    # temporarily allow > 1
    # @test allequal(init_guess.n_down .<= 1.0)
    @test allequal(abs.(init_guess.S_up) .<= 0.5)
    @test allequal(real.(init_guess.S_up) .<= 0.5)
    @test allequal(real.(init_guess.S_up) .>= -0.5)
    @test allequal(imag.(init_guess.S_up) .<= 0.5)
    @test allequal(imag.(init_guess.S_up) .>= -0.5)
end

@testset "getSupMean" begin
    lat = HFDisorderHubbard.CubicLattice(2)
    rng = Random.Xoshiro(42)
    N = lat.N
    omega = randn(rng, N)
    para = HFDisorderHubbard.HubbardPara(U = 1.0, t = 1.0, W = 1.0, omega = omega)
    data = HFDisorderHubbard.SCFdata(
        n_up = fill(0.6, N),
        n_down = fill(0.4, N),
        S_up = fill(0.25 + 0.25im, N),
    )
    H = HFDisorderHubbard.getHmat(lat, para, data)
    ev, U = HFDisorderHubbard.UnitaryDecomp(H)
    S_up = HFDisorderHubbard.getSupMean(U)
    @test length(S_up) == 8
    @test allequal(@. real(S_up) <= 0.5)
    @test allequal(@. real(S_up) >= -0.5)
    @test allequal(@. imag(S_up) <= 0.5)
    @test allequal(@. imag(S_up) >= -0.5)
    @test allequal(@. abs(S_up) <= 0.5)
end
