using HFDisorderHubbard
using Test
using Random
using LinearAlgebra

lat = HFDisorderHubbard.CubicLattice(2)
rng = Random.Xoshiro(42)
N = lat.N
omega = randn(rng, N)
para = HFDisorderHubbard.HubbardPara(
    U = 1.0,
    t = 1.0,
    W = 1.0,
    n_up = 0.6,
    n_down = 0.4,
    S_up = 1.0 - 1.0im,
    omega = omega,
)
H = HFDisorderHubbard.getHmat(lat, para)

@testset "Hubbard" begin
    @test allequal(
        omega .== [
            -0.36335748145177754,
            0.2517372155742292,
            -0.31498797116895605,
            -0.31125240132442067,
            0.8163067649323273,
            0.47673837983187795,
            -0.8595553820616212,
            -1.4692882055065464,
        ],
    )
    @test H[1, 1] == omega[1] + para.U * para.n_down
    @test H[1+N, 1+N] == omega[1] + para.U * para.n_up
    @test H[1, 2] == -2para.t # 2 is due to periodical boundary condition
    @test H[1+N, 2+N] == -2para.t
    @test H[1, 1+N] == -para.U * para.S_down
    @test H[1+N, 1] == conj(H[1, 1+N])
    @test allequal(H .== H')
end


@testset "Unitary Decomposition" begin
    ev, es = HFDisorderHubbard.UnitaryDecomp(lat, para)
    # test es is normalized
    # es is somehow complex valued
    @test isapprox(es * es', I, atol = 1e-8)
    @test allequal(@. ev isa Real)
end
