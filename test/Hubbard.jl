using HFDisorderHubbard
using Test
using Random
using LinearAlgebra

lat = HFDisorderHubbard.CubicLattice(2)
rng = Random.Xoshiro(42)
N = lat.N
omega = randn(rng, N)
para = HFDisorderHubbard.HubbardPara(U = 1.0, t = 1.0, W = 1.0, omega = omega)
data = HFDisorderHubbard.SCFdata(
    n_up = fill(0.6, N),
    n_down = fill(0.4, N),
    S_up = fill(1.0 + 1.0im, N),
)
H = HFDisorderHubbard.getHmat(lat, para, data)

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
    @test H[1, 1] == omega[1] + para.U * data.n_down[1]
    @test H[1+N, 1+N] == omega[1] + para.U * data.n_up[1]
    @test H[1, 2] == -2para.t # 2 is due to periodical boundary condition
    @test H[1+N, 2+N] == -2para.t
    @test H[1, 1+N] == -para.U * conj(data.S_up[1])
    @test H[1+N, 1] == conj(H[1, 1+N])
    @test allequal(H .== H')
end


@testset "Unitary Decomposition" begin
    ev, U = HFDisorderHubbard.UnitaryDecomp(H)
    @test allequal(ev .== sort(ev))
    @test size(U) == (2N, N)
    @inbounds for i = 1:N
        @test sum(abs2, U[:, i]) ≈ 1.0
    end
end
