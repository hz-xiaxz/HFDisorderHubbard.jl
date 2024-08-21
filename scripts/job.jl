using DrWatson
using HFDisorderHubbard
using Random

function run()
    rng = Random.Xoshiro(42)
    L = 2
    lat = CubicLattice(L)
    data = init(rng, lat)
    iter = 10
    steps = 10
    for i = 1:iter
        t = 1.0
        U = 1.0
        W = 1.0
        para = HubbardPara(t = t, U = U, W = W, omega = randn(rng, L^3))
        params = @strdict L U W i
        for _ = 1:steps
            flag = step!(data, lat, para)
            if flag
                break
            end
            dicted_data =
                Dict("n_up" => data.n_up, "n_down" => data.n_down, "S_up" => data.S_up)
            @tagsave(datadir("test", savename(params, "jld2")), dicted_data)
        end
    end
end

run()
firstsim = readdir(datadir("test"))[1]
@show wload(datadir("test", firstsim))
