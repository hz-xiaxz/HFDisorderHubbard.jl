using DrWatson
@quickactivate "HFDisorderHubbard"
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
        stored = Dict{String,SCFdata}()
        for j = 1:steps
            flag = step!(data, lat, para)
            stored[string(j)] = deepcopy(data)
            if flag
                @show "converged"
                break
            end
        end
        @tagsave(datadir("test", savename(params, "jld2")), stored)
    end
end

run()
using DataFrames
df = collect_results(datadir("test"))

# TODO: work_load and visualize
