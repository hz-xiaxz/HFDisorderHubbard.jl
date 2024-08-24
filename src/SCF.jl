"""
    getNMean(U::AbstractMatrix)

calculate `⟨n_{i}⟩ = ∑_{α} |⟨i|α ⟩|^2` where i ∈ 1:2L, α ∈ 1:N. 1:L labels spin-up operators, L+1:2L labels spin-down operators.

return: n_up, n_down
"""
function getNMean(U::AbstractMatrix)
    # eigenstates = U D U^†
    # φ_α^† = ∑_{I} U_{Iα} d_I^† only count U_{iα}
    # <n_i> = ∑_{α=1..N} |U_{iα}|^2
    N = size(U, 2)
    return [sum(abs2, U[i, :]) for i = 1:N], [sum(abs2, U[i, :]) for i = N+1:2N]
end

"""
    getSupMean(U::AbstractMatrix)
return: S_up
"""
function getSupMean(U::AbstractMatrix)
    # don't introduce ⟨β|S_i^+|α⟩ where |α⟩ and |β⟩ are different orbitals, the coeffecient on each is unapproachable
    # just calculate ⟨S_i^+⟩ = ∑_{α} (U_{α_i}^* × U_{i+N,α})
    N = size(U, 2)
    # Sup is of N long
    return [sum(@. conj(U[i, :]) * U[i+N, :]) for i = 1:N] # TODO:check this!
end


"""
    init(rng::AbstractRNG, lat::CubicLattice)
    initialize SCF guess
"""
function init(rng::AbstractRNG, lat::CubicLattice)
    N = lat.N
    n_up = rand(rng, N)
    left = N - sum(n_up)
    weights = rand(rng, N)
    n_down = left * weights / sum(weights)
    # init with the constraint of particle number
    # random init S_up with amplitude restricted to less than 1/2
    # use polar coordinate
    r = rand(rng, N) * 1 / 2
    θ = rand(rng, N) * 2π
    Sx = @. r * cos(θ)
    Sy = @. r * sin(θ)
    S_up = Sx + Sy * im
    return SCFdata(; n_up = n_up, n_down = n_down, S_up = S_up)
end

"""
    step!(data::SCFdata, lat::CubicLattice, para::HubbardPara)
    step on SCF, if converge, returns true.
"""
function step!(data::SCFdata, lat::CubicLattice, para::HubbardPara)
    H = getHmat(lat, para, data)
    eigenvalues, U = UnitaryDecomp(H)
    n_up, n_down = getNMean(U)
    S_up = getSupMean(U)
    flag = checkConverge(data, n_up, n_down, S_up)
    data.n_up = n_up
    data.n_down = n_down
    data.S_up = S_up
    return flag
end

"""
    checkConverge(data::SCFdata, n_up::Vector{Float64}, n_down::Vector{Float64}, S_up::Vector{Complex{Float64}})
check if SCF converges. Ref: Inui and Littlewood (1991)
"""
function checkConverge(
    data::SCFdata,
    n_up::Vector{Float64},
    n_down::Vector{Float64},
    S_up::Vector{Complex{Float64}},
    tol::Float64 = 1e-6,
)
    # check if the difference is smaller than 1e-6
    N = length(n_up)
    Sx = real.(S_up)
    Sy = imag.(S_up)
    Sz = @. n_up - n_down
    old_Sx = real.(data.S_up)
    old_Sy = imag.(data.S_up)
    old_Sz = @. data.n_up - data.n_down
    value = √(
        1 / N * sum(
            @.(
                abs2(n_up + n_down - data.n_up - data.n_down) +
                abs2(Sx - old_Sx) +
                abs2(Sy - old_Sy) +
                abs2(Sz - old_Sz)
            )
        ),
    )
    return value < tol
end
