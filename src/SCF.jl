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
    return [sum(@. conj(U[:, i]) * U[i+N, :]) for i = 1:N]
end

function init()
    # init disorder
    # init guess
end

"""
    sweep()
update the parameter
"""
function sweep()

end

function register()

end

"""
    checkConverge()
check if SCF converges
"""
function checkConverge() end
