"""
    getNMean(eigenstates::AbstractMatrix)

calculate `⟨n_{i}⟩ = ∑_{α} |⟨i|α ⟩|^2` where i ∈ 1:2L, α ∈ 1:N. 1:L labels spin-up operators, L+1:2L labels spin-down operators.
"""
function getNMean(eigenstates::AbstractMatrix)
    # eigenstates = U D U^†
    # φ_α^† = ∑_{I} U_{Iα} d_I^† only count U_{iα}
    # <n_i> = ∑_{α=1..N} |U_{iα}|^2
    return [sum(abs2, eigenstates[i, :]) for i = 1:2N]
end
