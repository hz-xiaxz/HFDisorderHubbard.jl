using Random
using LinearAlgebra

@doc raw"""
3D Anderson-Hubbard Model. Hamiltonian see index part.
---------------
    HubbardPara
* `t` : hopping amplitude
* `U` : Hubbard interaction
* `W` : disorder strength, onsite disorder energies are drawn from a Gaussian distribution in the interval `[-W/2, W/2]`
* `n_up` : `⟨n_{i↑}⟩ = ∑_{α} |⟨i↑|α⟩|^2`
* `n_down` : `⟨n_{i↓}⟩ = ∑_{α} |⟨i↓|α ⟩|^2` [1]
* `S_up` : `⟨S_{i}^+⟩ = ∑_{α} |⟨0|c_{i↓}^† c_{i↑}⟩|α⟩|^2`
* `S_down` : conjugation of `S_up`

[1] F.Fazileh et al. 2006 Physical Review Letters
"""
struct HubbardPara
    t::Float64
    U::Float64
    W::Float64
    n_up::Float64
    n_down::Float64
    S_up::Complex{Float64}
    S_down::Complex{Float64}
    omega::Vector{Float64}
    function HubbardPara(;
        t::Float64,
        U::Float64,
        W::Float64,
        n_up::Float64,
        n_down::Float64,
        S_up::Complex{Float64},
        omega::Vector{Float64},
    )
        S_down = conj(S_up)
        return new(t, U, W, n_up, n_down, S_up, S_down, omega)
    end
end


function getHmat(lat::CubicLattice, para::HubbardPara)
    N = lat.N
    @assert length(para.omega) == N
    # H is a complex-valued Hermitian matrix
    H = zeros(Complex{Float64}, 2N, 2N)
    for i = 1:N
        # tunneling
        for j in lat.neigh[i]
            # spin-up cases, not necessarily repeat the same for H[j, i]
            H[i, j] += -para.t # += is necessary
            # spin-down cases
            H[i+N, j+N] += -para.t
        end
        # hubbard interaction
        # n_{iσ} parts
        H[i, i] += para.U * para.n_down
        H[i+N, i+N] += para.U * para.n_up
        # S^{+} S^{-} parts
        # S^+
        H[i, i+N] += -para.U * para.S_down
        # S^-
        H[i+N, i] += -para.U * para.S_up
        # onsite disorder part
        H[i, i] += para.omega[i]
        H[i+N, i+N] += para.omega[i]
        # here we set ω_{i↑} = ω_{i↓} to conserve the rotation symmetry, which is not necessarily true in magnetic ordered state. Ref Pezzoli, Becca 2010
    end
    return H
end

"""
    UnitaryDecomp(lat::CubicLattice, para::HubbardPara)

Use Schur Decomposition to get `H = U^† T U` , where `T` is generally `T` is upper triangular but as `H` is Hermitian, `T` is diagonal (Note the `U` is not the same as the definition in Schur Decomposition).
"""
function UnitaryDecomp(lat::CubicLattice, para::HubbardPara)
    # maybe turn to ArnolidiMethod if bottlenecked
    H = getHmat(lat, para)
    F = schur(H)
    eigenvalues = diag(F.Schur)
    U = F.vectors
    return eigenvalues, U'
end
