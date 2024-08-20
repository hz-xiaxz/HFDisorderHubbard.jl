using ArnoldiMethod

@doc raw"""
3D Anderson-Hubbard Model

``
\begin{aligned}
H_{HF} &= -t \sum_{\langle i,j \rangle, \sigma} \left(  c_{i\sigma}^\dagger c_{j\sigma}+h.c.\right)\\
  &+ U \sum_i \left( n_{i\uparrow} \left< n_{i \downarrow} \right> + \left< n_{i\uparrow}\right> n_{i \downarrow} - \left<S_i^+ \right> S_i^-  - \left< S_i^-\right> S_i^+ \right) \\
& -U \sum_i \left(   \left< n_{i\uparrow}\right> \left< n_{i \downarrow} \right> - \left<S_i^+ \right>\left<S_i^- \right> \right) \\
& + \sum_{i} \left( \omega_i n_i \right)
\end{aligned}
``
---------------
    HubbardPara
* `t` : hopping amplitude
* `U` : Hubbard interaction
* `W` : disorder strength, onsite disorder energies are drawn from a Gaussian distribution in the interval `[-W/2, W/2]`
* `n_up` : average number of up spins
* `n_down` : average number of down spins
* `S_up` : average magnetization of up spins
* `S_down` : average magnetization of down spins

"""
struct HubbardPara
    t::Float64
    U::Float64
    W::Float64
    n_up::Float64
    n_down::Float64
    S_up::Float64
    S_down::Float64
end
