```@meta
CurrentModule = HFDisorderHubbard
```

# HFDisorderHubbard

Documentation for [HFDisorderHubbard](https://github.com/hz-xiaxz/HFDisorderHubbard.jl).

Implement Hatree Fock study of 3D Hubbard model with onsite disorder.

## Methodology

First, given the Hartree-Fock mean field Hamiltonian [1] with onsite disorder

$$
\begin{aligned}
H_{HF} &= -t \sum_{\langle i,j \rangle, \sigma} \left(  c_{i\sigma}^\dagger c_{j\sigma}+h.c.\right)\\
  &+ U \sum_i \left( n_{i\uparrow} \left< n_{i \downarrow} \right> + \left< n_{i\uparrow}\right> n_{i \downarrow} - \left<S_i^+ \right> S_i^-  - \left< S_i^-\right> S_i^+ \right) \\
& -U \sum_i \left(   \left< n_{i\uparrow}\right> \left< n_{i \downarrow} \right> - \left<S_i^+ \right>\left<S_i^- \right> \right) \\
& + \sum_{i} \left( \omega_i n_i \right)
\end{aligned}
$$

where $\omega_i$ is drawn from a Gaussian Distribution $N(0,W)$. One can easily diagonalize this Hamiltonian with the fact that $n_{i\sigma}=c^\dagger_{i\sigma}c_{i\sigma}$, $S^+_{i} = c^\dagger_{i\uparrow}c_{i\downarrow}$ and $S^-_{i} = c^\dagger_{i\downarrow}c_{i\uparrow}$ *Note: This matrix is Hermitian since* $\left( S^+ \right)^\dagger = S^- $

Second, generate a random initial guess for $\left<n_{i\uparrow} \right>$ , $\left<n_{i\downarrow} \right>$, $\left<S_i^+ \right>$ , $\left<S_i^- \right>$. Diagonalize the Hamiltonian and get the new $\left<n_{i\uparrow} \right>$ , $\left<n_{i\downarrow} \right>$, $\left<S_i^+ \right>$ , $\left<S_i^- \right>$. Repeat until convergence.

Third, check convergence...

Fourth, draw new onsite disorder and get average.

### Calculation of average

In case of half filling, we have $2N$ orbitals and $N$ electrons. The Fermi energy is the energy of the $N$th orbital.

Technically, we treat $\displaystyle \left<n_{i↑} \right> = \sum_{x}\left<\Phi_0|n_{i\uparrow} \left|x\left>\right<x \right|\Phi_0\right>$ where $\displaystyle \left<x|\Phi_0 \right>= \det \left( U_{R_j, \alpha} \right) $ is the slater determinant of certain configuration. $n_{i\uparrow}$ will be 1 or 0, thus $\displaystyle \left<n_{i\uparrow} \right> = \sum_{R_j \text{ having } n_{i\uparrow}} \left|   \det\left( U_{R_j,\alpha} \right)\right|^2 $. $R_j$ are drawn from $\begin{pmatrix} 2L-1\\ N-1\end{pmatrix}$ configurations, its computational cost is unaffordable, but can be achieved by Monte Carlo Sampling.

Let's see if we can make it more naive by introducing single particle approximation. Let's not construct a many-body state $\left| \Phi_0 \right>$. Regard each $\left|\alpha\right> = \phi_\alpha^\dagger \left|0\right>$ be orthogonal orbitals, they are one-particle orbitals. So we can calculate each $\left<n_{i\uparrow} \right>$ by $\displaystyle \left<n_{i\uparrow} \right> = \sum_{\alpha, \epsilon_{\alpha }<\epsilon_{F}} \left| \left   <i\uparrow |\alpha \right>\right|^2$.

Averages are defined as $\displaystyle \left< n_{i\sigma} \right> = \sum_{\alpha, \epsilon_{\alpha }<\epsilon_{F}} \left| \left   <i\sigma |\alpha \right>\right|^2$ so are $\left<S^+ \right>$ and $\left<S^- \right>$ [2] In our case orbitals can't be easily divided into spin-up and spin-down ones, thus the formalism is slightly different.

> TODO: check how good is the single-particle approximation

References:
> [1] M. Inui* and P. B. Littlewood (1991)
> [2] F. Fazileh et al. (2006)

## Contributors

```@raw html
<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
```
