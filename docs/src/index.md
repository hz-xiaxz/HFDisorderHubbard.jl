```@meta
CurrentModule = HFDisorderHubbard
```

# HFDisorderHubbard

Documentation for [HFDisorderHubbard](https://github.com/hz-xiaxz/HFDisorderHubbard.jl).

Implement Hatree Fock study of 3D Hubbard model with onsite disorder.

## Methodology [1]

1. Given the Hartree-Fock mean field Hamiltonian with onsite disorder

$$
\begin{aligned}
H_{HF} &= -t \sum_{\langle i,j \rangle, \sigma} \left(  c_{i\sigma}^\dagger c_{j\sigma}+h.c.\right)\\
  &+ U \sum_i \left( n_{i\uparrow} \left< n_{i \downarrow} \right> + \left< n_{i\uparrow}\right> n_{i \downarrow} - \left<S_i^+ \right> S_i^-  - \left< S_i^-\right> S_i^+ \right) \\
& -U \sum_i \left(   \left< n_{i\uparrow}\right> \left< n_{i \downarrow} \right> - \left<S_i^+ \right>\left<S_i^- \right> \right) \\
& + \sum_{i} \left( \omega_i n_i \right)
\end{aligned}
$$
where $\omega_i$ is drawn from a Gaussian Distribution $N(0,W)$
2. First generate a random initial guess for $\left<n_{i\uparrow} \right>$ , $\left<n_{i\downarrow} \right>$, $\left<S_i^+ \right>$ , $\left<S_i^- \right>$. Diagonalize the Hamiltonian and get the new $\left<n_{i\uparrow} \right>$ , $\left<n_{i\downarrow} \right>$, $\left<S_i^+ \right>$ , $\left<S_i^- \right>$. Repeat until convergence.
3. Check convergence...
4. Draw new onsite disorder and get average.

> [1] M. Inui* and P. B. Littlewood (1991)
>
## Contributors

```@raw html
<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
```
