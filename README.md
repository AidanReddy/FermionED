# FermionED

This repository contains a .jl file, FermionED.jl, with a few functions for constructing a numerical many-fermion Hamiltonian in the Julia programming language.

The code first constructs a many-body Fock space from a given set of single-particle states. Particle number, spin, and center-of-mass moemntum conservation are exploited as relevant. Each single particle state is assigned an integer index "i" valued between 1 and (total number of single-particle states). Each Fock state $\ket{I}$ is assigned an integer $I = \sum_{i \in occ} 2^i$ where "occ" is a set of occupied single particle states. The $i^{th}$ entry (zero-indexed, from right to left) of the bit string corresponding to the integer $I$ is 1 if $i \in occ$ and 0 if not. (Note that one could also 0-index the single-particle states, which is perhaps the more natural choice but not the one adopted here. Here, the single-particle states are 1-indexed.) For spinful systems with $N$, spin-up states ar given the indices 1 through N, and spin-down states are given the indices (N+1) through 2N. All manipulations use the functions "C" and "CDag" which implement second quantization entirely through arithmetic operations on the integers I. This methods allows for efficient construction of the many-body Hamiltonian and calculation of observables with a low memory footprint. 

The functions provided here only handle one-and two-body terms, but the "C" and "CDag" functions can be straightforwardly used to handle n-body terms. Explicitly, the functions are designed to handle Hamiltonians of the form
$H = T + V$ where $T = \sum_{ij,s_i s_j} t_{is_i,js_j} c^{\dagger}_{i,s_i} c_{j,s_j}$ is a ony-body term and $V = \frac{1}{2} \sum_{ijkl,s_js_js_ks_l} V_{isi,jsj,ks_k,ls_l} c^{\dagger}_{i,s_i} c^{\dagger}_{j,s_j} c_{k,s_k} c_{l,s_l}$ is a two-body term. The matrix elements $t_{i,s_i,j,s_j}$ and $V_{is_i,js_j,ks_k,ls_l}$ must be calculated independently for the physical Hamiltonian of interest.

The functions contained here build a many-body Hamiltonian in the SparseMatrixCSC format from the SparseArrays.jl package. This can be diagonalized using other packages such as KrylovKit.jl ( "eigsolve" function) or Arpack.jl ("eigs" function). 

I originally developed for the calculations published in Ref. [2], and have used it in several subsequent publications. I benefitted from Ref. [2].

This work is supported by the Air Force Office of Scientific Research (AFOSR) under Award No. FA9550-22-1-0432 and the Simons Investigator award from the Simons Foundation.

Dependencies: Combinatorics.jl, SparseArrays.jl

[1] Reddy, A.P., Alsallom, F., Zhang, Y., Devakul, T. and Fu, L., 2023. Fractional quantum anomalous Hall states in twisted bilayer MoTe2 and WSe2. Physical Review B, 108(8), p.085117.

[2] Jafari, S.A., 2008. Introduction to Hubbard model and exact diagonalization. arXiv preprint arXiv:0807.4878.