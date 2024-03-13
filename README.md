# FermionED

This repository contains a .jl file, FermionED, with a few functions for constructing a numerical many-fermion Hamiltonian in the Julia language.

The basic structure of the code is to construct a many-body Fock space from a given set of single-particle states. Particle number, spin, and center-of-mass moemntum conservation are exploited as relevant. Each single particle state is assigned an integer index "i" valued between 1 and (total number of single-particle states). Each Fock state $\ket{I}$ is assigned an integer $I = \sum_{i \in occ} 2^i$ where "occ" is a set of occupied single particle states. The $i^{th}$ entry (zero-indexed) of the bit string corresponding to the integer $I$ is 1 if $i \in occ$ and 0 if not. (Note that one could also 0-index the single-particle states, which is perhaps the more natural choice but not the one adopted here. Here, the single-particle states are 1-indexed.) For spinful systems with $N$, spin-up states ar given the indices 1 through N, and spin-down states are given the indices (N+1) through 2N. All manipulations use the functions "C" and "CDag" which implement second quantization entirely through arithmetic operations on the integers I. This methods allows for efficient construction of the many-body Hamiltonian and calculation of observables with a low memory footprint. 

The functions provided here only handle one-and two-body terms, but the "C" and "CDag" functions can be straightforwardly used to handle n-body terms. Explicitly, the functions are designed to handle Hamiltonians of the form

$ H = \sum_{ij,s_is_j} t_{i,s_i,j,s_j} c^{\dagger}_{i,s_i} c_{j,s_j} + \frac{1}{2} \sum_{ijkl,s_js_js_ks_l} V_{i,si,j,sj,k,si,l,sj} c^{\dagger}_{i,s_i} c^{\dagger}_{j,s_j} c_{k,sk} c_{l,sl}$.

The matrix elements $t_{i,s_i,j,s_j}$ and $V_{i,s_i,j,s_j,k,s_k,l,s_l}$ must be calculated externally for the physical Hamiltonian of interest.

The functions contained here build a many-body Hamiltonian in the SparseMatrixCSC format from the SparseArrays.jl package. This can be diagonalized using other packages such as KrylovKit.jl ( "eigsolve" function) or Arpack.jl ("eigs" function). 

I originally developed for the calculations published in Ref. [2], and have used it in several subsequent publications. I benefitted from Ref. [2].

This work is supported by the Air Force Office of Scientific Research (AFOSR) under Award No. FA9550-22-1-0432 and the Simons Investigator award from the Simons Foundation.

Dependencies: Combinatorics.jl, SparseArrays.jl

[1] Reddy, A.P., Alsallom, F., Zhang, Y., Devakul, T. and Fu, L., 2023. Fractional quantum anomalous Hall states in twisted bilayer MoTe2 and WSe2. Physical Review B, 108(8), p.085117.

[2] Jafari, S.A., 2008. Introduction to Hubbard model and exact diagonalization. arXiv preprint arXiv:0807.4878.