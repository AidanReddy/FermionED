# FermionED

This repository contains a .jl file, FermionED.jl, with a few functions for constructing a numerical many-fermion Hamiltonian in the Julia programming language.

The "exact diagonalization" method is a non-perturbative numerical method for studying quantum many-body systems. Quite simply, one constructs and (partially) diagonalizes a matrix representation of a many-body Hamiltonian. It relies on the presence of a finite-dimensional many-body Hilbert space. On a lattice system (e.g. one-band Hubbard model) this comes for free on a finite-system. In a continuum system, the many-body Hilbert space is infinite dimensional even on a finite size system. In this case, the method is variational, rather than "exact", and necessarily involves projection onto the Fock space defined by a set of single-particle states that are energetically isolated (e.g. a Bloch band or Landau level). This method has historically found great use in the study of the fractional quantum Hall effect partially-filled Landau levels and, more recently, the fractional quantum anomalous Hall effect in partially-filled Chern bands.

While this code if far from a complete functioning software package, I hope that it may find use as a starting point or reference for others developing their own exact diagonalization codes for many-fermion systems.

The code first constructs a many-body Fock space from a given set of single-particle states. Particle number, spin, and center-of-mass moemntum conservation are exploited as relevant. Each single particle state is assigned an integer index "i" valued between 1 and (total number of single-particle states). Each Fock state $\ket{I}$ is assigned an integer $I = \sum_{i \in occ} 2^i$ where "occ" is a set of occupied single particle states. The $i^{th}$ entry (zero-indexed, from right to left) of the bit string corresponding to the integer $I$ is 1 if $i \in occ$ and 0 if not. (Note that one could also 0-index the single-particle states, which is perhaps the more natural choice but not the one adopted here. Here, the single-particle states are 1-indexed.) For spinful systems with $N$, spin-up states are given the indices 1 through N and spin-down states are given the indices (N+1) through 2N. All manipulations use the functions "C" and "CDag" which implement second quantization entirely through arithmetic operations on the integers $I$. This methods allows for efficient construction of the many-body Hamiltonian and calculation of observables. 

The functions provided here only handle one-body and two-body terms. Explicitly, the functions are designed to handle Hamiltonians of the form
$H = T + V$ where $T = \sum_{ij,s_i s_j} t_{is_i,js_j} c_{i,s_i}^{\dagger} c_{j,s_j}$ is a ony-body term and $V = \frac{1}{2} \sum_{ijkl,s_js_js_ks_l} V_{is_i,js_j,ks_k,ls_l} c_{i,s_i}^{\dagger} c_{j,s_j}^{\dagger} c_{k,s_k} c_{l,s_l}$ is a two-body term. The matrix elements $t_{is_i,js_j}$ and $V_{is_i,js_j,ks_k,ls_l}$ must be calculated independently for the physical Hamiltonian at hand. The functions "C" and "CDag" functions could be straightforwardly used to handle n-body terms.

The functions contained here build a many-body Hamiltonian in the SparseMatrixCSC format from the SparseArrays.jl package. This sparse matrix can be partially diagonalized using other packages such as KrylovKit.jl ( "eigsolve" function) or Arpack.jl ("eigs" function), or converted to a dense matrix and fully diagonalized using LinearAlgebra.jl ("eigen" function).

I originally developed this code for the calculations published in Ref. [1] and have used it in several subsequent publications. I provide further discussion of the exact diagonalization method as applied to moiré superlattice systems in the Supplementary Material of Ref. [2-7]. In developing this code, I benefitted from Ref. [8]. I gratefully acknowledge helpful conversations from Trithep Devakul, Yang Zhang, Di Luo, and Donna Sheng. Finally, I thank Liang Fu for advising the research projects for which I developed and used this code. 

This work is supported by the Air Force Office of Scientific Research (AFOSR) under Award No. FA9550-22-1-0432 and the Simons Investigator award from the Simons Foundation.

Dependencies: Combinatorics.jl, SparseArrays.jl

[1] Reddy, A.P., Alsallom, F., Zhang, Y., Devakul, T. and Fu, L., 2023. Fractional quantum anomalous Hall states in twisted bilayer MoTe2 and WSe2. Physical Review B, 108(8), p.085117.

[2] Goldman, H., Reddy, A.P., Paul, N. and Fu, L., 2023. Zero-field composite Fermi liquid in twisted semiconductor bilayers. Physical Review Letters, 131(13), p.136501.

[3] Abouelkomsan, A., Reddy, A.P., Fu, L. and Bergholtz, E.J., 2023. Band mixing in the quantum anomalous Hall regime of twisted semiconductor bilayers. arXiv preprint arXiv:2309.16548.

[4] Reddy, A.P. and Fu, L., 2023. Toward a global phase diagram of the fractional quantum anomalous Hall effect. Physical Review B, 108(24), p.245159.

[5] Tan, T., Reddy, A.P., Fu, L. and Devakul, T., 2024. Designing topology and fractionalization in narrow gap semiconductor films via electrostatic engineering. arXiv preprint arXiv:2402.03085.

[6] Sheng, D. N., Reddy, A. P., Abouelkomsan, A., Bergholtz, E. J., & Fu, L. (2024). Quantum anomalous Hall crystal at fractional filling of moir\'e superlattices. arXiv preprint arXiv:2402.17832.

[7] Reddy, A.P., Paul, N., Abouelkomsan, A. and Fu, L., 2024. Non-Abelian fractionalization in topological minibands. arXiv preprint arXiv:2403.00059.

[8] Jafari, S.A., 2008. Introduction to Hubbard model and exact diagonalization. arXiv preprint arXiv:0807.4878.
