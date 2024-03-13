module ed
# This module contains several functions forming the core of an exact diagonalization calculation for a many-fermion system.

using Combinatorics
using SparseArrays

#FUNCTIONS GENERIC TO ALL MANY-FERMION SYSTEMS WITH A SINGLE TWO-LEVEL CONSERVED "SPIN" DEGREE OF FREEDOM.

#=
Definitions (not exhaustive):

n: total number of single-particle basis states, per spin
T: one-body part of many-body Hamiltonian
V: two-body interaction part of many-body Hamiltonian.
nU: number of "spin-up" fermions
nD: number of "spin-down" fermions
oneBody: Vector of one-body matrix elements, each of which is a tuple of ComplexF64 and a tuple {t_{i,si,j,sj}, {i, si, j, sj}}
twoBody: Vector of two-body matrix elements, each of which is a tuple of a ComplexF64 and a tuple {V_{i,si,j,sj,k,si,l,sj}, {i, si, j, sj, k, sk, l, sl}}
=#

#calculate sign coming from fermionic commutation relations
function calcSign(I::Int64, K::Int64)::Int64
    M = I-(I&(2K-1)) #occupied states to the left of (i,s)
    btwnCnt = count_ones(M) #counts all occupied to the left of (i,s)
    sign = (-1)^btwnCnt
    return sign
end

#annihilation operator
function C(i::Int64, s::Bool, I::Int64, sign::Int64, n::Int64)::Tuple{Int64, Int64}
    if sign == 0 #if sign is 0, the many body state is anihilated so can just proceed
        return 0, 0
    end
    K = (2^i * (s*2^(n) - (s-1))) #integer representation of single particle state
    if (I & K) != K # if (i,s) was already unoccupied, then anihilate many body state
        return 0, 0
    else
        L = I - K  #anihilate (i,s) in I if (i,s) is occupied in I
        sign *= calcSign(I, K)
        return L, sign
    end
end

#creation operator
function CDag(i::Int64, s::Bool, I::Int64, sign::Int64, n::Int64)::Tuple{Int64, Int64}
    if sign == 0 #if sign is 0, the many body state is anihilated so can just proceed
        return 0, 0
    end
    K = (2^i * (s*2^(n) - (s-1))) #integer representation of single particle state
    if I & K == K # if state is already occupied, anihilate many body state
        return 0, 0
    else
        L = I + K #create (i,s) in I
        sign *= calcSign(I, K)
        return L, sign
    end
end

#one-body part of hamiltonian
function calcT(oneBody::Vector{Tuple{ComplexF64, Tuple{Int64, Bool, Int64, Bool}}}, States::Vector{Int64}, n::Int64)::SparseMatrixCSC{ComplexF64, Int64}
    rows = Vector{Int64}()
    cols = Vector{Int64}()
    vals = Vector{ComplexF64}()
    for value in oneBody
        (t, (i, si, j, sj)) = value
        for IInd in eachindex(States)
            I = States[IInd]
            sign = 1
            J, sign = C(j, sj, I, sign, n)
            J, sign = CDag(i, si, J, sign, n)
            if sign == 0 #if State is anihilated
                continue
            end
            JInd = searchsortedfirst(States,J) #REQUIRES THAT States BE SORTED. see http://www.jlhub.com/julia/manual/en/function/searchsortedfirst
            if JInd == (length(States)+1) || States[JInd]!=J
                continue
            end
            push!(rows, JInd)
            push!(cols, IInd)
            push!(vals, sign*t)
        end
    end
    T = sparse(rows, cols, vals, length(States),length(States))
    return T
end

#two-body part of Hamiltonian
function calcV(twoBody::Vector{Tuple{ComplexF64, Tuple{Int64, Bool, Int64, Bool,Int64, Bool, Int64, Bool}}}, States::Vector{Int64}, n::Int64)::SparseMatrixCSC{ComplexF64, Int64}
    Rows = [Vector{Int64}() for _ in 1:Threads.nthreads()]
    Cols = [Vector{Int64}() for _ in 1:Threads.nthreads()]
    Vals = [Vector{ComplexF64}() for _ in 1:Threads.nthreads()]
    Threads.@threads for i in eachindex(twoBody)
        (v, (i, si, j, sj, l, sl, k, sk)) = twoBody[i]
        for IInd in eachindex(States)
            I = States[IInd]
            sign = 1
            J, sign = C(k, sk, I, sign, n)
            J, sign = C(l, sl, J, sign, n)
            J, sign = CDag(j, sj, J, sign, n)
            J, sign = CDag(i, si, J, sign, n)
            if sign == 0 #if State is anihilated
                continue
            end
            JInd = searchsortedfirst(States,J)
            if JInd == (length(States)+1) || States[JInd]!=J
                continue
            end
            push!(Rows[Threads.threadid()], JInd)
            push!(Cols[Threads.threadid()], IInd)
            push!(Vals[Threads.threadid()], sign*v/2) # factor of 1/2 is necessary to not double count interactions
        end
    end
    rows=reduce(vcat,Rows)
    cols=reduce(vcat,Cols)
    vals=reduce(vcat,Vals)
    #FREE UP RAM
    Rows=nothing
    Cols=nothing
    Vals=nothing
    GC.gc() #garbage collect
    V = sparse(rows, cols, vals, length(States),length(States))
    rows=nothing
    cols=nothing
    vals=nothing
    GC.gc() 
    return V
end


#generate spin projected fock state
function calcFockStates(ne::Int64, n::Int64)::Vector{Int64}
    States = Vector{Int64}()
    indVec = collect(1:n)
    occLists = Combinatorics.combinations(indVec, ne)
    for occList in occLists
        State = 0
        for i in occList
            State += 2^i
        end
        push!(States, State)
    end 
    sort!(States)
    return States
end

#occupation <c+_i c_i> numbers of single particle i
function calcocc(State::Vector{ComplexF64},n::Int64,Basis::Vector{Int64})::Array{Float64, 2}
    occ=zeros(Float64, (n,2))
    for i in 1:n, s in 0:1
        for IInd in eachindex(Basis)
            sign=1
            I=Basis[IInd]
            J,sign=C(i,Bool(s),I,sign,n)
            if sign != 0
                occ[i,s+1] += abs(State[IInd])^2
            end
        end
    end
    return occ
end

#one body density matrix, general
#D2[i,si,j,sj] = < c+_i,si c_j,sj >
function calcD1(State::Vector{ComplexF64},n::Int64, Basis::Vector{Int64})::Array{ComplexF64, 4}
    D1Chunks=[zeros(ComplexF64, (n,2,n,2)) for _ in 1:Threads.nthreads()]
    Threads.@threads for IInd in eachindex(Basis)
        for i in 1:n, si in 0:1, j in 1:n, sj in 0:1
            sign=1
            I=Basis[IInd]
            J,sign=C(j,Bool(sj),I,sign,n)
            J,sign=CDag(i,Bool(si),J,sign,n)
            JInd = searchsortedfirst(Basis,J)
            if isnothing(JInd)==false
                D1Chunks[Threads.threadid()][i,si+1,j,sj+1]+=conj(State[JInd])*State[IInd]*sign
            end
        end
    end
    D1=sum(D1Chunks)
    return D1
end


#two body density matrix, general
#D2[i,si,j,sj,k,sk,l,sl] = < c+_i,si c+_j,sj c_k,sk c_l,sl >
#note that in most cases, it is much more efficient to explot symmetries in calculating the two-body density matrix (see e.g. calcD2MomSpinConsOneBand below)
function calcD2(State::Vector{ComplexF64},n::Int64, Basis::Vector{Int64})::Array{ComplexF64, 8}
    println("num threads calcD2: $(Threads.nthreads())")
    D2Chunks=[zeros(ComplexF64, (n,2,n,2,n,2,n,2)) for _ in 1:Threads.nthreads()]
    Threads.@threads for IInd in eachindex(Basis)
        I=Basis[IInd]
        for i in 1:n, si in 0:1, j in 1:n, sj in 0:1, k in 1:n, sk in 0:1, l in 1:n, sl in 0:1
            sign=1
            J,sign=C(l,Bool(sl),I,sign,n)
            J,sign=C(k,Bool(sk),J,sign,n)
            J,sign=CDag(j,Bool(sj),J,sign,n)
            J,sign=CDag(i,Bool(si),J,sign,n)
            JInd = searchsortedfirst(Basis,J)
            if isnothing(JInd)==false
                D2Chunks[Threads.threadid()][i,si+1,j,sj+1,k,sk+1,l,sl+1]+=conj(State[JInd])*State[IInd]*sign
            end
        end
    end
    D2=sum(D2Chunks)
    return D2
end

#FUNCTIONS EXPLOTING CENTER-OF-MASS MOMENTUM CONSERVATION FOR 2D SYSTEMS WITH DISCRETE TRANSLATION SYMMETRY.

#= 

Momentum space setup and definitions:

On a finite two-dimensional torus defined by the vectors L1, L2, all allowed momenta can be written as an integer linear
combination of T1, T2 where Ti = 2\pi \epsilon_ij (L_j \times \hat{z}) /|L_1 \times L_2|. Consider a crystal with primitive vectors
a1, a2 placed on this torus in a commensurate way (i.e. such that the torus fits an integer number of crystal unit cells.
The "mesh" is a set of N_uc = |L_1 \times L_2|/|a_1 \times a_2| momentum points that are unique modulo bi where 
bi = 2\pi \epsilon_ij (L_j \times \hat{z}) /|L_1 \times L_2| is a primitive reciporcal lattice vector. 
In this file, "M" is the change of basis matrix from the bi basis to the Ti basis and B = M^{-1} is the inverse change of basis matrix.

"RL" refers to the reciprocal lattice, which is all possible integer linear combinations of b1, b2.

The functions can handle multiband systems (i.e. systemw with more than one one-body state per spin per mesh point.)

Additional definitions beyond those introduced previously (not exhaustive):

nk: length of mesh = N_uc


=#


function altmod(a,b)::Int64
    if a%b == 0
        return b
    else
        return a%b
    end
end

function sendtomesh(B::Matrix{Int64}, v::Vector{Int64})::Vector{Int64} #https://crypto.stackexchange.com/questions/29661/how-to-find-the-value-of-a-vector-modulo-a-basis-in-lattice-based-cryptography/29684#29684
    #B takes a vector in b basis and changes it to T basis. So the way this works is that we take a vector v from T basis to B basis, get the integer part of it (ie the RL part of it), then convert that vector back to T basis and subtract it from the original vector in T basis
    #This gives us the non-RL part of the original vector v in the T basis, i.e. the mesh part.
    BInv=inv(B)
    return v .- B*Int.(floor.(round.(BInv*v, digits=5))) ##all the rounding is for numerical stability, ie so 1.999999999999 doesnt floor to 1.0
end

#calculate Fock states, organized into seperate center-of-mass crystal momentum blocks
#Setting the option "capexcite" to true discards Fock states with more than "maxexcite" excitations out of the lowest band.
#This function can apply to a single band by setting n=nk=length(mesh).
function calcFockStatesKBlocksMultiband(nU::Int64, nD::Int64, n::Int64, M::Matrix{Float64}, mesh::Vector{Vector{Int64}}, maxexcite::Int64=1, capexcite::Bool=false)::Vector{Vector{Int64}}
    nk = length(mesh)
    States = Vector{Vector{Int64}}()
    for i in 1:nk
        push!(States, [])
    end
    B=Int.(round.(inv(M)))
    newmesh=Vector{Vector{Int64}}()
    for i in eachindex(mesh)
        push!(newmesh, sendtomesh(B,mesh[i]))
    end
    indVec = collect(1:n) #reason why I dont just define n as nk is so that this function can handle multiband case
    occListsU = Combinatorics.combinations(indVec, nU)
    occListsD = Combinatorics.combinations(indVec, nD)
    for occListU in occListsU, occListD in occListsD
        if capexcite
            excitecount=sum(occListU .> nk) + sum(occListD .> nk)
            if excitecount > maxexcite
                continue
            end
        end 
        State = 0
        kVec = [0,0]
        for iU in occListU
            kVec+=mesh[altmod(iU, nk)]
            State += 2^iU
        end
        for iD in occListD
            kVec+=mesh[altmod(iD,nk)]
            State += 2^iD * 2^n
        end
        k=findfirst(isequal(sendtomesh(B, kVec)), newmesh)
        push!(States[k], State)
    end
    for i in 1:nk
        sort!(States[i])
    end
    return States
end

#calculate Fock states, organized into seperate center-of-mass crystal momentum blocks
function calcFockStatesKBlocks(nU::Int64, nD::Int64, n::Int64, M::Matrix{Float64}, mesh::Vector{Vector{Int64}})::Vector{Vector{Int64}}
    nk = length(mesh)
    States = Vector{Vector{Int64}}()
    for i in 1:nk
        push!(States, [])
    end
    B=Int.(inv(M))
    newmesh=Vector{Vector{Int64}}()
    for i in eachindex(mesh)
        push!(newmesh, sendtomesh(B,mesh[i]))
    end
    indVec = collect(1:n) #reason why I dont just define n as nk is so that this function can handle multiband case
    occListsU = Combinatorics.combinations(indVec, nU)
    occListsD = Combinatorics.combinations(indVec, nD)
    for occListU in occListsU, occListD in occListsD
        State = 0
        kVec = [0,0]
        for iU in occListU
            kVec+=mesh[iU]  
            State += 2^iU
        end
        for iD in occListD
            kVec+=mesh[iD]
            State += 2^iD * 2^n
        end
        k=findfirst(isequal(sendtomesh(B, kVec)), newmesh)
        push!(States[k], State)
    end
    for i in 1:nk
        sort!(States[i])
    end
    return States
end

#one body density matrix, with explicit momentum and spin conservation for efficiency
function calcD1MomSpinCos(State::Vector{ComplexF64},bandList::Vector{Int64}, sVals::Vector{Int64}, mesh::Vector{Vector{Int64}}, M::Matrix{Float64}, Basis::Vector{Int64})::Array{ComplexF64, 4}
    #D1[k,s,m,n]
    println("num threads calcD1: $(Threads.nthreads())")
    nk = length(mesh)
    numBands=length(bandList)
    n=nk*numBands
    B=Int.(inv(M))
    newmesh=Vector{Vector{Int64}}()
    for i in eachindex(mesh)
        push!(newmesh, sendtomesh(B,mesh[i]))
    end
    D1Chunks=[zeros(ComplexF64, (nk,2,numBands,numBands)) for _ in 1:Threads.nthreads()] #[k,s,m,o]
    sValsShift=sVals.-1
    Threads.@threads for IInd in eachindex(Basis)
        for k in 1:nk, s in sValsShift, m in 1:numBands, o in 1:numBands
            sign=1
            I=Basis[IInd]
            J,sign=C(k+(o-1)*nk,Bool(s),I,sign,n)
            if sign==0 continue end
            J,sign=CDag(k+(m-1)*nk,Bool(s),J,sign,n)
            if sign==0 continue end
            JInd = searchsortedfirst(Basis,J)
            if JInd == (length(Basis)+1) || Basis[JInd]!=J continue end
            D1Chunks[Threads.threadid()][k,s+1,m,o]+=conj(State[JInd])*State[IInd]*sign
        end
    end
    D1=sum(D1Chunks)
    return D1
end

#two body density matrix, with explicit momentum and spin conservation for efficiency
function calcD2MomSpinConsOneBand(State::Vector{ComplexF64}, sVals::Vector{Int64}, mesh::Vector{Vector{Int64}}, M::Matrix{Float64}, Basis::Vector{Int64})::Array{ComplexF64, 5}
    println("num threads calcD2: $(Threads.nthreads())")
    nk = length(mesh)
    B=Int.(inv(M))
    newmesh=Vector{Vector{Int64}}()
    for i in eachindex(mesh)
        push!(newmesh, sendtomesh(B,mesh[i]))
    end
    D2Chunks=[zeros(ComplexF64, (nk,2,nk,2,nk)) for _ in 1:Threads.nthreads()] #[p,sp,k,sk,q]
    sValsShift=sVals.-1
    Threads.@threads for IInd in eachindex(Basis)
        I=Basis[IInd]
        for k in 1:nk, sk in sValsShift
            sign=1
            J1,sign1=C(k,Bool(sk),I,sign,nk)
            if sign1==0 continue end
            for p in 1:nk, sp in sValsShift
                J2,sign2=C(p,Bool(sp),J1,sign1,nk)
                if sign2==0 continue end
                for q in 1:nk
                    pMinusqBZ=findfirst(isequal(sendtomesh(B, mesh[p]-mesh[q])), newmesh)
                    kPlusqBZ=findfirst(isequal(sendtomesh(B, mesh[k]+mesh[q])), newmesh)
                    J3,sign3=CDag(pMinusqBZ,Bool(sp),J2,sign2,nk)
                    if sign3==0 continue end
                    J4,sign4=CDag(kPlusqBZ,Bool(sk),J3,sign3,nk)
                    if sign4==0 continue end
                    JInd = searchsortedfirst(Basis,J4)
                    if JInd == (length(Basis)+1) || Basis[JInd]!=J4 continue end
                    D2Chunks[Threads.threadid()][p,sp+1,k,sk+1,q]+=conj(State[JInd])*State[IInd]*sign4
                end
            end
        end
    end
    D2=sum(D2Chunks)
    return D2
end

end