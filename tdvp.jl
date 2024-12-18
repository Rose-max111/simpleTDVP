using ITensors, ITensorMPS
using ITensorTDVP
using KrylovKit

function hamiltonian(N::Int)
    ```
    Construct a 1D open boundary Heisenberg Hamiltonian with N sites, along with an initial MPS
    ```
    sites = siteinds("S=1/2", N)
    os = OpSum()
    for j=1:N-1
        os += "Sz",j,"Sz",j+1
        os += 1/2,"S+",j,"S-",j+1
        os += 1/2,"S-",j,"S+",j+1
    end
    H = MPO(os,sites)
    psi = MPS(s, n -> isodd(n) ? "Up" : "Dn")
    return H, psi
end

function applyH2(AC, lenv, M1, M2, renv)
    ```
    Apply the two-site operator M1*M2 along with left and right environment tensors to the center tensor AC
    
    AC: center tensor
    lenv: left environment tensor
    M1: left hamiltonian operator
    M2: right hamiltonian operator
    renv: right environment tensor
    ```
    return (((lenv * AC) * M1) * M2) * renv
end

function applyH1(AC, lenv, M1, renv)
    ```
    Apply the one-site operator M1 along with left and right environment tensors to the center tensor AC
    
    AC: center tensor
    lenv: left environment tensor
    M1: hamiltonian operator
    renv: right environment tensor
    ```
    return ((lenv * AC) * M1) * renv
end

function updateleftenv(FL, A, M)
    ```
    Update left environment tensor

    FL: previous left environment tensor
    A: Projected MPS tensor (Left orthogonalzed)
    M: Hamiltonian operator
    ```
    return ((FL*A)*M)*dag(prime(A))
end

function updaterightenv(FR, A, M)
    ```
    Update right environment tensor

    FR: previous right environment tensor
    A: Projected MPS tensor (Right orthogonalzed)
    M: Hamiltonian operator
    ```
    return ((FR*A)*M)*dag(prime(A))
end

function show_maxdim(psi::MPS)
    ```
    Show the maximum bound dimension of the MPS
    ```
    maxdim = 0
    for i in 1:length(psi)
        maxdim = max(maxdim, maximum(dim(inds(psi[i])[j]) for j in 1:length(inds(psi[i]))))
    end
    @show maxdim
end

function tdvp2sweep(psi::MPS, H::MPO, Δt::Float64, cutoff::Float64)
    ```
    Perform a 2site-TDVP sweep from left to right, then right to left. Return a final MPS |ψ(t+Δt)⟩

    psi: initial MPS |ψ(t)⟩
    H: Hamiltonian MPO
    Δt: time step
    cutoff: singular value cutoff
    ```
    psi = orthogonalize(psi, 1) # make initial MPS be right orthogonal
    n = length(psi)

    Env = Vector{ITensor}(undef, n+2) # environment tensors
    Env[n+2] = ITensor(1.)
    Env[1] = ITensor(1.)
    for i in n+1:-1:3
        Env[i] = updaterightenv(Env[i+1], psi[i-1], H[i-1]) # precompute right environment tensors
    end

    for i in 1:n-1 # sweep from left to right, each evolve Δt/2
        wf, info = exponentiate(x->noprime(applyH2(x, Env[i], H[i], H[i+1], Env[i+3])), -im*Δt/2, psi[i]*psi[i+1]) # evolve two-site MPS

        indsi = uniqueinds(psi[i], psi[i+1]) # select virtual indice (i-1, i) and physical indice (i) to be row indices
        U, S, V = svd(wf, indsi, cutoff=cutoff)

        psi[i] = U
        psi[i+1] = S * V
        
        Env[i+1] = updateleftenv(Env[i], psi[i], H[i]) # update left environment tensor
        if i!=n-1
            wf, info = exponentiate(x->noprime(applyH1(x, Env[i+1], H[i+1], Env[i+3])), +im*Δt/2, psi[i+1]) # evolve one-site MPS
            psi[i+1] = wf
        end
    end

    for i in n:-1:2 # sweep from right to left, each evolve Δt/2
        wf, info = exponentiate(x->noprime(applyH2(x, Env[i-1], H[i-1], H[i], Env[i+2])), -im*Δt/2, psi[i-1]*psi[i]) # evolve two-site MPS
        
        indsi = uniqueinds(psi[i-1], psi[i]) # select virtual indice (i, i+1) and physical indice (i) to be row indices
        U, S, V = svd(wf, indsi, cutoff=cutoff)

        psi[i-1] = U * S
        psi[i] = V

        Env[i+1] = updaterightenv(Env[i+2], psi[i], H[i]) # update right environment tensor
        if i!=2
            wf, info = exponentiate(x->noprime(applyH1(x, Env[i-1], H[i-1], Env[i+1])), +im*Δt/2, psi[i-1]) #evolve one-site MPS
            psi[i-1] = wf
        end
    end
    normalize!(psi)
    return psi
end


function ITensortdvpmain(psi::MPS, H::MPO, dt, Tmax)
    @info "ITensor TDVP results"
    c=div(length(psi),2)
    for t in 0.0:dt:Tmax
        Sz = expect(psi, "Sz"; sites=c)
        println("$t $Sz")
        # show_maxdim(psi)

        t ≈ Tmax && break
        psi = tdvp(-im*H, dt, psi; time_step = dt)
    end
    @info "ITensor TDVP results finished"
end

