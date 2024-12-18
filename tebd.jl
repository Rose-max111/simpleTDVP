using ITensors, ITensorMPS

function show_maxdim(psi::MPS)
    maxdim = 0
    for i in 1:length(psi)
        maxdim = max(maxdim, maximum(dim(inds(psi[i])[j]) for j in 1:length(inds(psi[i]))))
    end
    @show maxdim
end

function tebdmain(psi::MPS, tau::Float64, ttotal::Float64, cutoff::Float64)
    # Make an array of 'site' indices
    s = siteinds(psi)
    N = length(psi)

    # Make gates (1,2),(2,3),(3,4),...
    gates = ITensor[]
    for j in 1:(N - 1)
      s1 = s[j]
      s2 = s[j + 1]
      hj =
        op("Sz", s1) * op("Sz", s2) +
        1 / 2 * op("S+", s1) * op("S-", s2) +
        1 / 2 * op("S-", s1) * op("S+", s2)
      Gj = exp(-im * tau / 2 * hj)
      push!(gates, Gj)
    end
    # Include gates in reverse order too
    # (N,N-1),(N-1,N-2),...
    append!(gates, reverse(gates))

    c = div(N, 2) # center site

    # Compute and print <Sz> at each time step
    # then apply the gates to go to the next time
    @info "ITensor TEBD results"
    Sz_expect = Float64[]
    for t in 0.0:tau:ttotal
      Sz = expect(psi, "Sz"; sites=c)
      println("$t $Sz")
      push!(Sz_expect, Sz)
      
      # show_maxdim(psi)

      t≈ttotal && break

      psi = apply(gates, psi; cutoff)
      normalize!(psi)
    end
    @info "ITensor TEBD finished"
    return Sz_expect
    # time = collect(0.0:tau:ttotal)

  #   fig = Figure()
  #   ax = Axis(fig[1, 1])
  #   scatter!(ax, time, Sz_expect)
  #   fig
end

# main()