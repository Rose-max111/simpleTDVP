using CairoMakie

include("tdvp.jl")
include("tebd.jl")

function main(N::Int, dt::Float64, Tmax::Float64, cutoff::Float64)
    H, psi = hamiltonian(N)

    Sz_tebd = tebdmain(psi, dt, Tmax, cutoff)

    c=div(length(psi),2)
    Sz_tdvp = Float64[]
    for t in 0.0:dt:Tmax
        Sz = expect(psi, "Sz"; sites=c)
        println("$t $Sz")
        push!(Sz_tdvp, Sz)
        # show_maxdim(psi)
        t ≈ Tmax && break
        psi = tdvp2sweep(psi, H, dt, cutoff)
    end

    return Sz_tebd, Sz_tdvp, collect(0.0:dt:Tmax)
end

function plot_result(Sz_tebd, Sz_tdvp, time)
    fig = Figure(size=(1000, 500))
    ax1 = Axis(fig[1, 1], title = "Sz vs Time", xlabel = "Time", ylabel = "Sz")
    scatter!(ax1, time, Sz_tebd, label = "TEBD")
    lines!(ax1, time, Sz_tdvp, label = "TDVP", color = "red")

    ax2 = Axis(fig[1, 2], title = "Absolute Error", xlabel = "Time", ylabel = "Error", yticks = collect(0.0:1e-5:1e-4))
    scatter!(ax2, time, abs.(Sz_tebd .- Sz_tdvp), label = "Error")
    ylims!(ax2, 0.0, 1e-4)
    
    axislegend(ax1, position = :lt)
    axislegend(ax2, position = :lt)
    fig
    save("result.png", fig)
end

Sz_tebd, Sz_tdvp, time = main(40, 0.1, 5.0, 1E-8)
plot_result(Sz_tebd, Sz_tdvp, time)
