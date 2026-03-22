using Catalyst, OrdinaryDiffEqTsit5, DataFrames, CSV, LinearAlgebra
using ReactionNetworkImporters

# parameters
doplot = false
networkname = "NullRxs"
tf = 50
nsteps = 500

# BNG simulation data
datadir = joinpath(@__DIR__, "../data/nullrxs")
fname = joinpath(datadir, "birth-death.net")
gdatfile = joinpath(datadir, "birth-death.gdat")
print("getting gdat file...")
gdatdf = CSV.read(gdatfile, DataFrame; delim = ' ', ignorerepeated = true)
println("done")

# load the BNG reaction network
rnbng = loadrxnetwork(BNGNetwork(), fname)

# test metadata is set
@test Catalyst.has_u0_map(rnbng)
@test Catalyst.has_parameter_map(rnbng)
@test has_varstonames(rnbng)
@test has_groupstosyms(rnbng)

rn = complete(rnbng)
u0 = Float64[]
p = Float64[]
boprob = ODEProblem(rn, u0, (0.0, tf), p)

# note solvers run _much_ faster the second time
bsol = solve(boprob, Tsit5(), abstol = 1.0e-12, reltol = 1.0e-12, saveat = tf / nsteps)

if doplot
    plotlyjs()
    p1 = plot(gdatdf[!, :time], gdatdf[!, :Atot], label = :A_BNG)
    plot!(p1, bsol.t, bsol[1, :], label = :A_DE)
    display(p1)
    println("Err = ", norm(gdatdf[!, :Atot] - bsol[1, :], Inf))
end

@test all(bsol.t .== gdatdf[!, :time])
@test all(
    abs.(gdatdf[!, :Atot] - bsol[1, :]) .<
        max.(1.0e-8 .* abs.(gdatdf[!, :Atot]), 1.0e-12)
)
