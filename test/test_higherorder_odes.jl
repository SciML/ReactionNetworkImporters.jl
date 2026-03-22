using Catalyst, OrdinaryDiffEqTsit5, DataFrames, CSV, LinearAlgebra
using ReactionNetworkImporters

# parameters
doplot = false
networkname = "HigherOrder"
tf = 1.0e-2
nsteps = 2000

# BNG simulation data
datadir = joinpath(@__DIR__, "../data/higherorder")
fname = joinpath(datadir, "higherorder.net")
gdatfile = joinpath(datadir, "higherorder.gdat")
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
boprob = ODEProblem(rn, Float64[], (0.0, tf), Float64[])

@show observed(rn)

# BNG simulation data testing
@unpack A = rn
Asol = gdatdf[!, :A]

# note solvers run _much_ faster the second time
bsol = solve(boprob, Tsit5(), abstol = 1.0e-12, reltol = 1.0e-12, saveat = tf / nsteps);

# if doplot
#     plotlyjs()
#     p1 = plot(gdatdf[!, :time], gdatdf[!, :A], label = :A_BNG)
#     plot!(p1, bsol.t, bsol[A], label = :A_DE)
#     display(p1)
#     println("Err = ", norm(gdatdf[!, :A] - bsol[Aid, :], Inf))
# end

@test all(bsol.t .== gdatdf[!, :time])
@test all(abs.(gdatdf[!, :A] - bsol[A]) .< 1.0e-6 .* abs.(gdatdf[!, :A]))
