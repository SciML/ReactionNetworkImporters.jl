using Catalyst, OrdinaryDiffEqTsit5, DataFrames, CSV, LinearAlgebra
using ReactionNetworkImporters

# parameters
doplot = false
networkname = "Repressilator"
tf = 40000
nsteps = 400

# BNG simulation data
datadir = joinpath(@__DIR__, "../data/repressilator")
fname = joinpath(datadir, "Repressilator.net")
gdatfile = joinpath(datadir, "Repressilator.gdat")
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

# BNG simulation data testing
bngsol = gdatdf[!, :pTetR]

# note solvers run _much_ faster the second time
bsol = solve(boprob, Tsit5(), abstol = 1.0e-12, reltol = 1.0e-12, saveat = tf / nsteps);
@unpack pTetR = rn

# if doplot
#     plotlyjs()
#     p1 = plot(gdatdf[!, :time], gdatdf[!, :pTetR], label = :BNGPTetR)
#     plot!(p1, bsol.t, bsol[pTetR], label = :DEBPTetR)
#     display(p1)
# end

@test all(bsol.t .== gdatdf[!, :time])
@test all(abs.(gdatdf[!, :pTetR] - bsol[pTetR]) .< 1.0e-6 .* abs.(gdatdf[!, :pTetR]))
