# download input files from: https://gist.github.com/isaacsas/7648fe065f1884307a2906a0a48ea12b
# then add ReactionNetworkImporters package (unregistered)
#] 
# add https://github.com/isaacsas/ReactionNetworkImporters.jl.git

using DiffEqBase, DiffEqBiological, Plots, OrdinaryDiffEq, DataFrames, CSVFiles, LinearAlgebra
using ReactionNetworkImporters

# parameters
doplot = false
networkname = "testbcrbng"
tf = 40000
steps = 400

# BNG simulation data
datadir  = joinpath(@__DIR__,"../data/repressilator")
fname    = joinpath(datadir, "Repressilator.net")
gdatfile = joinpath(datadir, "Repressilator.gdat")
print("getting gdat file...")
gdatdf = DataFrame(load(File(format"CSV", gdatfile), header_exists=true, spacedelim=true))
println("done")

# load the BNG reaction network in DiffEqBio
prnbng = loadrxnetwork(BNGNetwork(),string(networkname,"bng"), fname)
rnbng = prnbng.rn; u₀ = prnbng.u₀; p = prnbng.p; 
u0 = u₀   
addodes!(rnbng; build_jac=false, build_symfuncs=false)
boprob = ODEProblem(rnbng, u0, (0.,tf), p)

# BNG simulation data testing
pTetRid = prnbng.groupstoids[:pTetR][1]
bngsol = gdatdf[:pTetR]

# note solvers run _much_ faster the second time 
bsol = solve(boprob, Tsit5(), abstol=1e-12, reltol=1e-12, saveat=tf/(4*10^2)); 

if doplot
    plotlyjs()
    plot(gdatdf[:time], gdatdf[:pTetR], label=:BNGPTetR)    
    plot!(bsol.t, bsol[pTetRid,:], label=:DEBPTetR)
end

@test all(bsol.t .== gdatdf[:time])
@test all(abs.(gdatdf[:pTetR] - bsol[pTetRid,:]) .< 1e-6*abs.(gdatdf[:pTetR]))