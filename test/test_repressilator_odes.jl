using DiffEqBase, DiffEqBiological, Plots, OrdinaryDiffEq, DataFrames, CSVFiles, LinearAlgebra
using ReactionNetworkImporters

# parameters
doplot = false
networkname = "Repressilator"
tf = 40000
nsteps = 400

# BNG simulation data
datadir  = joinpath(@__DIR__,"../data/repressilator")
fname    = joinpath(datadir, "Repressilator.net")
gdatfile = joinpath(datadir, "Repressilator.gdat")
print("getting gdat file...")
gdatdf = DataFrame(load(File(format"CSV", gdatfile), header_exists=true, spacedelim=true))
println("done")

# load the BNG reaction network in DiffEqBio
prnbng = loadrxnetwork(BNGNetwork(),string(networkname,"bng"), fname)
rnbng = prnbng.rn; u0 = prnbng.u₀; p = prnbng.p; 
addodes!(rnbng; build_jac=false, build_symfuncs=false)
boprob = ODEProblem(rnbng, u0, (0.,tf), p)

# BNG simulation data testing
pTetRid = prnbng.groupstoids[:pTetR][1]
bngsol = gdatdf[:pTetR]

# note solvers run _much_ faster the second time 
bsol = solve(boprob, Tsit5(), abstol=1e-12, reltol=1e-12, saveat=tf/nsteps); 

if doplot
    plotlyjs()
    p1 = plot(gdatdf[:time], gdatdf[:pTetR], label=:BNGPTetR)    
    plot!(p1, bsol.t, bsol[pTetRid,:], label=:DEBPTetR) 
    display(p1)    
end

@test all(bsol.t .== gdatdf[:time])
@test all( abs.(gdatdf[:pTetR] - bsol[pTetRid,:]) .< 1e-6 .* abs.(gdatdf[:pTetR]))

