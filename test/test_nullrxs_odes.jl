using DiffEqBase, DiffEqBiological, Plots, OrdinaryDiffEq, DataFrames, CSVFiles, LinearAlgebra
using ReactionNetworkImporters

# parameters
doplot = false
networkname = "NullRxs"
tf = 50
nsteps = 500

# BNG simulation data
datadir  = joinpath(@__DIR__,"../data/nullrxs")
fname    = joinpath(datadir, "birth-death.net")
gdatfile = joinpath(datadir, "birth-death.gdat")
print("getting gdat file...")
gdatdf = DataFrame(load(File(format"CSV", gdatfile), header_exists=true, spacedelim=true))
println("done")

# load the BNG reaction network in DiffEqBio
prnbng = loadrxnetwork(BNGNetwork(),string(networkname,"bng"), fname)
rnbng = prnbng.rn; u0 = prnbng.u₀; p = prnbng.p; 
addodes!(rnbng)
boprob = ODEProblem(rnbng, u0, (0.,tf), p)

# note solvers run _much_ faster the second time 
bsol = solve(boprob, Tsit5(), abstol=1e-12, reltol=1e-12, saveat=tf/nsteps); 

if doplot
    plotlyjs()
    p1 = plot(gdatdf[:time], gdatdf[:Atot], label=:A_BNG)    
    plot!(p1, bsol.t, bsol[1,:], label=:A_DE)
    display(p1)
    println("Err = ", norm(gdatdf[:Atot] - bsol[1,:],Inf))
end

@test all(bsol.t .== gdatdf[:time])
@test all(abs.(gdatdf[:Atot] - bsol[1,:]) .< max.(1e-8 .* abs.(gdatdf[:Atot]), 1e-12))

