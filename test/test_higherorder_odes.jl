using DiffEqBase, Catalyst, Plots, OrdinaryDiffEq, DataFrames, CSVFiles, LinearAlgebra
using ReactionNetworkImporters

# parameters
doplot = false
networkname = "HigherOrder"
tf = 1e-2
nsteps = 2000

# BNG simulation data
datadir  = joinpath(@__DIR__,"../data/higherorder")
fname    = joinpath(datadir, "higherorder.net")
gdatfile = joinpath(datadir, "higherorder.gdat")
print("getting gdat file...")
gdatdf = DataFrame(load(File{format"CSV"}(gdatfile), header_exists=true, spacedelim=true))
println("done")

# load the BNG reaction network in DiffEqBio
prnbng = loadrxnetwork(BNGNetwork(), fname)
rnbng = prnbng.rn; u0 = prnbng.u₀; p = prnbng.p; 
boprob = ODEProblem(rnbng, u0, (0.,tf), p)

# BNG simulation data testing
@unpack A = rnbng
Asol = gdatdf[!,:A]

# note solvers run _much_ faster the second time 
bsol = solve(boprob, Tsit5(), abstol=1e-12, reltol=1e-12, saveat=tf/nsteps); 

if doplot
    plotlyjs()
    p1 = plot(gdatdf[!,:time], gdatdf[!,:A], label=:A_BNG)    
    plot!(p1, bsol.t, bsol[A], label=:A_DE)
    display(p1)
    println("Err = ", norm(gdatdf[!,:A] - bsol[Aid,:],Inf))
end

@test all(bsol.t .== gdatdf[!,:time])
@test all(abs.(gdatdf[!,:A] - bsol[A]) .< 1e-6 .* abs.(gdatdf[!,:A]))

