# download input files from: https://gist.github.com/isaacsas/7648fe065f1884307a2906a0a48ea12b
# then add ReactionNetworkImporters package (unregistered)
#] 
# add https://github.com/isaacsas/ReactionNetworkImporters.jl.git

using DiffEqBase, DiffEqBiological, Plots, OrdinaryDiffEq, Sundials, DataFrames, CSVFiles, LinearAlgebra
using ReactionNetworkImporters
using TimerOutputs

# parameters
doplot = true
networkname = "testbcrbng"
tf = 10000.

# BNG simulation data
datadir  = joinpath(@__DIR__,"../data/bcr")
fname    = joinpath(datadir, "bcr.net")
cdatfile = joinpath(datadir, "bcr.cdat")
gdatfile = joinpath(datadir, "bcr.gdat")
print("getting cdat file...")
#cdatdf = CSV.File(cdatfile, delim=" ", ignorerepeated=true) |> DataFrame
cdatdf = DataFrame(load(File(format"CSV", cdatfile), header_exists=true, spacedelim=true) )
println("done")
print("getting gdat file...")
#gdatdf = CSV.File(gdatfile, delim=" ", ignorerepeated=true) |> DataFrame
gdatdf = DataFrame(load(File(format"CSV", gdatfile), header_exists=true, spacedelim=true) )
println("done")

const to = TimerOutput()
reset_timer!(to)

@timeit to "bionetgen" prnbng = loadrxnetwork(BNGNetwork(),string(networkname,"bng"), fname); 
rnbng = prnbng.rn; u₀ = prnbng.u₀; p = prnbng.p; 

u0 = u₀   #convert.(Float64, prn.u₀)

# get the bng reaction network
@timeit to "baddodes" addodes!(rnbng; build_jac=false, build_symfuncs=false)
@timeit to "bODEProb" boprob = ODEProblem(rnbng, u0, (0.,tf), p)
show(to)



# BNG simulation data testing
asykgroups = prnbng.groupstoids[:Activated_Syk]
asyksyms = findall(x -> x ∈ asykgroups, rnbng.syms_to_ints)
asynbng = zeros(length(cdatdf[:time]))
for sym in asyksyms
    global asynbng
    asynbng += cdatdf[sym]
end


# note solvers run _much_ faster the second time 
reset_timer!(to); @timeit to "BNG_CVODE_BDF" begin bsol = solve(boprob, CVODE_BDF(),dense=false, saveat=1., abstol=1e-8, reltol=1e-8); end; show(to)

basyk = sum(bsol[asykgroups,:], dims=1)
#vars = sum(bsol[prnbng.groupstoids[:Ig_alpha_P],:],dims=1)

if doplot
    plotlyjs()
    plot(gdatdf[:time][2:end], gdatdf[:Activated_Syk][2:end], xscale=:log10, label=:AsykGroup)
    plot!(cdatdf[:time][2:end], asynbng[2:end], xscale=:log10, label=:AsykSum)
    plot!(bsol.t[2:end], basyk'[2:end], label=:AsykDEBio, xscale=:log10)
    # plot!(bsol.t[2:end], vars'[2:end], label=:Ig_alpha_P, xscale=:log10)
end

# test the error
norm(gdatdf[:Activated_Syk] - asynbng, Inf)
norm(gdatdf[:Activated_Syk] - basyk', Inf)
norm(asynbng - basyk', Inf)