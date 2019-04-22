# load a BioNetGen defined BCR network, solve the ODEs
# plot a specific observable and compare to the BioNetGen solution

using DiffEqBase, DiffEqBiological, Plots, OrdinaryDiffEq, Sundials, DataFrames, CSVFiles, LinearAlgebra, SparseArrays
using ReactionNetworkImporters
using TimerOutputs

# parameters
doplot = true
networkname = "testbcrbng"
tf = 10000.
build_jac = true
sparse_jac = false
zeroout_jac = true

# BNG simulation data
datadir  = joinpath(@__DIR__,"../data/bcr")
fname    = joinpath(datadir, "bcr.net")
gdatfile = joinpath(datadir, "bcr.gdat")
print("getting gdat file...")
gdatdf = DataFrame(load(File(format"CSV", gdatfile), header_exists=true, spacedelim=true) )
println("done")

# we'll time the DiffEq solvers 
const to = TimerOutput()
reset_timer!(to)

# BioNetGen network
@timeit to "bionetgen" prnbng = loadrxnetwork(BNGNetwork(),string(networkname,"bng"), fname); 
rnbng = prnbng.rn; u0 = prnbng.u₀; p = prnbng.p; 
@timeit to "baddodes" addodes!(rnbng; build_jac=build_jac, zeroout_jac=zeroout_jac, sparse_jac=sparse_jac, build_symfuncs=false, build_paramjac=false)
@timeit to "bODEProb" boprob = ODEProblem(rnbng, u0, (0.,tf), p)
u = copy(u0);
du = similar(u);
@timeit to "f1" rnbng.f(du,u,p,0.)
@timeit to "f2" rnbng.f(du,u,p,0.)
if build_jac
    J = zeros(length(u),length(u))
    #J = similar(rnbng.odefun.jac_prototype)
    @timeit to "J1" rnbng.jac(J,u,p,0.)
    @timeit to "J2" rnbng.jac(J,u,p,0.)
end
show(to)
println()

# BNG simulation results for Activated Syk
asykgroups = prnbng.groupstoids[:Activated_Syk]
asyksyms = findall(x -> x ∈ asykgroups, rnbng.syms_to_ints)
# asynbng = zeros(length(gdatdf[:time]))
# for sym in asyksyms
#     global asynbng
#     asynbng += gdatdf[sym]
# end

# DiffEq solver 
#reset_timer!(to); 
#@timeit to "BNG_CVODE_BDF-LU-1" begin bsol = solve(boprob, CVODE_BDF(), saveat=1., abstol=1e-8, reltol=1e-8); end; show(to)
#@timeit to "BNG_CVODE_BDF-LU-2" begin bsol = solve(boprob, CVODE_BDF(), saveat=1., abstol=1e-8, reltol=1e-8); end; show(to)
#@timeit to "BNG_CVODE_BDF-GMRES-1" begin bsol = solve(boprob, CVODE_BDF(linear_solver=:GMRES), saveat=1., abstol=1e-8, reltol=1e-8); end; show(to)
#@timeit to "BNG_CVODE_BDF-GMRES-2" begin bsol = solve(boprob, CVODE_BDF(linear_solver=:GMRES), saveat=1., abstol=1e-8, reltol=1e-8); end; show(to)
# #reset_timer!(to); @timeit to "BNG_RODAS5_BDF" begin bsol2 = solve(boprob, rodas5(autodiff=false), saveat=1., abstol=1e-8, reltol=1e-8); end; show(to)
#@timeit to "KenCarp4-1" begin sol = solve(boprob,KenCarp4(autodiff=false,linsolve=LinSolveFactorize(lu)), saveat=1., abstol=1e-8, reltol=1e-8); end; show(to)
#@timeit to "KenCarp4-2" begin bsol = solve(boprob,KenCarp4(autodiff=false,linsolve=LinSolveFactorize(lu)), saveat=1., abstol=1e-8, reltol=1e-8); end; show(to)
@timeit to "KenCarp4-1" begin sol = solve(boprob,KenCarp4(autodiff=false), abstol=1e-8, reltol=1e-8, saveat=1.); end; show(to)
@timeit to "KenCarp4-2" begin bsol = solve(boprob,KenCarp4(autodiff=false), abstol=1e-8, reltol=1e-8, saveat=1.); end; show(to)

# Activated Syk from DiffEq
basyk = sum(bsol[asykgroups,:], dims=1)

if doplot
    plotlyjs()
    plot(gdatdf[:time][2:end], gdatdf[:Activated_Syk][2:end], xscale=:log10, label=:AsykGroup, linestyle=:dot)
#     # plot!(cdatdf[:time][2:end], asynbng[2:end], xscale=:log10, label=:AsykSum)
    plot!(bsol.t[2:end], basyk'[2:end], label=:AsykDEBio, xscale=:log10)
end

# test the error, note may be large in abs value though small relatively
# #norm(gdatdf[:Activated_Syk] - asynbng, Inf)
# #norm(asynbng - basyk', Inf)
norm(gdatdf[:Activated_Syk] - basyk', Inf)

# #@assert all(abs.(gdatdf[:Activated_Syk] - asynbng) .< 1e-6 * abs.(gdatdf[:Activated_Syk]))