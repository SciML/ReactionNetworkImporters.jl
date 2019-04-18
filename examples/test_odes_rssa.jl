# example using network from the RSSA Models
# download input files from: https://gist.github.com/isaacsas/7648fe065f1884307a2906a0a48ea12b

using DiffEqBase, DiffEqBiological, OrdinaryDiffEq, Sundials
using ReactionNetworkImporters
using TimerOutputs

# parameters
networkname = "tester"
tf = 10.

# input files, specify path 
datadir  = joinpath(@__DIR__,"../data/repressilator")
fname    = joinpath(datadir, "Repressilator.net")

const to = TimerOutput()
reset_timer!(to)

# get the reaction network
@timeit to "netgen" prn = loadrxnetwork(RSSAFile(), networkname, speciesf, rxsf; printrxs = false)
rn = prn.rn; initialpop = prn.u₀
@timeit to "addodes" addodes!(rn; build_jac=false, build_symfuncs=false, build_paramjac=false)
@timeit to "ODEProb" oprob = ODEProblem(rn,convert.(Float64,initialpop),(0.,tf))
show(to)
println()

# note solvers run faster the second time 

# I haven't been able to successfully solve the system with Rodas4/Rodas5
#reset_timer!(to); @timeit to "Rodas4" begin sol = solve(oprob,Rodas4(autodiff=false),dense=false,calck=false); end; show(to)
#reset_timer!(to); @timeit to "Rodas4" begin sol = solve(oprob,Rodas4(autodiff=false),dense=false,calck=false); end; show(to)
#reset_timer!(to); @timeit to "Rodas4" begin sol = solve(oprob,Rodas4(autodiff=false),dense=false); end;

# similarly CVODE_BDF with gmres is very slow and I haven't been able to complete a simulation
#reset_timer!(to); @timeit to "CVODE_BDF" begin sol = solve(oprob,CVODE_BDF(linear_solver=:GMRES),dense=false); end; show(to)
#reset_timer!(to); @timeit to "CVODE_BDF" begin sol = solve(oprob,CVODE_BDF(linear_solver=:GMRES),dense=false); end; show(to)

# CVODE_BDF with LU seems to finish, on my machine it takes upwards of 500 seconds
reset_timer!(to); @timeit to "CVODE_BDF" begin sol = solve(oprob,CVODE_BDF(),dense=false); end; show(to)
reset_timer!(to); @timeit to "CVODE_BDF" begin sol = solve(oprob,CVODE_BDF(),dense=false); end; show(to)

