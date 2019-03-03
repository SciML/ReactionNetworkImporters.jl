using DiffEqBase, DiffEqBiological, OrdinaryDiffEq, Sundials
using ReactionNetworkImporters
using TimerOutputs

# parameters
networkname = "testbcrbng"
tf          = 10000.
datadir     = joinpath(@__DIR__,"../data/bcr")
fname       = joinpath(datadir, "bcr.net")

const to = TimerOutput()
reset_timer!(to)

# get the reaction network
@timeit to "netgen" prnbng = loadrxnetwork(BNGNetwork(),string(networkname,"bng"), fname); 
rn = prnbng.rn; u0 = prnbng.u₀; p = prnbng.p; 
@timeit to "addodes" addodes!(rn; build_jac=false, build_symfuncs=false)
@timeit to "ODEProb" oprob = ODEProblem(rn, u0, (0.,tf), p)
show(to)

# note solvers run faster the second time 

# Rodas4
#reset_timer!(to); @timeit to "Rodas4" begin sol = solve(oprob,Rodas4(autodiff=false),dense=false,calck=false); end; show(to)
#reset_timer!(to); @timeit to "Rodas4" begin sol = solve(oprob,Rodas4(autodiff=false),dense=false,calck=false); end; show(to)
#reset_timer!(to); @timeit to "Rodas4" begin sol = solve(oprob,Rodas4(autodiff=false),dense=false); end;

# similarly CVODE_BDF with gmres 
#reset_timer!(to); @timeit to "CVODE_BDF" begin sol = solve(oprob,CVODE_BDF(linear_solver=:GMRES),dense=false); end; show(to)
#reset_timer!(to); @timeit to "CVODE_BDF" begin sol = solve(oprob,CVODE_BDF(linear_solver=:GMRES),dense=false); end; show(to)

# CVODE_BDF with LU 
reset_timer!(to); @timeit to "CVODE_BDF-1" begin sol = solve(oprob,CVODE_BDF(),dense=false); end; show(to)
reset_timer!(to); @timeit to "CVODE_BDF-2" begin sol = solve(oprob,CVODE_BDF(),dense=false); end; show(to)

