# example using network from the RSSA Models
# download input files from: https://gist.github.com/isaacsas/7648fe065f1884307a2906a0a48ea12b

using DiffEqBase, DiffEqBiological, OrdinaryDiffEq, Sundials
using ReactionNetworkImporters
using TimerOutputs

# parameters
networkname = "tester"
tf = 10.

# input files, specify path 
#speciesf = "SOME-PATH-TO/BCR_pop.txt"
#rxsf = "SOME-PATH-TO/BCR_rxn.txt"

const to = TimerOutput()
reset_timer!(to)

# get the reaction network
@timeit to "netgen" prn = loadrxnetwork(RSSAFile(), networkname, speciesf, rxsf; printrxs = false)
rn = prn.rn; initialpop = prn.u₀
@timeit to "addodes" addodes!(rn; build_jac=false, build_symfuncs=false)
@timeit to "ODEProb" oprob = ODEProblem(rn,convert.(Float64,initialpop),(0.,tf))
show(to)

# note solvers run faster the second time 
reset_timer!(to); @timeit to "Rodas4" begin sol = solve(oprob,Rodas4(autodiff=false),dense=false,calck=false); end; show(to)
reset_timer!(to); @timeit to "Rodas4" begin sol = solve(oprob,Rodas4(autodiff=false),dense=false,calck=false); end; show(to)
reset_timer!(to); @timeit to "Rodas4" begin sol = solve(oprob,Rodas4(autodiff=false),dense=false); end;

reset_timer!(to); @timeit to "CVODE_BDF" begin sol = solve(oprob,CVODE_BDF(linear_solver=:GMRES),dense=false); end; show(to)
reset_timer!(to); @timeit to "CVODE_BDF" begin sol = solve(oprob,CVODE_BDF(linear_solver=:GMRES),dense=false); end; show(to)

# note, restart REPL to see timing with first pass overhead for the next solver
reset_timer!(to); @timeit to "CVODE_BDF" begin sol = solve(oprob,CVODE_BDF(),dense=false); end; show(to)
reset_timer!(to); @timeit to "CVODE_BDF" begin sol = solve(oprob,CVODE_BDF(),dense=false); end; show(to)

