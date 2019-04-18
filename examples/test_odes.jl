using DiffEqBase, DiffEqBiological, OrdinaryDiffEq, Sundials
using ReactionNetworkImporters
using TimerOutputs

# parameters
networkname = "testbcrbng"
tf          = 10000.
datadir     = joinpath(@__DIR__,"../data/bcr")
fname       = joinpath(datadir, "bcr.net")
build_jac   = true
sparse_jac  = false

const to = TimerOutput()
reset_timer!(to)

# get the reaction network
@timeit to "netgen" prnbng = loadrxnetwork(BNGNetwork(),string(networkname,"bng"), fname); 
rn = prnbng.rn; u0 = prnbng.u₀; p = prnbng.p; 
@timeit to "addodes" addodes!(rn; build_jac=build_jac, sparse_jac=sparse_jac, build_symfuncs=false, build_paramjac=false)
@timeit to "ODEProb" oprob = ODEProblem(rn, u0, (0.,tf), p)
u = copy(u0);
du = similar(u);
@timeit to "f1" rn.f(du,u,p,0.)
@timeit to "f2" rn.f(du,u,p,0.)
if build_jac
    J = zeros(length(u),length(u))
    @timeit to "J1" rn.jac(J,u,p,0.)
    @timeit to "J2" rn.jac(J,u,p,0.)
end
show(to)
println()

# note solvers run faster the second time 

# Rodas4
#reset_timer!(to); @timeit to "Rodas4" begin sol = solve(oprob,Rodas4(autodiff=false),dense=false,calck=false); end; show(to)
#reset_timer!(to); @timeit to "Rodas4" begin sol = solve(oprob,Rodas4(autodiff=false),dense=false,calck=false); end; show(to)
#reset_timer!(to); @timeit to "Rodas4" begin sol = solve(oprob,Rodas4(autodiff=false),dense=false); end;

# similarly CVODE_BDF with gmres 
dense = true
calck = true
# @timeit to "CVODE_BDF-GMRES-1" begin sol = solve(oprob,CVODE_BDF(linear_solver=:GMRES),dense=dense, abstol=1e-8, reltol=1e-8); end; 
# @timeit to "CVODE_BDF-GMRES-2" begin sol = solve(oprob,CVODE_BDF(linear_solver=:GMRES),dense=dense, abstol=1e-8, reltol=1e-8); end; 

# CVODE_BDF with LU 
@timeit to "CVODE_BDF-1" begin sol = solve(oprob,CVODE_BDF(),dense=dense, abstol=1e-8, reltol=1e-8); end; show(to)
@timeit to "CVODE_BDF-2" begin sol = solve(oprob,CVODE_BDF(),dense=dense, abstol=1e-8, reltol=1e-8); end; show(to)
@timeit to "KenCarp4-1" begin sol = solve(oprob,KenCarp4(autodiff=false),dense=dense, calck=calck, abstol=1e-8, reltol=1e-8); end; show(to)
@timeit to "KenCarp4-2" begin sol = solve(oprob,KenCarp4(autodiff=false),dense=dense, calck=calck, abstol=1e-8, reltol=1e-8); end;

show(to)