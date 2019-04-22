# NEEDS DiffEqBiological master as of 4/20/19
using DiffEqBase, DiffEqBiological, OrdinaryDiffEq, Sundials, LinearAlgebra, SparseArrays
using ReactionNetworkImporters
using TimerOutputs

# parameters
const networkname = "testbcrbng"
const tf = 10000.

# finite difference
const build_jac = false
const sparse_jac = false
const zeroout_jac = false

# dense Jacobian
# const build_jac = true
# const sparse_jac = false
# const zeroout_jac = true

# sparse Jacobian
# const build_jac = true
# const sparse_jac = true
# const zeroout_jac = true


# BNG simulation data
datadir  = joinpath(@__DIR__,"../data/bcr")
fname    = joinpath(datadir, "bcr.net")

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

function solve_dense(boprob, u, p, build_jac)
    if build_jac
        J = zeros(length(u),length(u))
        @timeit to "J1" rnbng.jac(J,u,p,0.)
        @timeit to "J2" rnbng.jac(J,u,p,0.)
    end
    show(to)
    println()
    
    @timeit to "KenCarp4-LU-1" begin bsol = solve(boprob,KenCarp4(autodiff=false), abstol=1e-8, reltol=1e-8, saveat=1.); end; show(to)
    @timeit to "KenCarp4-LU-2" begin bsol = solve(boprob,KenCarp4(autodiff=false), abstol=1e-8, reltol=1e-8, saveat=1.); end; show(to)    
    
    bsol 
end

function solve_sparse(boprob, u, p)
    J = similar(rnbng.odefun.jac_prototype)
    @timeit to "J1" rnbng.jac(J,u,p,0.)
    @timeit to "J2" rnbng.jac(J,u,p,0.)
    show(to)
    println()
    
    @timeit to "KenCarp4-SLU-1" begin bsol = solve(boprob,KenCarp4(autodiff=false,linsolve=LinSolveFactorize(lu)), saveat=1., abstol=1e-8, reltol=1e-8); end; show(to)
    @timeit to "KenCarp4-SLU-2" begin bsol = solve(boprob,KenCarp4(autodiff=false,linsolve=LinSolveFactorize(lu)), saveat=1., abstol=1e-8, reltol=1e-8); end; show(to)

    bsol
end

bsol = sparse_jac ? solve_sparse(boprob, u, p) : solve_dense(boprob, u, p, build_jac);

println()
