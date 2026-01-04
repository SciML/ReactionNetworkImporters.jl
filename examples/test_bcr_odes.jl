# load a BioNetGen defined BCR network, solve the ODEs
# plot a specific observable and compare to the BioNetGen solution

using DiffEqBase, Catalyst, Plots, OrdinaryDiffEq, Sundials, DataFrames, CSVFiles,
    LinearAlgebra, SparseArrays
using ReactionNetworkImporters
using TimerOutputs

# parameters
doplot = false
tf = 10000.0
build_jac = false
sparse_jac = false
zeroout_jac = false

# BNG simulation data
datadir = joinpath(@__DIR__, "../data/bcr")
fname = joinpath(datadir, "bcr.net")
gdatfile = joinpath(datadir, "bcr.gdat")
print("getting gdat file...")
gdatdf = DataFrame(
    load(
        File(format"CSV", gdatfile), header_exists = true,
        spacedelim = true
    )
)
println("done")

# we'll time the DiffEq solvers
const to = TimerOutput()
reset_timer!(to)

# BioNetGen network
@timeit to "bionetgen" prnbng = loadrxnetwork(BNGNetwork(), fname);
show(to)
rnbng = prnbng.rn
@timeit to "bODESystem" bosys = convert(ODESystem, rnbng)
show(to)
@timeit to "bODEProb" boprob = ODEProblem(bosys, Float64[], (0.0, tf), Float64[])
show(to)
u = zeros(numspecies(rnbng))
p = zeros(length(parameters(rnbng)))
du = similar(u);
@timeit to "f1" boprob.f(du, u, p, 0.0)
@timeit to "f2" boprob.f(du, u, p, 0.0)
if build_jac
    J = zeros(length(u), length(u))
    #J = similar(rnbng.odefun.jac_prototype)
    @timeit to "J1" rnbng.jac(J, u, p, 0.0)
    @timeit to "J2" rnbng.jac(J, u, p, 0.0)
end
show(to)
println()

# DiffEq solver
#reset_timer!(to);
#@timeit to "BNG_CVODE_BDF-LU-1" begin bsol = solve(boprob, CVODE_BDF(), saveat=1., abstol=1e-8, reltol=1e-8); end; show(to)
#@timeit to "BNG_CVODE_BDF-LU-2" begin bsol = solve(boprob, CVODE_BDF(), saveat=1., abstol=1e-8, reltol=1e-8); end; show(to)
#@timeit to "BNG_CVODE_BDF-GMRES-1" begin bsol = solve(boprob, CVODE_BDF(linear_solver=:GMRES), saveat=1., abstol=1e-8, reltol=1e-8); end; show(to)
#@timeit to "BNG_CVODE_BDF-GMRES-2" begin bsol = solve(boprob, CVODE_BDF(linear_solver=:GMRES), saveat=1., abstol=1e-8, reltol=1e-8); end; show(to)
# #reset_timer!(to); @timeit to "BNG_RODAS5_BDF" begin bsol2 = solve(boprob, rodas5(autodiff=false), saveat=1., abstol=1e-8, reltol=1e-8); end; show(to)
#@timeit to "KenCarp4-1" begin sol = solve(boprob,KenCarp4(autodiff=false,linsolve=LinSolveFactorize(lu)), saveat=1., abstol=1e-8, reltol=1e-8); end; show(to)
#@timeit to "KenCarp4-2" begin bsol = solve(boprob,KenCarp4(autodiff=false,linsolve=LinSolveFactorize(lu)), saveat=1., abstol=1e-8, reltol=1e-8); end; show(to)
@timeit to "KenCarp4-1" begin
    sol = solve(
        boprob, KenCarp4(autodiff = false), abstol = 1.0e-8,
        reltol = 1.0e-8, saveat = 1.0
    )
end;
show(to);
@timeit to "KenCarp4-2" begin
    bsol = solve(
        boprob, KenCarp4(autodiff = false),
        abstol = 1.0e-8, reltol = 1.0e-8, saveat = 1.0
    )
end;
show(to);

# Activated Syk from DiffEq
@unpack Activated_Syk = rnbng

if doplot
    plotlyjs()
    plot(
        gdatdf[!, :time][2:end], gdatdf[!, :Activated_Syk][2:end], xscale = :log10,
        label = "AsykGroup", linestyle = :dot
    )
    plot!(
        bsol.t[2:end], sol[Activated_Syk][2:end], label = "Activated_Syk",
        xscale = :log10
    )
end

# test the error, note may be large in abs value though small relatively
# #norm(gdatdf[:Activated_Syk] - asynbng, Inf)
# #norm(asynbng - basyk', Inf)
norm(gdatdf[!, :Activated_Syk] - sol[Activated_Syk], Inf)

#@assert all(abs.(gdatdf[!,:Activated_Syk] - basyk') .< 1e-6 * abs.(gdatdf[:Activated_Syk]))
