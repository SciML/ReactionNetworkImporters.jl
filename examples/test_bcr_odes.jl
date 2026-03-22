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
@timeit to "bionetgen" rnbng = loadrxnetwork(BNGNetwork(), fname);
show(to)
rnbng = complete(rnbng)
@timeit to "bODESystem" bosys = ode_model(rnbng)
show(to)
@timeit to "bODEProb" boprob = ODEProblem(bosys, Float64[], (0.0, tf), Float64[])
show(to)
u = zeros(numspecies(rnbng))
p = zeros(length(parameters(rnbng)))
du = similar(u);
@timeit to "f1" boprob.f(du, u, p, 0.0)
@timeit to "f2" boprob.f(du, u, p, 0.0)
show(to)
println()

# DiffEq solver
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
norm(gdatdf[!, :Activated_Syk] - sol[Activated_Syk], Inf)
