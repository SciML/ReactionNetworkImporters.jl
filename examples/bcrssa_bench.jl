using ReactionNetworkImporters, Catalyst, ModelingToolkit, DiffEqJump
using Plots

datadir  = joinpath(@__DIR__,"../data/bcr_ssa")
fname    = joinpath(datadir, "BCRSSA.net")
tequil   = 30000.0

prn = loadrxnetwork(BNGNetwork(), fname)
rn = prn.rn
jsys = convert(JumpSystem, rn)
defs = ModelingToolkit.defaults(rn)
u0 = convert.(Int, ModelingToolkit.varmap_to_vars(ModelingToolkit.defaults(jsys), states(jsys)))
dprob = DiscreteProblem(jsys, u0, (0.0,tequil), [])
@assert eltype(dprob.u0) <: Int
jprob = JumpProblem(jsys, dprob, RSSACR(), save_positions=(false,false))

# try solving as a test
sol = solve(jprob, SSAStepper(), saveat=tequil/1000)

@unpack Activated_Syk = rn
plot(sol, vars=Activated_Syk)

u0 = sol[end]
tf = 20000.0
setdefaults!(rn, [:c => 3.0])
jsys = convert(JumpSystem, rn)
dprob = DiscreteProblem(jsys, u0, (0.0,tf), [])
@assert eltype(dprob.u0) <: Int
jprob = JumpProblem(jsys, dprob, RSSACR(), save_positions=(false,false))
sol = solve(jprob, SSAStepper(), saveat=tf/20000)
plot(sol, vars=Activated_Syk)