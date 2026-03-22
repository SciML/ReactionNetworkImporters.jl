using ReactionNetworkImporters, Catalyst, JumpProcesses
using Plots

datadir = joinpath(@__DIR__, "../data/bcr_ssa")
fname = joinpath(datadir, "BCRSSA.net")

rn = loadrxnetwork(BNGNetwork(), fname)
rn = complete(rn)

# post equilibration solve
tf = 20000.0
u0 = [s => Int(v) for (s, v) in pairs(Catalyst.get_u0_map(rn))]
jprob = JumpProblem(rn, u0, (0.0, tf), []; save_positions = (false, false))
sol = solve(jprob, SSAStepper(), saveat = tf / 20000)
@unpack Activated_Syk = rn
plot(sol, idxs = Activated_Syk)
