using ReactionNetworkImporters
using DiffEqBase, DiffEqBiological, DiffEqJump, OrdinaryDiffEq, Plots

using TimerOutputs

plotlyjs()

# Create a TimerOutput, this is the main type that keeps track of everything.
const to = TimerOutput()

# parameters
method = RSSA()
nsteps = 1e5                  # time to simulate to
networkname = "tester"
tf = 10.
#speciesf = path to species initial population file
#rxsf = path to reaction network file


# get the RSSA reaction network
reset_timer!(to)
@timeit to "netgen" prn = loadrxnetwork(RSSANetwork(), networkname, speciesf, rxsf; printrxs = false)
rn = prn.rn 
initialpop = prn.u₀
println("network parsed")

# get the BioNetGen reaction network
@timeit to "bionetgen" prnbng = loadrxnetwork(BNGNetwork(),string(networkname,"bng"), bngfname) 
rnbng = prnbng.rn; u₀ = prnbng.u₀; p = prnbng.p; shortsymstosyms = prnbng.symstonames;

u0 = round.(Int, u₀)
println(typeof(u0))

# one simulation
@timeit to "addjumps" addjumps!(rn,build_regular_jumps=false, minimal_jumps=true)
println("added jumps")
@timeit to "dprob" rprob = DiscreteProblem(rn, u0, (0.,tf))
@timeit to "jprob" rjprob = JumpProblem(rprob, method, rn; save_positions=(false,false))
@timeit to "solve" rsol = solve(rjprob, SSAStepper(), saveat=tf/1000)
show(to)


# bng file
reset_timer!(to)
@timeit to "addjumps" addjumps!(rnbng,build_regular_jumps=false, minimal_jumps=true)
@timeit to "dprob" bprob = DiscreteProblem(rnbng, u0, (0.,tf), p)
@timeit to "jprob" bjprob = JumpProblem(bprob, method, rnbng; save_positions=(false,false))
@timeit to "solve" bsol = solve(bjprob, SSAStepper(), saveat=tf/1000)
show(to)

# plot
#labs = [String(sym) for i=1:1, sym in rn.syms];
#plot(sol.t, sol[:,:]', linetype=:steppost, labels=labs)
#plot(sol.t, sol[:,:]', labels=labs)

rasyk = sum(rsol[prnbng.groupstoids[:Activated_Syk],:], dims=1)
basyk = sum(bsol[prnbng.groupstoids[:Activated_Syk],:], dims=1)
plot(rsol.t, rasyk', label=:rsol)
plot!(bsol.t, basyk', label=:bsol)


