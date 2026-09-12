using ReactionNetworkImporters, BenchmarkTools
using Catalyst, OrdinaryDiffEqTsit5

const SUITE = BenchmarkGroup()

datadir = joinpath(@__DIR__, "..", "data")

# =============================================================================
# Network file loading (BioNetGen .net → ReactionSystem)
# =============================================================================

SUITE["load"] = BenchmarkGroup()

SUITE["load"]["birth_death"] = @benchmarkable loadrxnetwork(
    BNGNetwork(), joinpath($datadir, "nullrxs", "birth-death.net")
)
SUITE["load"]["repressilator"] = @benchmarkable loadrxnetwork(
    BNGNetwork(), joinpath($datadir, "repressilator", "Repressilator.net")
)
SUITE["load"]["matrix_network"] = @benchmarkable loadrxnetwork(
    MatrixNetwork([1.0], reshape([1], 1, 1), reshape([0], 1, 1));
    name = :bench_matrix
)

# =============================================================================
# Problem generation from an imported network
# =============================================================================

rnbng = loadrxnetwork(BNGNetwork(), joinpath(datadir, "nullrxs", "birth-death.net"))
rn = complete(rnbng)

SUITE["problem"] = BenchmarkGroup()

SUITE["problem"]["ode_problem"] = @benchmarkable ODEProblem(
    $rn, Float64[], (0.0, 10.0), Float64[]
)
SUITE["problem"]["complete"] = @benchmarkable complete($rnbng)
