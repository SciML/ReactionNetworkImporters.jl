using ReactionNetworkImporters, Aqua
using Test
@testset "Aqua" begin
    Aqua.find_persistent_tasks_deps(ReactionNetworkImporters)
    Aqua.test_ambiguities(ReactionNetworkImporters, recursive = false)
    # The deps/weakdeps sub-checks pass; only the extras sub-check fails because
    # `Pkg` is in [extras] without a [compat] entry. Tracked in
    # https://github.com/SciML/ReactionNetworkImporters.jl/issues/175
    Aqua.test_deps_compat(ReactionNetworkImporters, check_extras = false)
    @test_broken false  # Aqua deps_compat extras: Pkg in [extras] lacks a [compat] entry — tracked in https://github.com/SciML/ReactionNetworkImporters.jl/issues/175
    Aqua.test_piracies(ReactionNetworkImporters)
    Aqua.test_project_extras(ReactionNetworkImporters)
    Aqua.test_stale_deps(ReactionNetworkImporters)
    Aqua.test_unbound_args(ReactionNetworkImporters)
    Aqua.test_undefined_exports(ReactionNetworkImporters)
end
