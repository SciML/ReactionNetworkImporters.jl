using ReactionNetworkImporters, Aqua
@testset "Aqua" begin
    Aqua.find_persistent_tasks_deps(ReactionNetworkImporters)
    Aqua.test_ambiguities(ReactionNetworkImporters, recursive = false)
    Aqua.test_deps_compat(ReactionNetworkImporters)
    Aqua.test_piracies(ReactionNetworkImporters)
    Aqua.test_project_extras(ReactionNetworkImporters)
    Aqua.test_stale_deps(ReactionNetworkImporters)
    Aqua.test_unbound_args(ReactionNetworkImporters)
    Aqua.test_undefined_exports(ReactionNetworkImporters)
end
