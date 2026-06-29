using SciMLTesting, ReactionNetworkImporters, Test

run_qa(
    ReactionNetworkImporters;
    explicit_imports = true,
    ei_kwargs = (;
        # `unwrap` is imported from `Symbolics` (its public, exported re-export) but
        # its `Base.which` owner is `SymbolicUtils`, where it is non-public; allow the
        # import from the public re-exporter.
        all_explicit_imports_via_owners = (; ignore = (:unwrap,)),
        all_explicit_imports_are_public = (; ignore = (:unwrap,)),
        all_qualified_accesses_are_public = (;
            ignore = (
                :U0Map,        # public in Catalyst (declared `public`; only invisible to EI on Julia 1.10)
                :ParameterMap, # public in Catalyst (declared `public`; only invisible to EI on Julia 1.10)
                :parse,        # Base.Meta.parse, non-public stdlib name
                :eval,         # Base.eval, non-public Base name; required to eval BNG exprs in a target module
            ),
        ),
    ),
    # Many names from the heavy DSL deps (Catalyst/Symbolics/SymbolicUtils/
    # DataStructures/SparseArrays) are used implicitly via `using`; making them all
    # explicit is a risky mass refactor. Tracked in SciML/ReactionNetworkImporters.jl#179.
    ei_broken = (:no_implicit_imports,)
)
