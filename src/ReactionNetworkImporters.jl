module ReactionNetworkImporters

using DataStructures, Catalyst, SparseArrays
using Symbolics: operation, unwrap

# creates a ModelingToolkit function-like Symbol
# can then do stuff like
# @parameters t
# S₁ = funcsym(S,1)
# u = S₁(t)
function funcsym(S::Symbol, t, args...)
    S = Symbol(S, args...)
    return (@variables $(S)(t))[1]
end

abstract type NetworkFileFormat end

# exported data types
#struct RSSANetwork <: NetworkFileFormat end
struct BNGNetwork <: NetworkFileFormat end

struct ParsedReactionNetwork
    "Catalyst Network"
    rn::ReactionSystem

    "Dict mapping initial condition symbolic variables to values."
    u₀::Any

    "Dict mapping parameter symbolic variables to values."
    p::Any

    "Dict mapping symbolic variable for species names to full string for species name"
    varstonames::Any

    "Dict from group name (as string) to corresponding symbolic variable"
    groupstosyms::Any
end
ParsedReactionNetwork(rn::ReactionSystem; u₀ = nothing, p = nothing, varstonames = nothing, groupstosyms = nothing) = ParsedReactionNetwork(rn,
                                                                                                                                            u₀,
                                                                                                                                            p,
                                                                                                                                            varstonames,
                                                                                                                                            groupstosyms)

export BNGNetwork, MatrixNetwork, ParsedReactionNetwork, ComplexMatrixNetwork

# parsers
include("parsing_routines_bngnetworkfiles.jl")
include("parsing_routines_matrixnetworks.jl")

export loadrxnetwork

end # module
