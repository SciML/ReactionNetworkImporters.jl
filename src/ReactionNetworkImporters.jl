module ReactionNetworkImporters

using DataStructures, Catalyst, SparseArrays
using Symbolics: operation, unwrap

# creates a ModelingToolkit function-like Symbol
# can then do stuff like
# @parameters t
# S₁ = funcsym(S,1)
# u = S₁(t)
function funcsym(S::Symbol, t, args...)
    S = Symbol(S,args...)
    (@variables $(S)(t))[1]
end

abstract type NetworkFileFormat end

# exported data types
#struct RSSANetwork <: NetworkFileFormat end
struct BNGNetwork <: NetworkFileFormat end

struct ParsedReactionNetwork    
    "Catalyst Network"
    rn::ReactionSystem

    "Initial Conditions"
    u₀

    "Parameters"
    p

    "Expressions for the Parameters"
    paramexprs

    "Dict from symbolic variable in species(rn) to full string for species name"
    varstonames

    "Dict from group name (as string) to corresponding symbolic variable"
    groupstosyms

end
ParsedReactionNetwork(rn::ReactionSystem, u₀; p=nothing, paramexprs=nothing, varstonames=nothing, groupstosyms=nothing) = 
                        ParsedReactionNetwork(rn, u₀, p, paramexprs, varstonames, groupstosyms)

export BNGNetwork, MatrixNetwork, ParsedReactionNetwork, ComplexMatrixNetwork

# parsers
#include("parsing_routines_rssafiles.jl")
include("parsing_routines_bngnetworkfiles.jl")
include("parsing_routines_matrixnetworks.jl")

export loadrxnetwork

end # module
