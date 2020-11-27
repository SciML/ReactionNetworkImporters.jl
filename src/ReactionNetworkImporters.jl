module ReactionNetworkImporters

using DataStructures, Catalyst, SparseArrays

# creates a ModelingToolkit function-like Symbol
# can then do stuff like
# @parameters t
# S₁ = funcsym(S,1)
# u = S₁(t)
function funcsym(S::Symbol, args...)
    Num(Variable{ModelingToolkit.FnType{Tuple{Any},Real}}(S,args...))
end

abstract type NetworkFileFormat end

# exported data types
#struct RSSANetwork <: NetworkFileFormat end
struct BNGNetwork <: NetworkFileFormat end
struct MatrixNetwork <: NetworkFileFormat end

struct ParsedReactionNetwork    
    "Catalyst Network"
    rn::ReactionSystem

    "Initial Conditions"
    u₀

    "Parameters"
    p

    "Expressions for the Parameters"
    paramexprs

    "Dict from `Variable` in species(rn) to full string for species name"
    varstonames

    "Dict from lumped species name (as string) to group of species ids"
    groupstoids

end
ParsedReactionNetwork(rn::ReactionSystem, u₀; p=nothing, paramexprs=nothing, varstonames=nothing, groupstoids=nothing) = 
                        ParsedReactionNetwork(rn, u₀, p, paramexprs, varstonames, groupstoids)

export RSSANetwork, BNGNetwork, MatrixNetwork, ParsedReactionNetwork

# parsers
#include("parsing_routines_rssafiles.jl")
include("parsing_routines_bngnetworkfiles.jl")
include("parsing_routines_matrixnetworks.jl")

export loadrxnetwork

end # module
