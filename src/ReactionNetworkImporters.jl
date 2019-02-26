module ReactionNetworkImporters

using DataStructures, DiffEqBiological

abstract type NetworkFileFormat end

# exported data types
struct RSSANetwork <: NetworkFileFormat end
struct BNGNetwork <: NetworkFileFormat end
struct ParsedReactionNetwork    
    "DiffEqBiological Network"
    rn

    "Initial Conditions"
    u₀

    "Parameters"
    p

    "Parameter Expressions"
    paramexprs

    "Dict from short sym in rn.syms to full sym for species name"
    symstonames

    "Dict from lumped species name (symbol) to group of species ids"
    groupstoids

    "string representation of the generated reaction_network macro"
    rnstr
end
ParsedReactionNetwork(rn, u₀; p=nothing, 
                              paramexprs=nothing,
                              symstonames=nothing, 
                              groupstoids=nothing,
                              rnstr=nothing) = 
                              ParsedReactionNetwork(rn, u₀, p, 
                                                    paramexprs, 
                                                    symstonames, 
                                                    groupstoids,
                                                    rnstr)

export RSSANetwork, BNGNetwork, ParsedReactionNetwork

# parsers
include("parsing_routines_rssafiles.jl")
include("parsing_routines_bngnetworkfiles.jl")

export loadrxnetwork

end # module
