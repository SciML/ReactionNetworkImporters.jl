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
    return only(@species $(S)(t))
end

abstract type NetworkFileFormat end

# exported data types
#struct RSSANetwork <: NetworkFileFormat end
struct BNGNetwork <: NetworkFileFormat end

"""
    ParsedReactionNetwork(rn::ReactionSystem; u0 = nothing, p = nothing, varstonames = nothing, groupstosyms = nothing)

A container for storing a parsed reaction network along with its associated metadata.

# Fields
- `rn`, a Catalyst `ReactionSystem`. **Note** this system is *not* marked complete by
  default (see the [Catalyst docs](https://catalyst.sciml.ai/) for a discussion of
  completeness of systems).
- `u0`, a `Dict` mapping initial condition symbolic variables to numeric values and/or
  symbolic expressions.
- `p`, a `Dict` mapping parameter symbolic variables to numeric values and/or symbolic
  expressions.
- `varstonames`, a `Dict` mapping the internal symbolic variable of a species used in the
  generated `ReactionSystem` to a `String` generated from the name in the .net file. This is
  necessary as BioNetGen can generate exceptionally long species names, involving characters
  that lead to malformed species names when used with `Catalyst`.
- `groupstosyms`, a `Dict` mapping the `String`s representing names for any groups defined
  in the BioNetGen file to the corresponding symbolic variable representing the
  `ModelingToolkit` symbolic observable associated with the group.
"""
struct ParsedReactionNetwork
    "Catalyst Network"
    rn::ReactionSystem

    "Dict mapping initial condition symbolic variables to values."
    u0::Any

    "Dict mapping parameter symbolic variables to values."
    p::Any

    "Dict mapping symbolic variable for species names to full string for species name"
    varstonames::Any

    "Dict from group name (as string) to corresponding symbolic variable"
    groupstosyms::Any
end
function ParsedReactionNetwork(
        rn::ReactionSystem; u0 = nothing, p = nothing,
        varstonames = nothing, groupstosyms = nothing
    )
    return ParsedReactionNetwork(rn, u0, p, varstonames, groupstosyms)
end

export BNGNetwork, MatrixNetwork, ParsedReactionNetwork, ComplexMatrixNetwork

# parsers
include("parsing_routines_bngnetworkfiles.jl")
include("parsing_routines_matrixnetworks.jl")

export loadrxnetwork

# Overload ensuring that u0 and u₀ can be used interchangeably.
# (introduced when the u₀ field was changed to u0)
# Should be deleted whenever u₀ is fully deprecated.

# Ensures that `prnbng.u₀` works.
function Base.getproperty(prnbng::ParsedReactionNetwork, name::Symbol)
    if name === :u₀
        return getfield(prnbng, :u0)
    else
        return getfield(prnbng, name)
    end
end

# Ensures that `prnbng.u₀ = ...` works.
function Base.setproperty!(prnbng::ParsedReactionNetwork, name::Symbol, x)
    if name === :u₀
        return setfield!(prnbng, :u0, x)
    else
        return setfield!(prnbng, name, x)
    end
end

end # module
