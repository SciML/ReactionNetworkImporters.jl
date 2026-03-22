module ReactionNetworkImporters

using DataStructures, Catalyst, SparseArrays
using Symbolics: operation, unwrap
using SymbolicUtils: hasmetadata, getmetadata, setmetadata

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

### System-Level Metadata Key Types and Accessors ###

"""
    VarsToNames

Metadata key for storing BNG species variable-to-full-name mappings on a
`ReactionSystem`. Stores a `Dict` mapping the internal symbolic variable of a
species to a `String` with its full name from the .net file.

See also: [`has_varstonames`](@ref), [`get_varstonames`](@ref), [`set_varstonames`](@ref)
"""
struct VarsToNames end

"""
    has_varstonames(rs::ReactionSystem)

Returns `true` if the `ReactionSystem` has a `VarsToNames` metadata entry.
"""
has_varstonames(rs::ReactionSystem) = hasmetadata(rs, VarsToNames)

"""
    get_varstonames(rs::ReactionSystem)

Returns the `VarsToNames` metadata from the `ReactionSystem`, or `nothing` if
not set.
"""
get_varstonames(rs::ReactionSystem) = getmetadata(rs, VarsToNames, nothing)

"""
    set_varstonames(rs::ReactionSystem, m)

Returns a **new** `ReactionSystem` with the `VarsToNames` metadata set to `m`.
The original system is not modified.
"""
set_varstonames(rs::ReactionSystem, m) = setmetadata(rs, VarsToNames, m)

"""
    GroupsToSyms

Metadata key for storing BNG group-name-to-symbol mappings on a
`ReactionSystem`. Stores a `Dict` mapping `String` group names to
corresponding symbolic observable variables.

See also: [`has_groupstosyms`](@ref), [`get_groupstosyms`](@ref), [`set_groupstosyms`](@ref)
"""
struct GroupsToSyms end

"""
    has_groupstosyms(rs::ReactionSystem)

Returns `true` if the `ReactionSystem` has a `GroupsToSyms` metadata entry.
"""
has_groupstosyms(rs::ReactionSystem) = hasmetadata(rs, GroupsToSyms)

"""
    get_groupstosyms(rs::ReactionSystem)

Returns the `GroupsToSyms` metadata from the `ReactionSystem`, or `nothing` if
not set.
"""
get_groupstosyms(rs::ReactionSystem) = getmetadata(rs, GroupsToSyms, nothing)

"""
    set_groupstosyms(rs::ReactionSystem, m)

Returns a **new** `ReactionSystem` with the `GroupsToSyms` metadata set to `m`.
The original system is not modified.
"""
set_groupstosyms(rs::ReactionSystem, m) = setmetadata(rs, GroupsToSyms, m)

export BNGNetwork, MatrixNetwork, ComplexMatrixNetwork
export VarsToNames, GroupsToSyms
export has_varstonames, get_varstonames, set_varstonames
export has_groupstosyms, get_groupstosyms, set_groupstosyms

# parsers
include("parsing_routines_bngnetworkfiles.jl")
include("parsing_routines_matrixnetworks.jl")

export loadrxnetwork

end # module
