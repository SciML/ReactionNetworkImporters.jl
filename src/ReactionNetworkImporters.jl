module ReactionNetworkImporters

import Catalyst
using Catalyst: Reaction, ReactionSystem, @parameters, @species, @variables
using HypergeometricFunctions: _₁F₁
using OrderedCollections: OrderedDict
using PrecompileTools: @compile_workload, @setup_workload
using SciMLPublic: @public
using SparseArrays: SparseMatrixCSC, nonzeros, nzrange, rowvals
using SymbolicUtils: getmetadata, hasmetadata, setmetadata, simplify
using Symbolics: Equation, Num, unwrap
using TermInterface: operation

function funcsym(S::Symbol, t, args...)
    S = Symbol(S, args...)
    return only(@species $(S)(t))
end

"""
    NetworkFileFormat

Abstract developer interface for reaction-network file-format selectors.

# Extension rules
- Define a concrete subtype for each supported external format.
- Extend [`loadrxnetwork`](@ref) for that subtype and the format's input.
- Return a Catalyst `ReactionSystem`; do not mutate parser-global state.
- Keep format-specific parser helpers private. Downstream packages should only
  rely on this type and `loadrxnetwork`.

# Fields
Concrete subtypes may store format-specific configuration. The built-in
selectors are fieldless.

# Examples
```jldoctest
using ReactionNetworkImporters

struct ExampleFormat <: ReactionNetworkImporters.NetworkFileFormat end
ExampleFormat() isa ReactionNetworkImporters.NetworkFileFormat

# output
true
```
"""
abstract type NetworkFileFormat end
@public NetworkFileFormat

# exported data types
#struct RSSANetwork <: NetworkFileFormat end
"""
    BNGNetwork()

File-format selector for BioNetGen `.net` reaction network files.

# Arguments
None.

# Keywords
None.

# Fields
None.

Pass `BNGNetwork()` to [`loadrxnetwork`](@ref) to parse a BioNetGen `.net` file
into a Catalyst `ReactionSystem`.

# Examples
```jldoctest
using ReactionNetworkImporters

BNGNetwork() isa ReactionNetworkImporters.NetworkFileFormat

# output
true
```
"""
struct BNGNetwork <: NetworkFileFormat end

### System-Level Metadata Key Types and Accessors ###

"""
    VarsToNames

Metadata key for storing BNG variable-to-full-name mappings on a
`ReactionSystem`. Stores a `Dict` mapping internal symbolic variables (both
dynamic species and constant-species parameters) to their full BNG name
strings from the .net file.

# Fields
None. This fieldless type is used as a metadata key.

See also: [`has_varstonames`](@ref), [`get_varstonames`](@ref), [`set_varstonames`](@ref)
"""
struct VarsToNames end

"""
    has_varstonames(rs::ReactionSystem)

Return whether `rs` has a `VarsToNames` metadata entry.

# Arguments
- `rs::ReactionSystem`: System to inspect.

# Returns
- `Bool`: `true` when the mapping exists.
"""
has_varstonames(rs::ReactionSystem) = hasmetadata(rs, VarsToNames)

"""
    get_varstonames(rs::ReactionSystem)

Return the `VarsToNames` mapping from `rs`, or `nothing` when it is absent.

# Arguments
- `rs::ReactionSystem`: System to inspect.

# Returns
- `Dict` or `nothing`: Mapping from system variables to BioNetGen names.
"""
get_varstonames(rs::ReactionSystem) = getmetadata(rs, VarsToNames, nothing)

"""
    set_varstonames(rs::ReactionSystem, m)

Return a new `ReactionSystem` with the `VarsToNames` metadata set to `m`.

# Arguments
- `rs::ReactionSystem`: System to update.
- `m`: Mapping from system variables to BioNetGen names.

# Returns
- `ReactionSystem`: Copy of `rs` carrying `m`; `rs` is not mutated.
"""
set_varstonames(rs::ReactionSystem, m) = setmetadata(rs, VarsToNames, m)

"""
    GroupsToSyms

Metadata key for storing BNG group-name-to-symbol mappings on a
`ReactionSystem`. Stores a `Dict` mapping `String` group names to
corresponding symbolic observable variables.

# Fields
None. This fieldless type is used as a metadata key.

See also: [`has_groupstosyms`](@ref), [`get_groupstosyms`](@ref), [`set_groupstosyms`](@ref)
"""
struct GroupsToSyms end

"""
    has_groupstosyms(rs::ReactionSystem)

Return whether `rs` has a `GroupsToSyms` metadata entry.

# Arguments
- `rs::ReactionSystem`: System to inspect.

# Returns
- `Bool`: `true` when the mapping exists.
"""
has_groupstosyms(rs::ReactionSystem) = hasmetadata(rs, GroupsToSyms)

"""
    get_groupstosyms(rs::ReactionSystem)

Return the `GroupsToSyms` mapping from `rs`, or `nothing` when it is absent.

# Arguments
- `rs::ReactionSystem`: System to inspect.

# Returns
- `Dict` or `nothing`: Mapping from BioNetGen group names to observables.
"""
get_groupstosyms(rs::ReactionSystem) = getmetadata(rs, GroupsToSyms, nothing)

"""
    set_groupstosyms(rs::ReactionSystem, m)

Return a new `ReactionSystem` with the `GroupsToSyms` metadata set to `m`.

# Arguments
- `rs::ReactionSystem`: System to update.
- `m`: Mapping from BioNetGen group names to observables.

# Returns
- `ReactionSystem`: Copy of `rs` carrying `m`; `rs` is not mutated.
"""
set_groupstosyms(rs::ReactionSystem, m) = setmetadata(rs, GroupsToSyms, m)

export BNGNetwork, MatrixNetwork, ComplexMatrixNetwork
export VarsToNames, GroupsToSyms
export has_varstonames, get_varstonames, set_varstonames
export has_groupstosyms, get_groupstosyms, set_groupstosyms

# parsers
include("parsing_routines_bngnetworkfiles.jl")
include("parsing_routines_matrixnetworks.jl")
include("precompilation.jl")

export loadrxnetwork

end # module
