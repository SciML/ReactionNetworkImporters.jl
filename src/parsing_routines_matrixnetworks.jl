"""
    MatrixNetwork(rateexprs, substoich, prodstoich; species = Any[], params = Any[], t = nothing)

Input representation for constructing a Catalyst reaction network from
substrate and product stoichiometry matrices.

# Arguments
- `rateexprs::AbstractVector`: One symbolic or numeric rate expression per reaction.
- `substoich::AbstractMatrix`: Species-by-reaction substrate coefficients.
- `prodstoich::AbstractMatrix{Int}`: Species-by-reaction product coefficients;
  it must have the same size as `substoich`.

# Keywords
- `species::AbstractVector = Any[]`: Symbolic species in row order. Empty uses
  generated species with the chosen independent variable.
- `params::AbstractVector = Any[]`: Symbolic parameters referenced by `rateexprs`.
- `t = nothing`: Independent variable. `nothing` uses `Catalyst.default_t()`.

# Fields
- `rateexprs`: Rate expressions, one for each reaction column.
- `substoich`: Substrate stoichiometry matrix.
- `prodstoich`: Product stoichiometry matrix.
- `species`: Symbolic species in matrix-row order.
- `params`: Symbolic parameters available to rate expressions.
- `t`: Independent variable or `nothing`.

# Examples
```jldoctest
using ReactionNetworkImporters

network = MatrixNetwork([1.0], reshape([1], 1, 1), reshape([0], 1, 1))
network isa MatrixNetwork

# output
true
```
"""
struct MatrixNetwork{S, T, U, V, W, X} <: NetworkFileFormat
    rateexprs::S
    substoich::T
    prodstoich::U
    species::V
    params::W
    t::X
end
function MatrixNetwork(
        rateexprs, substoich, prodstoich; species = Any[], params = Any[],
        t = nothing
    )
    return MatrixNetwork(rateexprs, substoich, prodstoich, species, params, t)
end

"""
    loadrxnetwork(mn::MatrixNetwork; name = gensym(:ReactionSystem))

Convert a `MatrixNetwork` into a Catalyst `ReactionSystem`.

# Arguments
- `mn::MatrixNetwork`: Input matrix representation satisfying `MatrixNetwork`'s
  dimensionality rules.

# Keywords
- `name::Symbol = gensym(:ReactionSystem)`: Name assigned to the resulting system.

# Returns
- `ReactionSystem`: An incomplete Catalyst system. Call `complete` before using
  it to construct a SciML problem.

# Rules
- `substoich` and `prodstoich` must be species-by-reaction matrices of equal size.
- Each rate expression corresponds to one reaction column.
- An empty `species` vector requests generated species symbols.

# Examples
```jldoctest
using ReactionNetworkImporters

network = MatrixNetwork([1.0], reshape([1], 1, 1), reshape([0], 1, 1))
system = loadrxnetwork(network; name = :decay)
nameof(system)

# output
:decay
```
"""
function loadrxnetwork(
        mn::MatrixNetwork{S, T, U, V, W, X};
        name = gensym(:ReactionSystem)
    ) where {
        S <: AbstractVector,
        T <: Matrix, U <: Matrix{Int},
        V <: AbstractVector,
        W <: AbstractVector, X <: Any,
    }
    sz = size(mn.substoich)
    @assert sz == size(mn.prodstoich)
    numspecs = sz[1]
    numrxs = sz[2]

    t = (mn.t === nothing) ? Catalyst.default_t() : mn.t
    species = isempty(mn.species) ? [funcsym(:S, t, i) for i in 1:numspecs] : mn.species

    # create the reactions
    # we need to create new vectors each time as the ReactionSystem
    # takes ownership of them
    rxs = Vector{Reaction}(undef, numrxs)
    for j in 1:numrxs
        subs = Any[]
        sstoich = Vector{eltype(mn.substoich)}()
        prods = Any[]
        pstoich = Vector{eltype(mn.prodstoich)}()

        # stoich
        for i in 1:numspecs
            scoef = mn.substoich[i, j]
            if (scoef > zero(scoef))
                push!(subs, species[i])
                push!(sstoich, scoef)
            end

            pcoef = mn.prodstoich[i, j]
            if (pcoef > zero(pcoef))
                push!(prods, species[i])
                push!(pstoich, pcoef)
            end
        end

        rxs[j] = Reaction(mn.rateexprs[j], subs, prods, sstoich, pstoich)
    end

    return ReactionSystem(rxs, t, species, mn.params; name)
end

# for sparse matrices
function loadrxnetwork(
        mn::MatrixNetwork{S, T, U, V, W, X};
        name = gensym(:ReactionSystem)
    ) where {
        S <: AbstractVector,
        T <: SparseMatrixCSC,
        U <:
        SparseMatrixCSC{Int, Int},
        V <: AbstractVector,
        W <: AbstractVector, X <: Any,
    }
    sz = size(mn.substoich)
    @assert sz == size(mn.prodstoich)
    numspecs = sz[1]
    numrxs = sz[2]

    t = (mn.t === nothing) ? Catalyst.default_t() : mn.t
    species = isempty(mn.species) ? [funcsym(:S, t, i) for i in 1:numspecs] : mn.species

    # create the reactions
    srows = rowvals(mn.substoich)
    svals = nonzeros(mn.substoich)
    prows = rowvals(mn.prodstoich)
    pvals = nonzeros(mn.prodstoich)
    rxs = Vector{Reaction}(undef, numrxs)
    for j in 1:numrxs
        subs = Any[]
        sstoich = Vector{eltype(mn.substoich)}()
        prods = Any[]
        pstoich = Vector{eltype(mn.prodstoich)}()

        for ir in nzrange(mn.substoich, j)
            i = srows[ir]
            scoef = svals[ir]
            if scoef > zero(scoef)
                push!(subs, species[i])
                push!(sstoich, scoef)
            end
        end

        for ir in nzrange(mn.prodstoich, j)
            i = prows[ir]
            pcoef = pvals[ir]
            if pcoef > zero(pcoef)
                push!(prods, species[i])
                push!(pstoich, pcoef)
            end
        end

        rxs[j] = Reaction(mn.rateexprs[j], subs, prods, sstoich, pstoich)
    end

    return ReactionSystem(rxs, t, species, mn.params; name)
end

"""
    ComplexMatrixNetwork(rateexprs, stoichmat, incidencemat; species = Any[], params = Any[], t = nothing)

Input representation for constructing a Catalyst reaction network from complex
stoichiometry and complex-incidence matrices.

# Arguments
- `rateexprs::AbstractVector`: One symbolic or numeric rate expression per reaction.
- `stoichmat::AbstractMatrix`: Species-by-complex nonnegative stoichiometric coefficients.
- `incidencemat::AbstractMatrix{Int}`: Complex-by-reaction matrix whose entries
  are `-1`, `0`, or `1`.

# Keywords
- `species::AbstractVector = Any[]`: Symbolic species in row order. Empty uses
  generated species with the chosen independent variable.
- `params::AbstractVector = Any[]`: Symbolic parameters referenced by `rateexprs`.
- `t = nothing`: Independent variable. `nothing` uses `Catalyst.default_t()`.

# Fields
- `rateexprs`: Rate expressions, one for each reaction column.
- `stoichmat`: Species-by-complex stoichiometry matrix.
- `incidencemat`: Complex-by-reaction incidence matrix.
- `species`: Symbolic species in matrix-row order.
- `params`: Symbolic parameters available to rate expressions.
- `t`: Independent variable or `nothing`.

# Examples
```jldoctest
using ReactionNetworkImporters

network = ComplexMatrixNetwork([1.0], [1 0], reshape([-1, 1], 2, 1))
network isa ComplexMatrixNetwork

# output
true
```
"""
struct ComplexMatrixNetwork{S, T, U, V, W, X} <: NetworkFileFormat
    rateexprs::S
    stoichmat::T
    incidencemat::U
    species::V
    params::W
    t::X
end
function ComplexMatrixNetwork(
        rateexprs, stoichmat, incidencemat; species = Any[],
        params = Any[], t = nothing
    )
    return ComplexMatrixNetwork(rateexprs, stoichmat, incidencemat, species, params, t)
end

"""
    loadrxnetwork(cmn::ComplexMatrixNetwork; name = gensym(:ReactionSystem))

Convert a `ComplexMatrixNetwork` into a Catalyst `ReactionSystem`.

# Arguments
- `cmn::ComplexMatrixNetwork`: Complex-matrix representation satisfying the
  documented stoichiometry and incidence rules.

# Keywords
- `name::Symbol = gensym(:ReactionSystem)`: Name assigned to the resulting system.

# Returns
- `ReactionSystem`: An incomplete Catalyst system. Call `complete` before using
  it to construct a SciML problem.

# Rules
- `stoichmat` is species-by-complex and has nonnegative entries.
- `incidencemat` is complex-by-reaction with entries in `(-1, 0, 1)`.
- Every reaction column needs one substrate complex (`-1`) and one product
  complex (`1`).

# Examples
```jldoctest
using ReactionNetworkImporters

network = ComplexMatrixNetwork([1.0], [1 0], reshape([-1, 1], 2, 1))
system = loadrxnetwork(network; name = :conversion)
nameof(system)

# output
:conversion
```
"""
function loadrxnetwork(
        cmn::ComplexMatrixNetwork{S, T, U, V, W, X};
        name = gensym(:ReactionSystem)
    ) where {
        S <: AbstractVector,
        T <: Matrix, U <: Matrix{Int},
        V <: AbstractVector,
        W <: AbstractVector, X <: Any,
    }
    numspecs, numcomp = size(cmn.stoichmat)
    @assert all(>=(0), cmn.stoichmat)
    @assert numcomp == size(cmn.incidencemat, 1)
    @assert all(∈([-1, 0, 1]), cmn.incidencemat)
    numrxs = size(cmn.incidencemat, 2)

    t = (cmn.t === nothing) ? Catalyst.default_t() : cmn.t
    species = isempty(cmn.species) ? [funcsym(:S, t, i) for i in 1:numspecs] : cmn.species

    rxs = Vector{Reaction}(undef, numrxs)
    sc_ind = argmin(cmn.incidencemat, dims = 1)  # cartesian indices of substrate complexes
    pc_ind = argmax(cmn.incidencemat, dims = 1)  # cartesian indices of products complexes
    for i in 1:numrxs

        # substrate index for i'th reaction in species(rn)
        ss_ind = findall(!iszero, @view cmn.stoichmat[:, sc_ind[i][1]])

        # products index for i'th reaction in species(rn)
        ps_ind = findall(!iszero, @view cmn.stoichmat[:, pc_ind[i][1]])

        if isempty(ss_ind) && !isempty(ps_ind)
            rxs[i] = Reaction(
                cmn.rateexprs[i], nothing, species[ps_ind],
                nothing, cmn.stoichmat[ps_ind, pc_ind[i][1]]
            )

        elseif !isempty(ss_ind) && isempty(ps_ind)
            rxs[i] = Reaction(
                cmn.rateexprs[i], species[ss_ind], nothing,
                cmn.stoichmat[ss_ind, sc_ind[i][1]], nothing
            )
        else
            rxs[i] = Reaction(
                cmn.rateexprs[i], species[ss_ind], species[ps_ind],
                cmn.stoichmat[ss_ind, sc_ind[i][1]],
                cmn.stoichmat[ps_ind, pc_ind[i][1]]
            )
        end
    end

    return ReactionSystem(rxs, t, species, cmn.params; name)
end

# for sparse matrices version
function loadrxnetwork(
        cmn::ComplexMatrixNetwork{S, T, U, V, W, X};
        name = gensym(:ReactionSystem)
    ) where {
        S <: AbstractVector,
        T <: SparseMatrixCSC,
        U <:
        SparseMatrixCSC{Int, Int},
        V <: AbstractVector,
        W <: AbstractVector, X <: Any,
    }
    numspecs, numcomp = size(cmn.stoichmat)
    @assert all(>=(0), cmn.stoichmat)
    @assert numcomp == size(cmn.incidencemat, 1)
    @assert all(∈([-1, 0, 1]), cmn.incidencemat)
    numrxs = size(cmn.incidencemat, 2)

    t = (cmn.t === nothing) ? Catalyst.default_t() : cmn.t
    species = isempty(cmn.species) ? [funcsym(:S, t, i) for i in 1:numspecs] : cmn.species

    rxs = Vector{Reaction}(undef, numrxs)
    sc_ind = argmin(cmn.incidencemat, dims = 1)  # cartesian indices of substrate complexes
    pc_ind = argmax(cmn.incidencemat, dims = 1)  # cartesian indices of products complexes
    rows = rowvals(cmn.stoichmat)
    vals = nonzeros(cmn.stoichmat)
    for i in 1:numrxs
        # substrate index for i'th reaction in species(rn)
        ss_ind = @view rows[nzrange(cmn.stoichmat, sc_ind[i][1])]
        # products index for i'th reaction in species(rn)
        ps_ind = @view rows[nzrange(cmn.stoichmat, pc_ind[i][1])]

        if isempty(ss_ind) && !isempty(ps_ind)
            rxs[i] = Reaction(
                cmn.rateexprs[i], nothing, species[ps_ind],
                nothing, vals[nzrange(cmn.stoichmat, pc_ind[i][1])]
            )

        elseif !isempty(ss_ind) && isempty(ps_ind)
            rxs[i] = Reaction(
                cmn.rateexprs[i], species[ss_ind], nothing,
                vals[nzrange(cmn.stoichmat, sc_ind[i][1])], nothing
            )
        else
            rxs[i] = Reaction(
                cmn.rateexprs[i], species[ss_ind], species[ps_ind],
                vals[nzrange(cmn.stoichmat, sc_ind[i][1])],
                vals[nzrange(cmn.stoichmat, pc_ind[i][1])]
            )
        end
    end

    return ReactionSystem(rxs, t, species, cmn.params; name)
end
