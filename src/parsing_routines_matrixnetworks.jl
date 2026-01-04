"""
Given substrate and product stoichiometric matrices, rate-expressions and lists
of species and parameters, return a ReactionSystem that describes the chemical
reaction network.

Assumed that the substrate and product stoichiometry matrices are stored as
numspecies by numrxs matrices, with entry (i,j) giving the stoichiometric
coefficient of species i within rx j.
"""
struct MatrixNetwork{S, T, U, V, W, X} <: NetworkFileFormat
    """The symbolic expressions for each reaction rate."""
    rateexprs::S

    """substoich[i,j] = is the substrate stoichiometric coefficient matrix"""
    substoich::T

    """prodstoich[i,j] is the product stoichiometric coefficient matrix"""
    prodstoich::U

    """species in the network """
    species::V

    """ Parameters """
    params::W

    """independent variable, time """
    t::X
end
function MatrixNetwork(
        rateexprs, substoich, prodstoich; species = Any[], params = Any[],
        t = nothing
    )
    return MatrixNetwork(rateexprs, substoich, prodstoich, species, params, t)
end

# for dense matrices
"""
    loadrxnetwork(mn::MatrixNetwork; name = gensym(:ReactionSystem))

Converts a `MatrixNetwork` into a `ParsedReactionNetwork` by constructing a
`ReactionSystem`.

# Arguments
- `mn::MatrixNetwork`: A `MatrixNetwork` object containing the stoichiometric matrices, rate
  expressions, species, and parameters.
- `name::Symbol`: (Optional) Name for the resulting `ReactionSystem`. Defaults to a
  generated symbol.

# Returns
A `ParsedReactionNetwork` 

# Notes
- The `MatrixNetwork` must have substrate (`substoich`) and product (`prodstoich`)
  stoichiometric matrices of the same size.
- The stoichiometric matrices are assumed to be `numspecies x numrxs`, where each entry `(i,
  j)` represents the stoichiometric coefficient of species `i` in reaction `j`.
- If the `species` field in `MatrixNetwork` is empty, species symbols are automatically
  generated.
- The `t` field specifies the independent variable (e.g., time). If not provided, a default
  time variable is used.

# Example
```julia
using Catalyst

# Define a MatrixNetwork
rateexprs = [1.0, 2.0]
substoich = [1 0; 0 1]
prodstoich = [0 1; 1 0]
species = [@species A(t), B(t)]
params = [@parameters k1, k2]
mn = MatrixNetwork(rateexprs, substoich, prodstoich; species = species, params = params)

# Convert to a ParsedReactionNetwork
parsed_network = loadrxnetwork(mn, name = :MyReactionSystem)
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

    return ParsedReactionNetwork(ReactionSystem(rxs, t, species, mn.params; name = name))
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

    return ParsedReactionNetwork(ReactionSystem(rxs, t, species, mn.params; name = name))
end

"""
Given complex stoichiometric matrix ,incidence matrix ,rate-expressions and list of species,parameters
, return a ReactionSystem that describes the chemical reaction network

Notes:

  - The column of complex stoichiometric represents composition of reaction complexes,
    with positive entries of size num_of_species by num_of_complexes, where
    the non-zero positive entries in the k'th column denote stoichiometric
    coefficients of the species participating in the k'th reaction complex.

  - The complex incidence matrix, is number of complexes by number of reactions with
    Bᵢⱼ = -1, if the i'th complex is the substrate of the j'th reaction,
    1, if the i'th complex is the product of the j'th reaction,
    0, otherwise
"""
struct ComplexMatrixNetwork{S, T, U, V, W, X} <: NetworkFileFormat
    """The symbolic expressions for each reaction rate."""
    rateexprs::S

    """stoichmat[i,j] = is the stoichiometric coefficient in the j'th reaction for the i'th species"""
    stoichmat::T

    """incidencemat[i,j] is the incidence matrix with incidencemat[i,j] = -1, if the i'th complex
    is the substrate of the j'th reaction, 1, if the i'th complex is the product
    of the j'th reaction, 0, otherwise"""
    incidencemat::U  # all elements are integer always in incidence matrix

    """species in the network """
    species::V

    """ Parameters """
    params::W

    """independent variable, time """
    t::X
end
function ComplexMatrixNetwork(
        rateexprs, stoichmat, incidencemat; species = Any[],
        params = Any[], t = nothing
    )
    return ComplexMatrixNetwork(rateexprs, stoichmat, incidencemat, species, params, t)
end

# for Dense matrices version
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

    return ParsedReactionNetwork(ReactionSystem(rxs, t, species, cmn.params; name = name))
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

    return ParsedReactionNetwork(ReactionSystem(rxs, t, species, cmn.params; name = name))
end
