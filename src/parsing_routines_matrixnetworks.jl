"""
Given substrate and product stoichiometric matrices ,rate-expressions and list of species,parameters
, return a ReactionSystem that describes the chemical reaction network

Assumed that the substrate and product stoichiometry matrices
are stored as numspecies by numrxs matrices, with entry (i,j)
giving the stoichiometric coefficient of species i within rx j.
"""

struct MatrixNetwork{S ,T ,U ,V ,W,X} <: NetworkFileFormat
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
MatrixNetwork(rateexprs,substoich,prodstoich; species=Any[],params=Any[],t=nothing) =
           MatrixNetwork(rateexprs,substoich,prodstoich,species,params,t)

# for dense matrices
function loadrxnetwork(mn::MatrixNetwork{S,T,U,V,W,X}; name = gensym(:ReactionSystem)) where {S <: AbstractVector,
        T <: Matrix, U <: Matrix{Int}, V <: AbstractVector, W <: AbstractVector,X <: Any}

    sz = size(mn.substoich)
    @assert sz == size(mn.prodstoich)
    numspecs= sz[1]
    numrxs = sz[2]

    # create the network
    rn = make_empty_network(;name = name)
    t = (mn.t === nothing) ? (@variables t)[1] : mn.t

    # create the species if none passed in
    species = isempty(mn.species) ? [funcsym(:S,t,i) for i = 1:numspecs] : mn.species
    foreach(s -> addspecies!(rn, s, disablechecks=true), species)

    # create the parameters
    foreach(p -> addparam!(rn, p, disablechecks=true), mn.params)

    # create the reactions
    # we need to create new vectors each time as the ReactionSystem
    # takes ownership of them
    for j = 1:numrxs
        subs    = Any[]
        sstoich = Vector{eltype(mn.substoich)}()
        prods   = Any[]
        pstoich = Vector{eltype(mn.prodstoich)}()

        # stoich
        for i = 1:numspecs
            scoef = mn.substoich[i,j]
            if (scoef > zero(scoef))
                push!(subs, species[i])
                push!(sstoich, scoef)
            end

            pcoef = mn.prodstoich[i,j]
            if (pcoef > zero(pcoef))
                push!(prods, species[i])
                push!(pstoich, pcoef)
            end
        end

        addreaction!(rn, Reaction(mn.rateexprs[j], subs, prods, sstoich, pstoich))
    end

    ParsedReactionNetwork(rn, nothing)
end

# for sparse matrices
function loadrxnetwork(mn::MatrixNetwork{S,T,U,V,W,X}; name = gensym(:ReactionSystem)) where {S<:AbstractVector,
        T<:SparseMatrixCSC,U<:SparseMatrixCSC{Int,Int},V<:AbstractVector, W<:AbstractVector,X <: Any}
    sz = size(mn.substoich)
    @assert sz == size(mn.prodstoich)
    numspecs= sz[1]
    numrxs = sz[2]

    # create the network
    rn = make_empty_network(;name = name)
    t = (mn.t === nothing) ? (@variables t)[1] : mn.t

    # create the species if none passed in
    species = isempty(mn.species) ? [funcsym(:S,t,i) for i = 1:numspecs] : mn.species
    foreach(s -> addspecies!(rn, s, disablechecks=true), species)

    # create the parameters
    foreach(p -> addparam!(rn, p, disablechecks=true), mn.params)

    # create the reactions
    srows = rowvals(mn.substoich)
    svals = nonzeros(mn.substoich)
    prows = rowvals(mn.prodstoich)
    pvals = nonzeros(mn.prodstoich)
    for j = 1:numrxs
        subs    = Any[]
        sstoich = Vector{eltype(mn.substoich)}()
        prods   = Any[]
        pstoich = Vector{eltype(mn.prodstoich)}()

        for ir in nzrange(mn.substoich, j)
           i     = srows[ir]
           scoef = svals[ir]
           if scoef > zero(scoef)
                push!(subs, species[i])
                push!(sstoich, scoef)
           end
        end

        for ir in nzrange(mn.prodstoich, j)
            i     = prows[ir]
            pcoef = pvals[ir]
            if pcoef > zero(pcoef)
                push!(prods, species[i])
                push!(pstoich, pcoef)
            end
        end

        addreaction!(rn, Reaction(mn.rateexprs[j], subs, prods, sstoich, pstoich))
     end

    ParsedReactionNetwork(rn, nothing)
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

struct ComplexMatrixNetwork{S ,T ,U ,V ,W,X} <: NetworkFileFormat
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
ComplexMatrixNetwork(rateexprs,stoichmat,incidencemat; species=Any[],params=Any[],t=nothing) =
           ComplexMatrixNetwork(rateexprs,stoichmat,incidencemat,species,params,t)



# for Dense matrices version
function loadrxnetwork(cmn::ComplexMatrixNetwork{S,T,U,V,W,X}; name = gensym(:ReactionSystem)) where {S <: AbstractVector,
                T <: Matrix, U <: Matrix{Int}, V <: AbstractVector,
                                        W <: AbstractVector,X <: Any}
   numspecs, numcomp = size(cmn.stoichmat)
   @assert all(>=(0),cmn.stoichmat)
   @assert numcomp == size(cmn.incidencemat, 1)
   @assert all(∈([-1,0,1]),cmn.incidencemat)
   numrxs = size(cmn.incidencemat, 2)

   t = (cmn.t === nothing) ? (@variables t)[1] : cmn.t
   species = isempty(cmn.species) ? [funcsym(:S,t,i) for i = 1:numspecs] : cmn.species

   rn = Vector{Reaction}(undef,numrxs)
   sc_ind = argmin(cmn.incidencemat, dims=1)  # cartesian indices of substrate complexes
   pc_ind = argmax(cmn.incidencemat, dims=1)  # cartesian indicies of products complexes

   for i ∈ 1:numrxs

       # substrate index for i'th reaction in species(rn)
       ss_ind = findall(!iszero, @view cmn.stoichmat[:,sc_ind[i][1]])

       # products index for i'th reaction in species(rn)
       ps_ind = findall(!iszero, @view cmn.stoichmat[:,pc_ind[i][1]])

       if isempty(ss_ind) && !isempty(ps_ind)
           rn[i] = Reaction(cmn.rateexprs[i], nothing, species[ps_ind],
                       nothing, cmn.stoichmat[ps_ind,pc_ind[i][1]])

       elseif !isempty(ss_ind) &&  isempty(ps_ind)
           rn[i] = Reaction(cmn.rateexprs[i], species[ss_ind], nothing,
                       cmn.stoichmat[ss_ind,sc_ind[i][1]], nothing)
       else
           rn[i] = Reaction(cmn.rateexprs[i], species[ss_ind],species[ps_ind],
                       cmn.stoichmat[ss_ind,sc_ind[i][1]],
                       cmn.stoichmat[ps_ind,pc_ind[i][1]])
       end
   end
   rs = ReactionSystem(rn,t,species, cmn.params;name = name)
   return ParsedReactionNetwork(rs,nothing)
end

# for sparse matrices version
function loadrxnetwork(cmn::ComplexMatrixNetwork{S,T,U,V,W,X}; name = gensym(:ReactionSystem)) where {S<:AbstractVector,
                T<:SparseMatrixCSC,U<:SparseMatrixCSC{Int,Int},V<:AbstractVector,
                                                    W<:AbstractVector,X <: Any}
   numspecs, numcomp = size(cmn.stoichmat)
   @assert all(>=(0),cmn.stoichmat)
   @assert numcomp == size(cmn.incidencemat, 1)
   @assert all(∈([-1,0,1]),cmn.incidencemat)
   numrxs = size(cmn.incidencemat, 2)

   t = (cmn.t === nothing) ? (@variables t)[1] : cmn.t
   species = isempty(cmn.species) ? [funcsym(:S,t,i) for i = 1:numspecs] : cmn.species

   rn = Vector{Reaction}(undef,numrxs)
   sc_ind = argmin(cmn.incidencemat, dims=1)  # cartesian indices of substrate complexes
   pc_ind = argmax(cmn.incidencemat, dims=1)  # cartesian indicies of products complexes

   rows = rowvals(cmn.stoichmat)
   vals = nonzeros(cmn.stoichmat)
   for i ∈ 1:numrxs
       # substrate index for i'th reaction in species(rn)
       ss_ind = @view rows[nzrange(cmn.stoichmat, sc_ind[i][1])]
       # products index for i'th reaction in species(rn)
       ps_ind = @view rows[nzrange(cmn.stoichmat, pc_ind[i][1])]

       if isempty(ss_ind) && !isempty(ps_ind)
           rn[i] = Reaction(cmn.rateexprs[i], nothing, species[ps_ind],
                       nothing, vals[nzrange(cmn.stoichmat, pc_ind[i][1])])

       elseif !isempty(ss_ind) &&  isempty(ps_ind)
           rn[i] = Reaction(cmn.rateexprs[i], species[ss_ind], nothing,
                       vals[nzrange(cmn.stoichmat, sc_ind[i][1])], nothing)
       else
           rn[i] = Reaction(cmn.rateexprs[i], species[ss_ind],species[ps_ind],
                       vals[nzrange(cmn.stoichmat, sc_ind[i][1])],
                       vals[nzrange(cmn.stoichmat, pc_ind[i][1])])
       end
   end
   rs = ReactionSystem(rn,t,species, cmn.params;name = name)
   return ParsedReactionNetwork(rs,nothing)
end
