"""
    For loading networks with stoichiometry stored in matrices.
    Assumed that the substrate and product stoichiometry matrices
    are stored as numspecies by numrxs matrices, with entry (i,j)
    giving the stoichiometric coefficient of species i within rx j.
"""


""" 
    For dense matrices
"""
function loadrxnetwork(::MatrixNetwork, 
                        rateexprs::AbstractVector, 
                        substoich::AbstractMatrix, 
                        prodstoich::AbstractMatrix; 
                        species::AbstractVector=Any[], 
                        params::AbstractVector=Any[])

    sz = size(substoich)
    @assert sz == size(prodstoich)
    numspecs = sz[1]
    numrxs = sz[2]

    # create the network
    rn = make_empty_network()        
    t  = ModelingToolkit.get_iv(rn)

    # create the species if none passed in        
    isempty(species) && (species = [funcsym(:S,t,i) for i=1:numspecs])
    foreach(s -> addspecies!(rn, s, disablechecks=true), species)

    # create the parameters
    foreach(p -> addparam!(rn, p, disablechecks=true), params)

    # create the reactions
    # we need to create new vectors each time as the ReactionSystem
    # takes ownership of them
    for j = 1:numrxs
        subs    = Any[]
        sstoich = Vector{eltype(substoich)}()
        prods   = Any[]
        pstoich = Vector{eltype(prodstoich)}()
    
        # stoich
        for i = 1:numspecs
            scoef = substoich[i,j]
            if (scoef > zero(scoef)) 
                push!(subs, species[i])
                push!(sstoich, scoef)
            end

            pcoef = prodstoich[i,j]
            if (pcoef > zero(pcoef)) 
                push!(prods, species[i])
                push!(pstoich, pcoef)
            end
        end

        addreaction!(rn, Reaction(rateexprs[j], subs, prods, sstoich, pstoich))
    end

    ParsedReactionNetwork(rn, nothing)
end

""" 
    For sparse matrices
"""
function loadrxnetwork(ft::MatrixNetwork,
                        rateexprs::AbstractVector, 
                        substoich::SparseMatrixCSC, 
                        prodstoich::SparseMatrixCSC; 
                        species::AbstractVector=Any[], 
                        params::AbstractVector=Any[])

    sz = size(substoich)
    @assert sz == size(prodstoich)
    numspecs = sz[1]
    numrxs = sz[2]

    # create the network
    rn = make_empty_network()
    t  = ModelingToolkit.get_iv(rn)

    # create the species if none passed in
    isempty(species) && (species = [funcsym(:S,t,i) for i=1:numspecs])
    foreach(s -> addspecies!(rn, s, disablechecks=true), species)

    # create the parameters
    foreach(p -> addparam!(rn, p, disablechecks=true), params)

    # create the reactions
    srows = rowvals(substoich)
    svals = nonzeros(substoich)
    prows = rowvals(prodstoich)
    pvals = nonzeros(prodstoich)
    for j = 1:numrxs
        subs    = Any[]
        sstoich = Vector{eltype(substoich)}()
        prods   = Any[]
        pstoich = Vector{eltype(prodstoich)}()

        for ir in nzrange(substoich, j)
           i     = srows[ir]
           scoef = svals[ir]
           if scoef > zero(scoef)
                push!(subs, species[i])
                push!(sstoich, scoef)
           end
        end

        for ir in nzrange(prodstoich, j)
            i     = prows[ir]
            pcoef = pvals[ir]
            if pcoef > zero(pcoef)
                push!(prods, species[i])
                push!(pstoich, pcoef)
            end
        end

        addreaction!(rn, Reaction(rateexprs[j], subs, prods, sstoich, pstoich))
     end

    ParsedReactionNetwork(rn, nothing)
end


"""
    For loading networks with complex stoichiomatix ,incidence matrix.

Notes:
- Assumed that column of complex stoichiomatrix represents composition of reaction complexes,
  with positive entries of size num_of_species x num_of_complexes, where
  the non-zero positive entries in the kth column denote stoichiometric
  coefficients of the species participating in the kth reaction complex.

- The complex incidence matrix, is number of complexes by number of reactions with
  Bᵢⱼ = -1, if the i'th complex is the substrate of the j'th reaction,
         1, if the i'th complex is the product of the j'th reaction,
         0, otherwise
"""

struct ComplexMatrixNetwork{S <: AbstractVector,T <: Matrix,
                U <:Matrix{Int} ,V <: AbstractVector,W <: AbstractVector}
    rateexprs::S
    cxstoichmat::T
    in_mat::U  # all elements are integer always in incidence matrix
    species::V
    params::W
end
ComplexMatrixNetwork(rateexprs,cxstoichmat,in_mat; species=Any[],params=Any[]) =
           ComplexMatrixNetwork(rateexprs,cxstoichmat,in_mat,species,params)

struct ComplexSparseMatrixNetwork{S <: AbstractVector,T <: SparseMatrixCSC,
                U <: SparseMatrixCSC{Int,Int} ,V <: AbstractVector,W <: AbstractVector}
    rateexprs::S
    cxstoichmat::T
    in_mat::U  # all elements are integer always in incidence matrix
    species::V
    params::W
end
ComplexSparseMatrixNetwork(rateexprs,cxstoichmat,in_mat;species=Any[],params=Any[]) =
           ComplexSparseMatrixNetwork(rateexprs,cxstoichmat,in_mat,species,params)


# for Dense matrices version
function loadrxnetwork(cmn::T) where {T <: ComplexMatrixNetwork}
    numspecs, numcomp = size(cmn.cxstoichmat)
    @assert all(>=(0),cmn.cxstoichmat)
    @assert numcomp == size(cmn.in_mat, 1)
    @assert all(∈([-1,0,1]),cmn.in_mat)
    numrxs = size(cmn.in_mat, 2)

    ModelingToolkit.@parameters t

    isempty(cmn.species) ? species = [funcsym(:S,t,i) for i = 1:numspecs] : species = cmn.species

    rn = Vector{Reaction}(undef,numrxs)
    sc_ind = argmin(cmn.in_mat, dims=1)  # cartesian indices of substrate complexes
    pc_ind = argmax(cmn.in_mat, dims=1)  # cartesian indicies of products complexes

    for i ∈ 1:numrxs

        # substrate index for i'th reaction in species(rn)
        ss_ind = findall(!iszero, @view cmn.cxstoichmat[:,sc_ind[i][1]])

        # products index for i'th reaction in species(rn)
        ps_ind = findall(!iszero, @view cmn.cxstoichmat[:,pc_ind[i][1]])

        if isempty(ss_ind) && !isempty(ps_ind)
            rn[i] = Reaction(cmn.rateexprs[i], nothing, species[ps_ind],
                        nothing, cmn.cxstoichmat[ps_ind,pc_ind[i][1]])

        elseif !isempty(ss_ind) &&  isempty(ps_ind)
            rn[i] = Reaction(cmn.rateexprs[i], species[ss_ind], nothing,
                        cmn.cxstoichmat[ss_ind,sc_ind[i][1]], nothing)
        else
            rn[i] = Reaction(cmn.rateexprs[i], species[ss_ind],species[ps_ind],
                        cmn.cxstoichmat[ss_ind,sc_ind[i][1]],
                        cmn.cxstoichmat[ps_ind,pc_ind[i][1]])
        end
    end
    @named rs = ReactionSystem(rn,t,species, cmn.params)
    return ParsedReactionNetwork(rs,nothing)
end

# for sparse matrices version
function loadrxnetwork(cmn::T) where {T <: ComplexSparseMatrixNetwork}
    numspecs, numcomp = size(cmn.cxstoichmat)
    @assert all(>=(0),cmn.cxstoichmat)
    @assert numcomp == size(cmn.in_mat, 1)
    @assert all(∈([-1,0,1]),cmn.in_mat)
    numrxs = size(cmn.in_mat, 2)

    cp = cmn.cxstoichmat;   inmat= cmn.in_mat

    ModelingToolkit.@parameters t

    isempty(cmn.species) ? species = [funcsym(:S,t,i) for i = 1:numspecs] : species = cmn.species

    rn = Vector{Reaction}(undef,numrxs)
    sc_ind = findall(x -> x == -1, inmat)  # cartesian indices of substrate complexes
    pc_ind = findall(x -> x == 1, inmat)  # cartesian indicies of products complexes


    for i ∈ 1:numrxs
        # substrate index for i'th reaction in species(rn)
        ss_ind = rowvals(cp)[cp.colptr[sc_ind[i][1]]:cp.colptr[sc_ind[i][1]+1]-1]
        # products index for i'th reaction in species(rn)
        ps_ind = rowvals(cp)[cp.colptr[pc_ind[i][1]]:cp.colptr[pc_ind[i][1]+1]-1]

        if isempty(ss_ind) && !isempty(ps_ind)
            rn[i] = Reaction(cmn.rateexprs[i], nothing, species[ps_ind],
                        nothing, cp.nzval[cp.colptr[pc_ind[i][1]]:cp.colptr[pc_ind[i][1]+1]-1])

        elseif !isempty(ss_ind) &&  isempty(ps_ind)
            rn[i] = Reaction(cmn.rateexprs[i], species[ss_ind], nothing,
                        cp.nzval[cp.colptr[sc_ind[i][1]]:cp.colptr[sc_ind[i][1]+1]-1], nothing)
        else
            rn[i] = Reaction(cmn.rateexprs[i], species[ss_ind],species[ps_ind],
                        cp.nzval[cp.colptr[sc_ind[i][1]]:cp.colptr[sc_ind[i][1]+1]-1],
                        cp.nzval[cp.colptr[pc_ind[i][1]]:cp.colptr[pc_ind[i][1]+1]-1])
        end
    end
    @named rs = ReactionSystem(rn,t,species, cmn.params)
    return ParsedReactionNetwork(rs,nothing)
end
