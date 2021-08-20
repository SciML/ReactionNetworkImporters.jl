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
"""
for dense matrices
"""
function LoadReacCompNetwork(rateexprs::AbstractVector,
                        compstoichmat::AbstractMatrix,
                        incidencemat::AbstractMatrix;
                        species::AbstractVector = Any[],
                        params::AbstractVector = Any[])

    numspecs, numcomp = size(compstoichmat)
    @assert all(>=(0),compstoichmat)
    @assert numcomp == size(incidencemat, 1)
    @assert all(∈([-1,0,1]),incidencemat)
    numrxs = size(incidencemat, 2)

    ModelingToolkit.@parameters t
    isempty(species) && (species = [funcsym(:S,t,i) for i = 1:numspecs])

    rn = Vector{Reaction}(undef,numrxs)
    sub_indices = argmin(incidencemat, dims=1)     # cartesian indices of substrate complexes
    prod_indices = argmax(incidencemat, dims=1)    # cartesian indicies of products complexes

    for i ∈ 1:numrxs
        subStoichInd = getindex.(findall(!iszero, compstoichmat[:,sub_indices[i][1]]))
        prodStoichInd = getindex.(findall(!iszero, compstoichmat[:,prod_indices[i][1]]))

        if subStoichInd == Int64[] &&  prodStoichInd != Int64[]
            rn[i] = Reaction(rateexprs[i], nothing, species[prodStoichInd],
                        nothing, compstoichmat[prodStoichInd,prod_indices[i][1]])

        elseif subStoichInd != Int64[] &&  prodStoichInd == Int64[]
            rn[i] = Reaction(rateexprs[i], species[subStoichInd], nothing,
                        compstoichmat[subStoichInd,sub_indices[i][1]], nothing)
        else
            rn[i] = Reaction(rateexprs[i], species[subStoichInd],species[prodStoichInd],
                        compstoichmat[subStoichInd,sub_indices[i][1]],
                        compstoichmat[prodStoichInd,prod_indices[i][1]])
        end
    end

    @named rs = ReactionSystem(rn, t, species, params)
end

"""
for sparse matrices
"""
function LoadReacCompNetwork(rateexprs::AbstractVector,
                        compstoichmat::SparseMatrixCSC,
                        incidencemat::SparseMatrixCSC;
                        species::AbstractVector = Any[],
                        params::AbstractVector = Any[])

    numspecs, numcomp = size(compstoichmat)
    @assert all(>=(0),compstoichmat)
    @assert numcomp == size(incidencemat, 1)
    @assert all(∈([-1,0,1]),incidencemat)
    numrxs = size(incidencemat, 2)

    ModelingToolkit.@parameters t
    isempty(species) && (species = [funcsym(:S,t,i) for i = 1:numspecs])

    rn = Vector{Reaction}(undef,numrxs)

    sub_indices = argmin(incidencemat, dims=1)     # cartesian indices of substrate complexes
    prod_indices = argmax(incidencemat, dims=1)    # cartesian indicies of products complexes

    for i ∈ 1:numrxs
        subStoichInd = getindex.(findall(!iszero, compstoichmat[:,sub_indices[i][1]]))
        prodStoichInd = getindex.(findall(!iszero, compstoichmat[:,prod_indices[i][1]]))

        if subStoichInd == Int64[] &&  prodStoichInd != Int64[]
            rn[i] = Reaction(rateexprs[i], nothing, species[prodStoichInd],
                        nothing, compstoichmat[prodStoichInd,prod_indices[i][1]].nzval)

        elseif subStoichInd != Int64[] &&  prodStoichInd == Int64[]
            rn[i] = Reaction(rateexprs[i], species[subStoichInd], nothing,
                        compstoichmat[subStoichInd,sub_indices[i][1]].nzval, nothing)
        else
            rn[i] = Reaction(rateexprs[i], species[subStoichInd],species[prodStoichInd],
                        compstoichmat[subStoichInd,sub_indices[i][1]].nzval,
                        compstoichmat[prodStoichInd,prod_indices[i][1]].nzval)
        end
    end

    @named rs = ReactionSystem(rn, t, species, params)
end

