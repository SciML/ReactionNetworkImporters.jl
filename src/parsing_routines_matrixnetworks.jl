"""
    For loading networks with stoichiometry stored in matrices.
    Assumed that the substrate and product stoichiometry matrices
    are stored as numspecies by numrxs matrices, with entry (i,j)
    giving the stoichiometric coefficient of species i within rx j.
"""


""" 
    For dense matrices
"""
function loadrxnetwork(ft::MatrixNetwork, networkname::String, 
                        rateexprs::AbstractVector, 
                        substoich::AbstractMatrix, 
                        prodstoich::AbstractMatrix; 
                        species=Symbol[], 
                        params=Symbol[])

    sz = size(substoich)
    @assert sz == size(prodstoich)
    numspecs = sz[1]
    numrxs = sz[2]

    # create the network
    rn = eval(Meta.parse("@empty_reaction_network $networkname"))

    # create the species if none passed in
    isempty(species) && (species = [Symbol("S",i) for i=1:numspecs])
    foreach(i -> addspecies!(rn, species[i]), 1:numspecs)    

    # create the parameters
    foreach(param -> addparam!(rn, param), params)

    # create the reactions
    subs = Vector{Pair{Symbol,Int}}()
    prods = Vector{Pair{Symbol,Int}}()
    for j = 1:numrxs
        empty!(subs)
        empty!(prods)

        # stoich
        for i = 1:numspecs
            scoef = substoich[i,j]
            (scoef > 0) && push!(subs, species[i] => scoef)

            pcoef = prodstoich[i,j]
            (pcoef > 0) && push!(prods, species[i] => pcoef)
        end

        addreaction!(rn, rateexprs[j], Tuple(subs), Tuple(prods))
    end

    ParsedReactionNetwork(rn, nothing)
end

""" 
    For sparse matrices
"""
function loadrxnetwork(ft::MatrixNetwork, networkname::String, 
                        rateexprs::AbstractVector, 
                        substoich::SparseMatrixCSC, 
                        prodstoich::SparseMatrixCSC; 
                        species=Symbol[], 
                        params=Symbol[])

    sz = size(substoich)
    @assert sz == size(prodstoich)
    numspecs = sz[1]
    numrxs = sz[2]

    # create the network
    rn = eval(Meta.parse("@empty_reaction_network $networkname"))

    # create the species if none passed in
    isempty(species) && (species = [Symbol("S",i) for i=1:numspecs])
    foreach(i -> addspecies!(rn, species[i]), 1:numspecs)    

    # create the parameters
    foreach(param -> addparam!(rn, param), params)

    # create the reactions
    subs = Vector{Pair{Symbol,Int}}()
    prods = Vector{Pair{Symbol,Int}}()
    srows = rowvals(substoich)
    svals = nonzeros(substoich)
    prows = rowvals(prodstoich)
    pvals = nonzeros(prodstoich)
    for j = 1:numrxs
        empty!(subs)
        empty!(prods)

        for ir in nzrange(substoich, j)
           i     = srows[ir]
           scoef = svals[ir]
           (scoef > 0) && push!(subs, species[i] => scoef)
        end

        for ir in nzrange(prodstoich, j)
            i     = prows[ir]
            pcoef = pvals[ir]
            (pcoef > 0) && push!(prods, species[i] => pcoef)
         end 

         addreaction!(rn, rateexprs[j], Tuple(subs), Tuple(prods))
     end

    ParsedReactionNetwork(rn, nothing)
end
