"""
    For loading networks with stoichiometry stored in matrices.
    Assumed that the substrate and product stoichiometry matrices
    are stored as numspecies by numrxs matrices, with entry (i,j)
    giving the stoichiometric coefficient of species i within rx j.
"""


""" 
    For dense matrices
"""
function loadrxnetwork(ft::MatrixNetwork, networkname::String, params::Vector{Symbol}, 
                       rateexprs::AbstractVector, substoich::AbstractMatrix, prodstoich::AbstractMatrix)

    sz = size(substoich)
    @assert sz == size(prodstoich)
    numspecs = sz[1]
    numrxs = sz[2]
    

    # create the network
    rn = eval(Meta.parse("@empty_reaction_network $networkname"))

    # create the species
    foreach(i -> addspecies!(rn, Symbol("S",i)), 1:numspecs)

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
            (scoef > 0) && push!(subs, Symbol("S",i) => scoef)

            pcoef = prodstoich[i,j]
            (pcoef > 0) && push!(prods, Symbol("S",i) => pcoef)
        end

        addreaction!(rn, rateexprs[j], Tuple(subs), Tuple(prods))
    end

    ParsedReactionNetwork(rn, nothing)
end
