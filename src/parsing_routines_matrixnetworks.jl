"""
    For loading networks with stoichiometry stored in matrices.
    Assumed that the substrate and product stoichiometry matrices
    are stored as numspecies by numrxs matrices, with entry (i,j)
    giving the stoichiometric coefficient of species i within rx j.
"""


""" 
    For dense matrices
"""

function loadrxnetwork(ft::MatrixNetwork, networkname::String, rateexprs::AbstractVector, 
                        substoich::AbstractMatrix, prodstoich::AbstractMatrix) 
    
    loadrxnetwork(ft, networkname, Symbol[], rateexprs, substoich, prodstoich)
end

function loadrxnetwork(ft::MatrixNetwork, networkname::String, params::Vector{Symbol}, 
                       rateexprs::AbstractVector, substoich::AbstractMatrix, prodstoich::AbstractMatrix)

    sz = size(substoich)
    @assert sz == size(prodstoich)
    numspecs = sz[1]

    #create species names
    speciessyms = Vector{Symbol}(undef, numspecs)
    foreach(i -> speciessyms[i] = Symbol("S",i), 1:numspecs)

    loadrxnetwork(ft, networkname, speciessyms, params, rateexprs, substoich, prodstoich)
end


function loadrxnetwork(ft::MatrixNetwork, networkname::String, speciessyms::Vector{Symbol}, 
                        params::Vector{Symbol}, rateexprs::AbstractVector, 
                        substoich::AbstractMatrix, prodstoich::AbstractMatrix)

    sz = size(substoich)
    @assert sz == size(prodstoich)
    numspecs = sz[1]
    @assert sz[1] == length(speciessyms)
    numrxs = sz[2]

    # create the network
    rn = eval(Meta.parse("@empty_reaction_network $networkname"))

    # create the species
    foreach(i -> addspecies!(rn, speciessyms[i]), 1:numspecs)

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
            (scoef > 0) && push!(subs, speciessyms[i] => scoef)

            pcoef = prodstoich[i,j]
            (pcoef > 0) && push!(prods, speciessyms[i] => pcoef)
        end

        addreaction!(rn, rateexprs[j], Tuple(subs), Tuple(prods))
    end

    ParsedReactionNetwork(rn, nothing)
end
