"""
Parsing routines for BioNetGen .net files. Note, this only handles a subset
of the BioNetGen .net file format. In particular, any functions in the file
that use non-Julia expressions (i.e. like conditional logic) or time-dependent
features, will definitely not work. 
"""

"""
Seek to a line matching `start_string`.
"""
function seek_to_block(lines, idx, start_string)
    while lines[idx] != start_string
        idx += 1        
        (idx > length(lines)) && error("Block: ", start_string, " was never found.")
    end
    idx += 1
end

# for swapping out BioNetGen symbols
const REPLACEMENT_DICT = Dict( :ln => :log )
function recursive_replace!(expr::Any, replace_requests::Dict{Symbol,Symbol}=REPLACEMENT_DICT)
    if expr isa Symbol
        haskey(replace_requests,expr) && return replace_requests[expr]
    elseif expr isa Expr
        for i = 1:length(expr.args)
            expr.args[i] = recursive_replace!(expr.args[i], replace_requests)
        end
    end
    return expr
end

const PARAM_BLOCK_START = "begin parameters"
const PARAM_BLOCK_END = "end parameters"
function parse_params(ft::BNGNetwork, lines, idx)

    idx = seek_to_block(lines, idx, PARAM_BLOCK_START)

    # parse params
    pvals = []
    ptoids = OrderedDict{Symbol,Int}()
    while lines[idx] != PARAM_BLOCK_END
        vals = split(lines[idx])
        pidx = parse(Int,vals[1])
        psym = Symbol(vals[2])
        ptoids[psym] = pidx    

        # value could be an expression
        pval = Meta.parse(vals[3])
        (pval isa Expr) && recursive_replace!(pval)

        push!(pvals,pval)
        idx += 1
        (idx > length(lines)) && error("Block: ", PARAM_BLOCK_END, " was never found.")
    end

    ptoids,pvals,idx
end

const CHARS_TO_STRIP = ",~!()."
stripchars = (s) -> replace(s, Regex("[$CHARS_TO_STRIP]") => "")
hasstripchars = (s) -> occursin(Regex("[$CHARS_TO_STRIP]"), s)
const SPECIES_BLOCK_START = "begin species"
const SPECIES_BLOCK_END = "end species"
const MAX_SYM_LEN = 8
function parse_species(ft::BNGNetwork, lines, idx)

    idx = seek_to_block(lines, idx, SPECIES_BLOCK_START)

    # parse species
    shortsymstoids = OrderedDict{Symbol,Int}()
    shortsymstosyms = OrderedDict{Symbol,Symbol}()
    u0exprs   = Vector{Any}()
    while lines[idx] != SPECIES_BLOCK_END
        vals = split(lines[idx])
        sidx = parse(Int,vals[1])
        sym = Symbol(vals[2])
        shortsym = length(vals[2]) < MAX_SYM_LEN ? Symbol(stripchars(vals[2])) : Symbol(string("S$sidx"))
        shortsymstoids[shortsym] = sidx
        shortsymstosyms[shortsym] = sym
        push!(u0exprs,  Meta.parse(vals[3]))
        idx += 1
        (idx > length(lines)) && error("Block: ", SPECIES_BLOCK_END, " was never found.")
    end

    u0exprs,shortsymstoids,shortsymstosyms,idx
end

const REACTIONS_BLOCK_START = "begin reactions"
const REACTIONS_BLOCK_END = "end reactions"
function parse_reactions!(rn, ft, lines, idx, idstosyms)

    idx = seek_to_block(lines, idx, REACTIONS_BLOCK_START)
    cntdict = Dict{Int,Int}()
    while lines[idx] != REACTIONS_BLOCK_END
        vals        = split(lines[idx])        
        reactantids = (parse(Int,rid) for rid in split(vals[2],","))
        productids  = (parse(Int,pid) for pid in split(vals[3],","))        
        rateexpr    = Meta.parse(vals[4])
        (typeof(rateexpr) <: Expr) && recursive_replace!(rateexpr)

        any(iszero,reactantids) && (@assert length(reactantids)==1 "Found more than one reactant with id 0")
        any(iszero,productids) && (@assert length(productids)==1 "Found more than one product with id 0")

        # reactants and correct for higher-order rate rescalings by BioNetGen
        empty!(cntdict)
        foreach(rid -> (rid>0) && (haskey(cntdict,rid) ? (cntdict[rid] += 1) : (cntdict[rid]=1)), reactantids)
        scalefactor = isempty(cntdict) ? 0 : prod(factorial, values(cntdict))
        (!iszero(scalefactor)) && (rateexpr = :($scalefactor * $rateexpr))
        rstoich = Tuple(idstosyms[rid] => cnt for (rid,cnt) in cntdict)

        # product stoichiometry
        empty!(cntdict)
        foreach(pid -> (pid>0) && (haskey(cntdict,pid) ? (cntdict[pid] += 1) : (cntdict[pid]=1)), productids)
        pstoich = Tuple(idstosyms[pid] => cnt for (pid,cnt) in cntdict)

        # create the reaction
        addreaction!(rn, rateexpr, rstoich, pstoich)

        idx += 1
        (idx > length(lines)) && error("Block: ", REACTIONS_BLOCK_END, " was never found.")
    end

    idx
end

const GROUPS_BLOCK_START = "begin groups"
const GROUPS_BLOCK_END = "end groups"
function parse_groups(ft::BNGNetwork, lines, idx, idstoshortsyms, rn)
    
    idx = seek_to_block(lines, idx, GROUPS_BLOCK_START)
    namestoids = Dict{Symbol,Vector{Int}}()
    while lines[idx] != GROUPS_BLOCK_END
        vals = split(lines[idx])
        name = Symbol(vals[2])

        # map from BioNetGen id to reaction_network id
        ids = [rn.syms_to_ints[idstoshortsyms[parse(Int,val)]] for val in split(vals[3],",")]        

        namestoids[name] = ids
        idx += 1
        (idx > length(lines)) && error("Block: ", GROUPS_BLOCK_END, " was never found.")
    end

    namestoids,idx
end

function exprs_to_nums(ptoids, pvals, u0exprs)
    p    = zeros(Float64, length(ptoids))
    pmod = Module()
    for (psym,pid) in ptoids
        p[pid] = Base.eval(pmod, :($psym = $(pvals[pid])))
    end

    u0 = zeros(Float64, length(u0exprs))
    for (i,u0expr) in enumerate(u0exprs)
        u0[i] = Base.eval(pmod, :($u0expr))
    end

    p,u0
end


# for parsing a subset of the BioNetGen .net file format
function loadrxnetwork(ft::BNGNetwork, networkname::String, rxfilename; kwargs...)

    file  = open(rxfilename, "r");
    lines = readlines(file)
    idx   = 1
    
    print("Parsing parameters...")
    ptoids,pvals,idx = parse_params(ft, lines, idx)
    println("done")

    rn = eval(Meta.parse("@empty_reaction_network $networkname"))

    print("Adding parameters...")
    foreach(psym -> addparam!(rn, psym), keys(ptoids))
    println("done")

    print("Parsing species...")
    u0exprs,shortsymstoids,shortsymstosyms,idx = parse_species(ft, lines, idx)
    println("done")

    # map from species id to short sym
    print("Adding species...")
    idstoshortsyms = Vector{Symbol}(undef,length(shortsymstoids))
    for (k,v) in shortsymstoids
        idstoshortsyms[v] = k
        addspecies!(rn, k)
    end
    @assert species(rn) == idstoshortsyms "species(rn) noteq to idstoshortsyms"
    println("done")

    print("Parsing and adding reactions...")
    idx = parse_reactions!(rn, ft, lines, idx, idstoshortsyms)
    println("done")

    print("Parsing groups...")
    groupstoids,idx = parse_groups(ft, lines, idx, idstoshortsyms, rn)
    println("done")

    close(file)

    # get numeric values for parameters and u₀ 
    p,u₀ = exprs_to_nums(ptoids, pvals, u0exprs)
    @assert all( speciesmap(rn)[sym] == i for (i,sym) in enumerate(idstoshortsyms) )

    ParsedReactionNetwork(rn, u₀; p = p, paramexprs = pvals, symstonames = shortsymstosyms, groupstoids = groupstoids)
end