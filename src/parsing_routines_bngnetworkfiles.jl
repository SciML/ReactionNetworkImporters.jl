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
const REPLACEMENT_DICT = Dict(:ln => :log)
function recursive_replace!(expr::Any,
                            replace_requests::Dict{Symbol, Symbol} = REPLACEMENT_DICT)
    if expr isa Symbol
        haskey(replace_requests, expr) && return replace_requests[expr]
    elseif expr isa Expr
        for i in 1:length(expr.args)
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
    ptoids = OrderedDict{Symbol, Int}()
    while lines[idx] != PARAM_BLOCK_END
        vals = split(lines[idx])
        pidx = parse(Int, vals[1])
        psym = Symbol(vals[2])
        ptoids[psym] = pidx

        # value could be an expression
        pval = Meta.parse(vals[3])
        (pval isa Expr) && recursive_replace!(pval)

        push!(pvals, pval)
        idx += 1
        (idx > length(lines)) && error("Block: ", PARAM_BLOCK_END, " was never found.")
    end

    ptoids, pvals, idx
end

const CHARS_TO_STRIP = ",~!()."
stripchars = (s) -> replace(s, Regex("[$CHARS_TO_STRIP]") => "")
hasstripchars = (s) -> occursin(Regex("[$CHARS_TO_STRIP]"), s)
const SPECIES_BLOCK_START = "begin species"
const SPECIES_BLOCK_END = "end species"
const MAX_SYM_LEN = 8
function make_shortsym(name, sidx)
    length(name) < MAX_SYM_LEN ? Symbol(stripchars(name)) : Symbol(string("S$sidx"))
end
function parse_species(ft::BNGNetwork, lines, idx)
    idx = seek_to_block(lines, idx, SPECIES_BLOCK_START)

    # parse species
    shortsymstoids = OrderedDict{Symbol, Int}()
    shortsymstosyms = OrderedDict{Symbol, Symbol}()
    u0exprs = Vector{Any}()
    while lines[idx] != SPECIES_BLOCK_END
        vals = split(lines[idx])
        sidx = parse(Int, vals[1])
        sym = Symbol(vals[2])
        shortsym = make_shortsym(vals[2], sidx)
        shortsymstoids[shortsym] = sidx
        shortsymstosyms[shortsym] = sym
        push!(u0exprs, Meta.parse(vals[3]))
        idx += 1
        (idx > length(lines)) && error("Block: ", SPECIES_BLOCK_END, " was never found.")
    end

    u0exprs, shortsymstoids, shortsymstosyms, idx
end

const REACTIONS_BLOCK_START = "begin reactions"
const REACTIONS_BLOCK_END = "end reactions"
function parse_reactions!(ft::BNGNetwork, specs, ps, t, lines, idx, idstosyms, opmod)
    idx = seek_to_block(lines, idx, REACTIONS_BLOCK_START)
    cntdict = Dict{Int, Int}()
    rxs = Reaction[]
    while lines[idx] != REACTIONS_BLOCK_END
        vals = split(lines[idx])
        reactantids = (parse(Int, rid) for rid in split(vals[2], ","))
        productids = (parse(Int, pid) for pid in split(vals[3], ","))
        rateexpr = Meta.parse(vals[4])

        # create a ModelingToolkitExpr from rateexpr
        if rateexpr isa Expr
            recursive_replace!(rateexpr)
            rate = Base.eval(opmod, rateexpr)
        elseif rateexpr isa Symbol
            rate = (@parameters $rateexpr)[1]
        else
            rate = rateexpr
        end

        any(iszero, reactantids) &&
            (@assert length(reactantids)==1 "Found more than one reactant with id 0")
        any(iszero, productids) &&
            (@assert length(productids)==1 "Found more than one product with id 0")

        # reactants and correct for higher-order rate rescalings by BioNetGen
        empty!(cntdict)
        foreach(rid -> (rid > 0) &&
                    (haskey(cntdict, rid) ? (cntdict[rid] += 1) : (cntdict[rid] = 1)),
                reactantids)
        scalefactor = isempty(cntdict) ? 0 : prod(factorial, values(cntdict))
        (!iszero(scalefactor)) && (rate = simplify(scalefactor * rate))
        rspecs = Vector{Any}(undef, length(cntdict))
        rstoich = Vector{Int}(undef, length(cntdict))
        for (i, (rid, cnt)) in enumerate(cntdict)
            rspecs[i] = funcsym(idstosyms[rid], t)
            rstoich[i] = cnt
        end

        # product stoichiometry
        empty!(cntdict)
        foreach(pid -> (pid > 0) &&
                    (haskey(cntdict, pid) ? (cntdict[pid] += 1) : (cntdict[pid] = 1)),
                productids)
        pspecs = Vector{Any}(undef, length(cntdict))
        pstoich = Vector{Int}(undef, length(cntdict))
        for (i, (pid, cnt)) in enumerate(cntdict)
            pspecs[i] = funcsym(idstosyms[pid], t)
            pstoich[i] = cnt
        end

        # create the reaction
        push!(rxs, Reaction(rate, rspecs, pspecs, rstoich, pstoich))

        idx += 1
        (idx > length(lines)) && error("Block: ", REACTIONS_BLOCK_END, " was never found.")
    end

    rxs, idx
end

const GROUPS_BLOCK_START = "begin groups"
const GROUPS_BLOCK_END = "end groups"
function parse_groups(ft::BNGNetwork, lines, idx, shortsymstosyms, idstoshortsyms, specs, t)
    idx = seek_to_block(lines, idx, GROUPS_BLOCK_START)
    namestosyms = Dict{String, Any}()
    obseqs = Equation[]
    while lines[idx] != GROUPS_BLOCK_END
        vals = split(lines[idx])
        name = Symbol(vals[2])
        obs = (@variables $name($t))[1]
        namestosyms[vals[2]] = obs

        # map from BioNetGen id to reaction_network id
        rhs = Num(0)
        groupterms = split(vals[3], ",")
        for groupterm in groupterms
            splitterm = split(groupterm, "*")
            if length(splitterm) == 1
                ssym = idstoshortsyms[parse(Int, splitterm[1])]
                rhs += funcsym(ssym, t)
            elseif length(splitterm) == 2
                ssym = idstoshortsyms[parse(Int, splitterm[2])]
                rhs += parse(Int, splitterm[1]) * funcsym(ssym, t)
            else
                error("Don't know how to handle term: $groupterm, appearing in a Group")
            end
        end
        push!(obseqs, Equation(obs, rhs))

        idx += 1
        (idx > length(lines)) && error("Block: ", GROUPS_BLOCK_END, " was never found.")
    end

    obseqs, namestosyms, idx
end

function exprs_to_defs(opmod, ptoids, pvals, specs, u0exprs)
    pmap = Dict()
    for (psym, pid) in ptoids
        pvar = getproperty(opmod, psym)
        parsedval = pvals[pid]
        if (parsedval isa Expr) || (parsedval isa Symbol)
            pval = Base.eval(opmod, :($parsedval))
        else
            @assert parsedval isa Number
            pval = parsedval
        end
        push!(pmap, pvar => pval)
    end

    u0map = Dict()
    for (i, u0expr) in enumerate(u0exprs)
        uvar = specs[i]
        if (u0expr isa Expr) || (u0expr isa Symbol)
            u0val = Base.eval(opmod, :($u0expr))
        else
            @assert u0expr isa Number
            u0val = u0expr
        end
        push!(u0map, uvar => u0val)
    end

    union(pmap, u0map), pmap, u0map
end

# for parsing a subset of the BioNetGen .net file format
function loadrxnetwork(ft::BNGNetwork, rxfilename; name = gensym(:ReactionSystem),
                       kwargs...)
    file = open(rxfilename, "r")
    lines = readlines(file)
    idx = 1
    t = Catalyst.default_t()

    print("Parsing parameters...")
    ptoids, pvals, idx = parse_params(ft, lines, idx)
    println("done")

    print("Creating parameters...")
    ps = [(@parameters $psym)[1] for psym in keys(ptoids)]
    println("done")

    print("Parsing species...")
    u0exprs, shortsymstoids, shortsymstosyms, idx = parse_species(ft, lines, idx)
    println("done")

    # map from species id to short sym
    print("Creating species...")
    idstoshortsyms = Vector{Symbol}(undef, length(shortsymstoids))
    specs = []
    for (k, v) in shortsymstoids
        idstoshortsyms[v] = k
        push!(specs, funcsym(k, t))
    end
    @assert all(s -> nameof(operation(unwrap(s[1]))) == s[2], zip(specs, idstoshortsyms)) "species ≂̸ to idstoshortsyms"
    println("done")

    # we evaluate all the parameters and species names in a module
    print("Creating species and parameters for evaluating expressions...")
    opmod = Module()
    Base.eval(opmod, :(using Catalyst))
    Base.eval(opmod, :(t = Catalyst.default_t()))
    for p in ps
        psym = nameof(p)
        Base.eval(opmod, :($(psym) = $p))
    end
    for s in specs
        ssym = nameof(operation(unwrap(s)))
        Base.eval(opmod, :($(ssym) = $s))
    end
    println("done")

    print("Parsing and adding reactions...")
    rxs, idx = parse_reactions!(ft, specs, ps, t, lines, idx, idstoshortsyms, opmod)
    println("done")

    print("Parsing groups...")
    obseqs, groupstosyms, idx = parse_groups(ft, lines, idx, shortsymstosyms,
                                             idstoshortsyms, specs, t)
    println("done")

    close(file)

    # setup default values / expressions for params and initial conditions
    defmap, pmap, u0map = exprs_to_defs(opmod, ptoids, pvals, specs, u0exprs)

    # build the model
    rn = ReactionSystem(rxs, t, specs, ps; name = name, observed = obseqs,
                        defaults = defmap, kwargs...)

    # get numeric values for parameters and u0
    sm = speciesmap(rn)
    @assert all(sm[funcsym(sym, t)] == i for (i, sym) in enumerate(idstoshortsyms))

    ParsedReactionNetwork(rn; u0 = u0map, p = pmap, varstonames = shortsymstosyms,
                          groupstosyms = groupstosyms)
end
