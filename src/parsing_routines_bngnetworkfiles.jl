"""
Parsing routines for BioNetGen .net files. Supports parameters, species,
reactions, groups (observables), functions (concentration/time-dependent rate
expressions), fixed species (`\$`-prefixed boundary conditions), `Sat` typed
rate laws, compartment-aware species naming, and BNG built-in functions
(`if`, `time`, `min`, `max`, `mratio`, etc.).

Unsupported features (error if encountered):
- `MM` and `Hill` typed rate laws
- `tfun()` tabular interpolation functions
"""

# ─────────────────────────────────────────────────────────────────────────────
# Block discovery
# ─────────────────────────────────────────────────────────────────────────────

"""
    find_block_boundaries(lines)

Scan all lines once and return a `Dict{String, Tuple{Int,Int}}` mapping each
`"begin <name>"` string to `(content_start, content_end)` line indices, where
`content_start` is the first line after the `begin` marker and `content_end`
is the line containing the corresponding `end` marker.
"""
function find_block_boundaries(lines)
    boundaries = Dict{String, Tuple{Int, Int}}()
    i = 1
    while i <= length(lines)
        stripped = strip(lines[i])
        if startswith(stripped, "begin ")
            block_start = String(stripped)
            end_marker = replace(block_start, "begin " => "end ", count = 1)
            content_start = i + 1
            j = content_start
            while j <= length(lines) && strip(lines[j]) != end_marker
                j += 1
            end
            (j > length(lines)) &&
                error("Block '$block_start' starting at line $i was never closed.")
            boundaries[block_start] = (content_start, j)
            i = j + 1
        else
            i += 1
        end
    end
    return boundaries
end

# ─────────────────────────────────────────────────────────────────────────────
# BNG → Julia symbol replacements and expression preprocessing
# ─────────────────────────────────────────────────────────────────────────────

const REPLACEMENT_DICT = Dict{Symbol, Any}(
    :ln => :log,
    :_pi => :pi,
    :_e => :(MathConstants.e)
)

function recursive_replace!(
        expr::Any,
        replace_requests::Dict{Symbol, Any} = REPLACEMENT_DICT
    )
    if expr isa Symbol
        haskey(replace_requests, expr) && return replace_requests[expr]
    elseif expr isa Expr
        for i in 1:length(expr.args)
            expr.args[i] = recursive_replace!(expr.args[i], replace_requests)
        end
    end
    return expr
end

"""
    preprocess_bng_expr(s::AbstractString)

Apply string-level preprocessing to convert BNG expression syntax to Julia
before calling `Meta.parse`. Currently handles:
- `if(cond, true, false)` → `ifelse(cond, true, false)`
"""
function preprocess_bng_expr(s::AbstractString)
    return replace(s, r"\bif\(" => "ifelse(")
end

"""
    parse_bng_expr(s::AbstractString)

Parse a BNG expression string into a Julia expression, applying BNG-to-Julia
preprocessing and symbol replacements.
"""
function parse_bng_expr(s::AbstractString)
    expr = Meta.parse(preprocess_bng_expr(s))
    (expr isa Expr) && recursive_replace!(expr)
    return expr
end

"""
    eval_bng_expr(opmod::Module, expr; extra_bindings=nothing)

Evaluate a parsed BNG expression in the given module context. Handles raw
`Expr`, `Symbol` (looked up in `extra_bindings` then `opmod` first, otherwise
created as a new parameter), and numeric literals.

`extra_bindings` is an optional `Dict{Symbol, Any}` of additional name→value
mappings that take precedence over `opmod`. This is needed because `Base.eval`
bindings created inside compiled functions are not visible to `isdefined` until
the next world age.
"""
function eval_bng_expr(opmod::Module, expr; extra_bindings = nothing)
    if expr isa Expr
        # substitute any extra_bindings into the expression before eval
        if extra_bindings !== nothing
            expr = _substitute_bindings(expr, extra_bindings)
        end
        return Base.eval(opmod, expr)
    elseif expr isa Symbol
        if extra_bindings !== nothing && haskey(extra_bindings, expr)
            return extra_bindings[expr]
        elseif isdefined(opmod, expr)
            return Base.eval(opmod, expr)
        else
            return (@parameters $expr)[1]
        end
    else
        return expr
    end
end

"""
Recursively substitute symbols in an expression using the given bindings dict.
"""
function _substitute_bindings(expr, bindings::Dict{Symbol, Any})
    if expr isa Symbol
        return get(bindings, expr, expr)
    elseif expr isa Expr
        new_args = map(a -> _substitute_bindings(a, bindings), expr.args)
        return Expr(expr.head, new_args...)
    else
        return expr
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Parameters
# ─────────────────────────────────────────────────────────────────────────────

const PARAM_BLOCK_START = "begin parameters"

function parse_params(ft::BNGNetwork, lines, boundaries)
    haskey(boundaries, PARAM_BLOCK_START) ||
        error("Required block '$PARAM_BLOCK_START' not found.")
    (start_idx, end_idx) = boundaries[PARAM_BLOCK_START]

    pvals = []
    ptoids = OrderedDict{Symbol, Int}()
    for idx in start_idx:(end_idx - 1)
        stripped = strip(lines[idx])
        isempty(stripped) && continue
        vals = split(stripped)
        pidx = parse(Int, vals[1])
        psym = Symbol(vals[2])
        ptoids[psym] = pidx

        pval = parse_bng_expr(vals[3])
        push!(pvals, pval)
    end

    return ptoids, pvals
end

# ─────────────────────────────────────────────────────────────────────────────
# Species
# ─────────────────────────────────────────────────────────────────────────────

const CHARS_TO_STRIP = ",~!().'`"
const STRIP_REGEX = Regex("[$CHARS_TO_STRIP]")
stripchars(s) = replace(s, STRIP_REGEX => "")

const SPECIES_BLOCK_START = "begin species"
const MAX_SYM_LEN = 8
const COMPARTMENT_PREFIX_RE = r"^@([A-Za-z_]\w*)::(.+)$"

"""
    extract_species_name(raw_name::AbstractString, sidx::Int)

Process a raw BNG species name, handling compartment prefixes (`@Comp::Name`)
and `\$`-prefixed fixed species. Returns `(shortsym, fullsym, is_fixed)`.
"""
function extract_species_name(raw_name::AbstractString, sidx::Int)
    name = raw_name
    is_fixed = false

    # strip $ prefix (fixed/constant species)
    if startswith(name, "\$")
        is_fixed = true
        name = name[2:end]
    end

    # handle compartment prefix @Comp::Name
    m = match(COMPARTMENT_PREFIX_RE, name)
    if m !== nothing
        comp = m.captures[1]
        species_part = m.captures[2]
        # check for $ after compartment prefix: @Comp::$Name
        if startswith(species_part, "\$")
            is_fixed = true
            species_part = species_part[2:end]
        end
        stripped = stripchars(species_part)
        shortsym = Symbol(comp * "_" * stripped)
        fullsym = Symbol(raw_name)
    else
        stripped = stripchars(name)
        if length(stripped) < MAX_SYM_LEN
            shortsym = Symbol(stripped)
        else
            shortsym = Symbol("S$sidx")
        end
        fullsym = Symbol(raw_name)
    end

    return shortsym, fullsym, is_fixed
end

function parse_species(ft::BNGNetwork, lines, boundaries)
    haskey(boundaries, SPECIES_BLOCK_START) ||
        error("Required block '$SPECIES_BLOCK_START' not found.")
    (start_idx, end_idx) = boundaries[SPECIES_BLOCK_START]

    shortsymstoids = OrderedDict{Symbol, Int}()
    shortsymstosyms = OrderedDict{Symbol, Symbol}()
    u0exprs = Vector{Any}()
    isfixed = BitVector()

    for idx in start_idx:(end_idx - 1)
        stripped = strip(lines[idx])
        isempty(stripped) && continue
        vals = split(stripped)
        sidx = parse(Int, vals[1])

        shortsym, fullsym, fixed = extract_species_name(vals[2], sidx)
        shortsymstoids[shortsym] = sidx
        shortsymstosyms[shortsym] = fullsym
        push!(u0exprs, parse_bng_expr(vals[3]))
        push!(isfixed, fixed)
    end

    return u0exprs, shortsymstoids, shortsymstosyms, isfixed
end

# ─────────────────────────────────────────────────────────────────────────────
# Groups (observables)
# ─────────────────────────────────────────────────────────────────────────────

const GROUPS_BLOCK_START = "begin groups"

function parse_groups(ft::BNGNetwork, lines, boundaries, shortsymstosyms,
                      idstovars, t)
    haskey(boundaries, GROUPS_BLOCK_START) ||
        error("Required block '$GROUPS_BLOCK_START' not found.")
    (start_idx, end_idx) = boundaries[GROUPS_BLOCK_START]

    namestosyms = Dict{String, Any}()
    obseqs = Equation[]
    syms_set = Set(keys(shortsymstosyms))

    for idx in start_idx:(end_idx - 1)
        stripped = strip(lines[idx])
        isempty(stripped) && continue
        vals = split(stripped)
        name = Symbol(vals[2])

        # skip groups whose name matches an existing species symbol or the time variable
        (name in syms_set || name == :t) && continue
        obs = (@variables $name($t))[1]
        namestosyms[vals[2]] = obs

        # build RHS as sum of species
        rhs = Num(0)
        groupterms = split(vals[3], ",")
        for groupterm in groupterms
            splitterm = split(groupterm, "*")
            if length(splitterm) == 1
                sid = parse(Int, splitterm[1])
                rhs += idstovars[sid]
            elseif length(splitterm) == 2
                coeff = parse(Int, splitterm[1])
                sid = parse(Int, splitterm[2])
                rhs += coeff * idstovars[sid]
            else
                error("Don't know how to handle term: $groupterm, appearing in a Group")
            end
        end
        push!(obseqs, Equation(obs, rhs))
    end

    return obseqs, namestosyms
end

# ─────────────────────────────────────────────────────────────────────────────
# Functions block
# ─────────────────────────────────────────────────────────────────────────────

const FUNCTIONS_BLOCK_START = "begin functions"

"""
    parse_functions!(ft, lines, boundaries, opmod, extra_bindings)

Parse the optional `begin functions` block. Returns `extra_bindings` (mutated)
with function name symbols mapped to their symbolic values.

Because `Base.eval` bindings created inside compiled functions are not visible
to `isdefined` until the next world age, we track function bindings in
`extra_bindings` (a `Dict{Symbol, Any}`) and pass them through to
`eval_bng_expr` for resolution.
"""
function parse_functions!(ft::BNGNetwork, lines, boundaries, opmod,
                          extra_bindings::Dict{Symbol, Any})
    haskey(boundaries, FUNCTIONS_BLOCK_START) || return extra_bindings
    (start_idx, end_idx) = boundaries[FUNCTIONS_BLOCK_START]

    # collect all function names first (for substitution in expressions)
    known_func_names = String[]
    for idx in start_idx:(end_idx - 1)
        stripped = strip(lines[idx])
        isempty(stripped) && continue
        vals = split(stripped)
        push!(known_func_names, replace(vals[2], "()" => ""))
    end

    for idx in start_idx:(end_idx - 1)
        stripped = strip(lines[idx])
        isempty(stripped) && continue
        vals = split(stripped)

        # vals[1] = index, vals[2] = "name()", vals[3:end] = expression tokens
        fname = replace(vals[2], "()" => "")

        # rejoin expression (may have spaces) and strip trailing comments
        fexpr_str = join(vals[3:end], " ")
        comment_pos = findfirst('#', fexpr_str)
        if comment_pos !== nothing
            fexpr_str = strip(fexpr_str[1:(comment_pos - 1)])
        end

        # replace known zero-arg function calls "funcname()" with "funcname"
        # so they resolve as variable lookups, not function calls
        for kf in known_func_names
            fexpr_str = replace(fexpr_str, Regex("\\b$(kf)\\(\\)") => kf)
        end

        fexpr = parse_bng_expr(fexpr_str)
        fval = eval_bng_expr(opmod, fexpr; extra_bindings)

        fsym = Symbol(fname)
        extra_bindings[fsym] = fval
    end

    return extra_bindings
end

# ─────────────────────────────────────────────────────────────────────────────
# Reactions
# ─────────────────────────────────────────────────────────────────────────────

const REACTIONS_BLOCK_START = "begin reactions"
const TYPED_RATE_KEYWORDS = Set(["Ele", "Sat", "MM", "Hill"])

function parse_reactions!(ft::BNGNetwork, t, lines, boundaries, idstovars, opmod;
                         extra_bindings = nothing)
    haskey(boundaries, REACTIONS_BLOCK_START) ||
        error("Required block '$REACTIONS_BLOCK_START' not found.")
    (start_idx, end_idx) = boundaries[REACTIONS_BLOCK_START]

    cntdict = Dict{Int, Int}()
    rxs = Reaction[]

    for idx in start_idx:(end_idx - 1)
        stripped = strip(lines[idx])
        isempty(stripped) && continue
        vals = split(stripped)
        reactantids = [parse(Int, rid) for rid in split(vals[2], ",")]
        productids = [parse(Int, pid) for pid in split(vals[3], ",")]

        # ── determine rate expression ──
        rate_token = vals[4]
        if rate_token in TYPED_RATE_KEYWORDS
            rate = _parse_typed_rate(rate_token, vals, reactantids,
                                     idstovars, cntdict, opmod;
                                     extra_bindings)
        else
            # ordinary elementary rate expression
            rateexpr = parse_bng_expr(rate_token)
            rate = eval_bng_expr(opmod, rateexpr; extra_bindings)
        end

        # ── validate null species (id 0 must be the sole entry if present) ──
        if any(iszero, reactantids)
            @assert length(reactantids) == 1 "Null species (id 0) must be the sole reactant, got: $reactantids"
        end
        if any(iszero, productids)
            @assert length(productids) == 1 "Null species (id 0) must be the sole product, got: $productids"
        end

        # ── reactant stoichiometry + factorial correction ──
        empty!(cntdict)
        for rid in reactantids
            (rid > 0) && (cntdict[rid] = get(cntdict, rid, 0) + 1)
        end
        scalefactor = isempty(cntdict) ? 0 : prod(factorial, values(cntdict))
        (!iszero(scalefactor)) && (rate = simplify(scalefactor * rate))

        rspecs = Vector{Any}(undef, length(cntdict))
        rstoich = Vector{Int}(undef, length(cntdict))
        for (i, (rid, cnt)) in enumerate(cntdict)
            rspecs[i] = idstovars[rid]
            rstoich[i] = cnt
        end

        # ── product stoichiometry ──
        empty!(cntdict)
        for pid in productids
            (pid > 0) && (cntdict[pid] = get(cntdict, pid, 0) + 1)
        end
        pspecs = Vector{Any}(undef, length(cntdict))
        pstoich = Vector{Int}(undef, length(cntdict))
        for (i, (pid, cnt)) in enumerate(cntdict)
            pspecs[i] = idstovars[pid]
            pstoich[i] = cnt
        end

        push!(rxs, Reaction(rate, rspecs, pspecs, rstoich, pstoich))
    end

    return rxs
end

"""
    _parse_typed_rate(rate_type, vals, reactantids, idstovars, cntdict, opmod)

Handle typed rate keywords (`Ele`, `Sat`, `MM`, `Hill`) in reaction lines.
Returns the rate expression to use (mass-action multiplication is still
applied by the caller for `Sat` and `Ele`).
"""
function _parse_typed_rate(rate_type, vals, reactantids, idstovars, cntdict, opmod;
                          extra_bindings = nothing)
    if rate_type == "Ele"
        # explicit elementary: use next token as rate expression
        rateexpr = parse_bng_expr(vals[5])
        return eval_bng_expr(opmod, rateexpr; extra_bindings)

    elseif rate_type == "Sat"
        # extract parameter tokens (everything between "Sat" and #comment)
        param_tokens = String[]
        for i in 5:length(vals)
            startswith(vals[i], "#") && break
            push!(param_tokens, vals[i])
        end

        # verify stoichiometry = 1 for all reactants
        empty!(cntdict)
        for rid in reactantids
            (rid > 0) && (cntdict[rid] = get(cntdict, rid, 0) + 1)
        end
        for (rid, cnt) in cntdict
            cnt > 1 && error(
                "Sat rate law with reactant stoichiometry > 1 is not supported. " *
                "BNG's stochastic Sat uses binomial coefficients in the saturation " *
                "term, which cannot be represented as a single symbolic expression " *
                "valid for both ODE and stochastic simulation. Consider reformulating " *
                "as a function block in the .bngl source."
            )
        end

        # resolve parameters
        params = [eval_bng_expr(opmod, parse_bng_expr(tok); extra_bindings) for tok in param_tokens]
        kcat = params[1]
        Km_vals = params[2:end]

        # build rate constant: kcat / prod(Kmi + Ri) for substrates
        # substrates are the first len(Km) reactants; rest are non-substrates
        # collect unique reactant ids in order of first appearance
        unique_rids = Int[]
        for rid in reactantids
            (rid > 0 && !(rid in unique_rids)) && push!(unique_rids, rid)
        end

        rate = kcat
        for (i, Km) in enumerate(Km_vals)
            i <= length(unique_rids) || break
            rate = rate / (Km + idstovars[unique_rids[i]])
        end
        return rate

    elseif rate_type in ("MM", "Hill")
        error(
            "Typed rate law '$rate_type' is not currently supported for .net import. " *
            "MM uses a tQSSA formula that cannot be decomposed into rate constant + " *
            "mass action, and Hill has binomial coefficient issues for stoichiometry " *
            "> 1. If you need support for this rate type, please open an issue at " *
            "the ReactionNetworkImporters.jl repository."
        )
    else
        error("Unknown typed rate keyword: '$rate_type'")
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Parameter / initial condition evaluation
# ─────────────────────────────────────────────────────────────────────────────

function exprs_to_defs(opmod, ptoids, pvals, idstovars, u0exprs, ps, isfixed;
                       extra_bindings = nothing)
    pmap = Dict()
    psym_to_pvar = Dict(nameof(p) => p for p in ps)
    for (psym, pid) in ptoids
        pvar = psym_to_pvar[psym]
        parsedval = pvals[pid]
        pval = eval_bng_expr(opmod, parsedval; extra_bindings)
        push!(pmap, pvar => pval)
    end

    u0map = Dict()
    for (i, u0expr) in enumerate(u0exprs)
        uvar = idstovars[i]
        u0val = eval_bng_expr(opmod, u0expr; extra_bindings)
        if isfixed[i]
            push!(pmap, uvar => u0val)
        else
            push!(u0map, uvar => u0val)
        end
    end

    return union(pmap, u0map), pmap, u0map
end

# ─────────────────────────────────────────────────────────────────────────────
# Main entry point
# ─────────────────────────────────────────────────────────────────────────────

"""
    loadrxnetwork(ft::BNGNetwork, rxfilename; name = gensym(:ReactionSystem), verbose = true, kwargs...)

Parse a BioNetGen `.net` file and construct a Catalyst `ReactionSystem`.

# Arguments
- `ft::BNGNetwork`: Indicates the file to be parsed is a BioNetGen ".net" file.
- `rxfilename::String`: Path to the `.net` file to be parsed.
- `name::Symbol`: (Optional) Name for the resulting `ReactionSystem`. Defaults to a
  generated symbol.
- `verbose::Bool`: (Optional) If `true`, prints detailed progress information during
  parsing. Defaults to `true`.
- `kwargs...`: Additional keyword arguments passed to the `ReactionSystem` constructor.

# Returns
A Catalyst `ReactionSystem` (not marked as complete). Initial conditions, parameter
values, and BNG-specific mappings are stored as system metadata:
- `Catalyst.get_u0_map(rn)`: `Dict` mapping species to initial condition values.
- `Catalyst.get_parameter_map(rn)`: `Dict` mapping parameters to their values.
- `get_varstonames(rn)`: `Dict` mapping internal symbolic variables (both dynamic
  species and constant-species parameters) to full BNG name strings.
- `get_groupstosyms(rn)`: `Dict` mapping BNG group name strings to observable symbols.

# Supported features
- Parameter expressions (including `min`, `max`, `if`/`ifelse`, `mratio`)
- Species with initial conditions (including `\$`-prefixed constant species)
- Compartment-prefixed species (`@Comp::Name` → `Comp_Name` symbol)
- Reactions with elementary and `Sat` typed rate laws
- Function blocks (concentration-dependent and time-dependent rate expressions)
- Groups (observables as weighted sums of species)

# Unsupported features (will error)
- `MM` and `Hill` typed rate laws
- `tfun()` tabular interpolation functions

# Example
```julia
rn = loadrxnetwork(BNGNetwork(), "path/to/network.net"; verbose = true)
rn = complete(rn)
```
"""
function loadrxnetwork(
        ft::BNGNetwork, rxfilename; name = gensym(:ReactionSystem),
        verbose = true, kwargs...
    )
    lines = readlines(rxfilename)
    t = Catalyst.default_t()

    # ── Step 1: discover all blocks ──
    verbose && print("Scanning blocks...")
    boundaries = find_block_boundaries(lines)
    verbose && println("done")

    # ── Step 2: parse parameters ──
    verbose && print("Parsing parameters...")
    ptoids, pvals = parse_params(ft, lines, boundaries)
    verbose && println("done")

    verbose && print("Creating parameters...")
    ps = [(@parameters $psym)[1] for psym in keys(ptoids)]
    verbose && println("done")

    # ── Step 3: parse species (with fixed-species + compartment detection) ──
    verbose && print("Parsing species...")
    u0exprs, shortsymstoids, shortsymstosyms, isfixed = parse_species(
        ft, lines, boundaries)
    verbose && println("done")

    # ── Step 4: build unified id-to-variable map ──
    verbose && print("Creating variables...")
    n_species = length(shortsymstoids)
    idstovars = Vector{Any}(undef, n_species)
    idstoshortsyms = Vector{Symbol}(undef, n_species)
    dynamic_specs = Any[]
    constant_specs = Any[]

    for (shortsym, sid) in shortsymstoids
        idstoshortsyms[sid] = shortsym
        if isfixed[sid]
            var = only(@parameters $shortsym [isconstantspecies = true])
            idstovars[sid] = var
            push!(constant_specs, var)
        else
            var = funcsym(shortsym, t)
            idstovars[sid] = var
            push!(dynamic_specs, var)
        end
    end
    verbose && println("done")

    # ── Step 5: set up opmod with parameters, variables, time, and built-ins ──
    verbose && print("Setting up expression evaluation module...")
    opmod = Module()
    Base.eval(opmod, :(using Catalyst))
    Base.eval(opmod, :(t = Catalyst.default_t()))
    Base.eval(opmod, :(time() = Catalyst.default_t()))

    # register mratio using HypergeometricFunctions
    Base.eval(opmod, quote
        using HypergeometricFunctions: _₁F₁
        mratio(a, b, z) = _₁F₁(a + 1, b + 1, z) / _₁F₁(a, b, z)
    end)

    for p in ps
        psym = nameof(p)
        Base.eval(opmod, :($psym = $p))
    end
    for (sid, var) in enumerate(idstovars)
        if isfixed[sid]
            ssym = nameof(var)
        else
            ssym = nameof(operation(unwrap(var)))
        end
        Base.eval(opmod, :($ssym = $var))
    end
    verbose && println("done")

    # ── Step 6: parse groups (before functions — functions reference groups) ──
    verbose && print("Parsing groups...")
    obseqs, groupstosyms = parse_groups(
        ft, lines, boundaries, shortsymstosyms, idstovars, t)
    verbose && println("done")

    # Build extra_bindings for group RHS expressions and function values.
    # We use a Dict rather than Base.eval because bindings created via Base.eval
    # inside compiled functions are not visible to isdefined until the next world age.
    extra_bindings = Dict{Symbol, Any}()
    for obseq in obseqs
        gsym = nameof(operation(unwrap(obseq.lhs)))
        extra_bindings[gsym] = obseq.rhs
    end

    # ── Step 7: parse functions (optional block) ──
    verbose && print("Parsing functions...")
    parse_functions!(ft, lines, boundaries, opmod, extra_bindings)
    verbose && println("done")

    # ── Step 8: parse reactions ──
    verbose && print("Parsing and adding reactions...")
    rxs = parse_reactions!(ft, t, lines, boundaries, idstovars, opmod;
                           extra_bindings)
    verbose && println("done")

    # ── Step 9: build ReactionSystem ──
    # evaluate parameter values and initial conditions
    all_ps = vcat(ps, constant_specs)
    defmap, pmap, u0map = exprs_to_defs(
        opmod, ptoids, pvals, idstovars, u0exprs, all_ps, isfixed;
        extra_bindings)

    # build metadata
    metadata = [Catalyst.U0Map => u0map, Catalyst.ParameterMap => pmap,
        VarsToNames => shortsymstosyms, GroupsToSyms => groupstosyms]
    rn = ReactionSystem(
        rxs, t, dynamic_specs, all_ps; name, observed = obseqs,
        initial_conditions = defmap,
        metadata, kwargs...
    )

    return rn
end
