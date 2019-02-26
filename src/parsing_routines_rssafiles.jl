function parse_species(ft::RSSANetwork, fname)
    specs_ic = Dict{Symbol,Int}()

    open(fname, "r") do file
        for line in eachline(file)
            var,val = split(line, "=")            
            specs_ic[Symbol(strip(var))] = parse(Int64, val)
        end
    end

    specs_ic
end

function parse_reactions(ft::RSSANetwork, fname)
    rxstrs  = Vector{String}()
    rxrates = Vector{Float64}()

    open(fname, "r") do file
        for line in eachline(file)
            rxstr,rate = split(line, ",")
            rxstr = replace(strip(rxstr), "-" => "--")
            rxstr = replace(strip(rxstr), "_" => "0")
            push!(rxstrs, strip(rxstr))
            push!(rxrates, parse(Float64, rate))
        end
    end

    rxstrs, rxrates
end

function build_rxnetwork(ft::RSSANetwork, networkname, rxstrs, rxrates; printrxs=false, kwargs...)
    
    # string representing the network
    rxiobuf = IOBuffer()
    write(rxiobuf, "@min_reaction_network $(networkname) begin\n")
    for (i,rxstr) in enumerate(rxstrs)
        rate = rxrates[i]
        write(rxiobuf, "\t $rate, $rxstr \n")
    end
    write(rxiobuf, "end\n")
    rnstr = String(take!(rxiobuf))

    if printrxs
        print(rnstr)
    end

    # build the network using DiffEqBiological
    rn = eval( Meta.parse(rnstr) )

    rn,rnstr
end

function get_init_condit(ft::RSSANetwork, rn, specs_ic)
    icvec = zeros(Int64, length(specs_ic))    
    for (i,sym) in enumerate(rn.syms)
        icvec[i] = specs_ic[sym]
    end

    icvec
end


# for parsing the simple format from the book by Thanh et al.
function loadrxnetwork(ft::RSSANetwork, networkname, specs_ic_file, rxs_file; kwargs...)

    # parse initial conditions
    specs_ic = parse_species(ft, specs_ic_file)

    # parse reaction network
    rxstrs,rxrates = parse_reactions(ft, rxs_file)

    # build the DiffEqBiological representation of the network
    rn,rnstr = build_rxnetwork(ft, networkname, rxstrs, rxrates; kwargs...)
    initialpop = get_init_condit(ft, rn, specs_ic)

    ParsedReactionNetwork(rn, initialpop; rnstr=rnstr)
end