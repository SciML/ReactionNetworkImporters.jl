using Catalyst, OrdinaryDiffEqTsit5, DataFrames, CSV, LinearAlgebra
using ReactionNetworkImporters
using ReactionNetworkImporters: find_block_boundaries

# ─────────────────────────────────────────────────────────────────────────────
# Helper: load .gdat and compare ODE solution against BNG reference data
# ─────────────────────────────────────────────────────────────────────────────
function load_gdat(path)
    # BNG .gdat files may have a '#' prefix on the header line; strip it
    lines = readlines(path)
    if startswith(lines[1], "#")
        lines[1] = lstrip(lines[1], ['#', ' '])
    end
    CSV.read(IOBuffer(join(lines, '\n')), DataFrame; delim = ' ', ignorerepeated = true)
end

function test_ode_vs_gdat(rn_incomplete, gdatfile, obs_col::Symbol;
                          tf, nsteps, atol = 1e-6, rtol = 1e-12)
    gdatdf = load_gdat(gdatfile)
    rn = complete(rn_incomplete)
    prob = ODEProblem(rn, Float64[], (0.0, tf), Float64[])
    sol = solve(prob, Tsit5(); abstol = 1e-12, reltol = rtol,
                saveat = tf / nsteps)
    @test all(sol.t .== gdatdf[!, :time])
    bng_vals = gdatdf[!, obs_col]
    # use relative tolerance where values are non-tiny, absolute otherwise
    cat_vals = sol[getproperty(rn, obs_col)]
    for (i, (bng, cat)) in enumerate(zip(bng_vals, cat_vals))
        if abs(bng) > 1e-10
            @test abs(bng - cat) < atol * abs(bng)
        else
            @test abs(bng - cat) < 1e-10
        end
    end
end

# ═════════════════════════════════════════════════════════════════════════════
# Test: block boundary discovery
# ═════════════════════════════════════════════════════════════════════════════
@testset "Block boundary discovery" begin
    lines = [
        "# comment",
        "substanceUnits(\"Number\");",
        "begin parameters",
        "  1 k 0.5",
        "end parameters",
        "begin species",
        "  1 A() 100",
        "end species",
    ]
    boundaries = find_block_boundaries(lines)
    @test haskey(boundaries, "begin parameters")
    @test haskey(boundaries, "begin species")
    @test !haskey(boundaries, "begin functions")
    # content lines are between begin+1 and end-1
    pstart, pend = boundaries["begin parameters"]
    @test pstart == 4
    @test pend == 5
end

# ═════════════════════════════════════════════════════════════════════════════
# Test: michment.net — function block referencing observables
# ═════════════════════════════════════════════════════════════════════════════
@testset "michment.net (function block)" begin
    datadir = joinpath(@__DIR__, "../data/michment")
    fname = joinpath(datadir, "michment.net")

    rn = loadrxnetwork(BNGNetwork(), fname; verbose = false)
    @test Catalyst.has_u0_map(rn)
    @test Catalyst.has_parameter_map(rn)
    @test has_varstonames(rn)
    @test has_groupstosyms(rn)

    @test length(species(rn)) == 5
    @test length(reactions(rn)) == 6

    # verify function-based rate resolves to symbolic expression, not a bare parameter
    rx_rates = [rx.rate for rx in reactions(rn)]
    # michment() = kcat/(Km+Sa0) should appear as a division expression
    @test any(r -> !Symbolics.issym(Symbolics.unwrap(r)), rx_rates)

    # verify ODE can be built and solved (gdat ICs don't match .net ICs so no comparison)
    rn_c = complete(rn)
    prob = ODEProblem(rn_c, Float64[], (0.0, 1.0), Float64[])
    sol = solve(prob, Tsit5(); abstol = 1e-12, reltol = 1e-12)
    @test sol.retcode == ReturnCode.Success
end

# ═════════════════════════════════════════════════════════════════════════════
# Test: CaOscillate_Func.net — multiple functions referencing observables
# ═════════════════════════════════════════════════════════════════════════════
@testset "CaOscillate_Func.net (function block)" begin
    datadir = joinpath(@__DIR__, "../data/CaOscillate_Func")
    fname = joinpath(datadir, "CaOscillate_Func.net")
    gdatfile = joinpath(datadir, "CaOscillate_Func.gdat")

    rn = loadrxnetwork(BNGNetwork(), fname; verbose = false)
    @test length(species(rn)) == 4
    @test length(reactions(rn)) == 8

    test_ode_vs_gdat(rn, gdatfile, :Ga; tf = 50.0, nsteps = 500, atol = 1e-4)
end

# ═════════════════════════════════════════════════════════════════════════════
# Test: CaOscillate_Sat.net — Sat typed rate law
# ═════════════════════════════════════════════════════════════════════════════
@testset "CaOscillate_Sat.net (Sat typed rate)" begin
    datadir = joinpath(@__DIR__, "../data/CaOscillate_Sat")
    fname = joinpath(datadir, "CaOscillate_Sat.net")
    gdatfile = joinpath(datadir, "CaOscillate_Sat.gdat")

    rn = loadrxnetwork(BNGNetwork(), fname; verbose = false)
    @test length(species(rn)) == 4
    @test length(reactions(rn)) == 8

    # Sat and Func versions of CaOscillate should produce identical results
    test_ode_vs_gdat(rn, gdatfile, :Ga; tf = 50.0, nsteps = 500, atol = 1e-4)
end

# ═════════════════════════════════════════════════════════════════════════════
# Test: motor.net — complex exponential functions
# ═════════════════════════════════════════════════════════════════════════════
@testset "motor.net (complex function expressions)" begin
    datadir = joinpath(@__DIR__, "../data/motor")
    fname = joinpath(datadir, "motor.net")
    gdatfile = joinpath(datadir, "motor.gdat")

    rn = loadrxnetwork(BNGNetwork(), fname; verbose = false)
    @test length(species(rn)) == 4
    @test length(reactions(rn)) == 4

    test_ode_vs_gdat(rn, gdatfile, :CheYp; tf = 0.2, nsteps = 20, atol = 1e-6)
end

# ═════════════════════════════════════════════════════════════════════════════
# Test: test_time.net — time-dependent functions
# ═════════════════════════════════════════════════════════════════════════════
@testset "test_time.net (time-dependent functions)" begin
    datadir = joinpath(@__DIR__, "../data/test_time")
    fname = joinpath(datadir, "test_time.net")
    gdatfile = joinpath(datadir, "test_time.gdat")

    rn = loadrxnetwork(BNGNetwork(), fname; verbose = false)
    @test length(species(rn)) == 3
    @test length(reactions(rn)) == 3

    test_ode_vs_gdat(rn, gdatfile, :x_exact; tf = 10.0, nsteps = 100, atol = 1e-4)
end

# ═════════════════════════════════════════════════════════════════════════════
# Test: test_fixed.net — $-prefixed constant species
# ═════════════════════════════════════════════════════════════════════════════
@testset "test_fixed.net (constant species)" begin
    datadir = joinpath(@__DIR__, "../data/test_fixed")
    fname = joinpath(datadir, "test_fixed.net")

    rn = loadrxnetwork(BNGNetwork(), fname; verbose = false)

    # $A, $Degr, $Trash are constant species (parameters), rest are dynamic
    @test length(species(rn)) == 9   # 12 total - 3 constant = 9 dynamic

    # constant species should appear in parameter map
    pmap = Catalyst.get_parameter_map(rn)
    pmap_names = Set(nameof(k) for (k, _) in pmap)
    @test :A in pmap_names
    @test :Degr in pmap_names
    @test :Trash in pmap_names

    # constant species should appear in VarsToNames
    vtn = get_varstonames(rn)
    vtn_shortsyms = Set(keys(vtn))
    @test :A in vtn_shortsyms
    @test :Degr in vtn_shortsyms
    @test :Trash in vtn_shortsyms

    # verify we can build and solve an ODE
    rn_c = complete(rn)
    prob = ODEProblem(rn_c, Float64[], (0.0, 60.0), Float64[])
    sol = solve(prob, Tsit5(); abstol = 1e-12, reltol = 1e-12)
    @test sol.retcode == ReturnCode.Success
end

# ═════════════════════════════════════════════════════════════════════════════
# Test: test_assignment.net — nested if() in parameter expressions
# ═════════════════════════════════════════════════════════════════════════════
@testset "test_assignment.net (if/ifelse in parameters)" begin
    datadir = joinpath(@__DIR__, "../data/test_assignment")
    fname = joinpath(datadir, "test_assignment.net")

    rn = loadrxnetwork(BNGNetwork(), fname; verbose = false)
    @test length(species(rn)) == 4
    # just verify it loads without error — the if() preprocessing works
end

# ═════════════════════════════════════════════════════════════════════════════
# Test: test_mratio.net — mratio + min/max in parameter expressions
# ═════════════════════════════════════════════════════════════════════════════
@testset "test_mratio.net (mratio/min/max builtins)" begin
    datadir = joinpath(@__DIR__, "../data/test_mratio")
    fname = joinpath(datadir, "test_mratio.net")
    gdatfile = joinpath(datadir, "test_mratio.gdat")

    rn = loadrxnetwork(BNGNetwork(), fname; verbose = false)
    @test length(species(rn)) == 3
    @test length(reactions(rn)) == 2

    test_ode_vs_gdat(rn, gdatfile, :C_obs; tf = 10.0, nsteps = 1000, atol = 1e-4)
end

# ═════════════════════════════════════════════════════════════════════════════
# Test: test_sbml_flat.net — compartment-prefixed species
# ═════════════════════════════════════════════════════════════════════════════
@testset "test_sbml_flat.net (compartment-prefixed species)" begin
    datadir = joinpath(@__DIR__, "../data/test_sbml_flat")
    fname = joinpath(datadir, "test_sbml_flat.net")

    rn = loadrxnetwork(BNGNetwork(), fname; verbose = false)
    @test length(species(rn)) == 5
    @test length(reactions(rn)) == 6

    # verify compartment-aware naming: @cell::A____() should produce cell_A____
    vtn = get_varstonames(rn)
    @test vtn !== nothing
    short_names = Set(String(k) for k in keys(vtn))
    @test "cell_A____" in short_names
end

# ═════════════════════════════════════════════════════════════════════════════
# Test: MM/Hill error handling
# ═════════════════════════════════════════════════════════════════════════════
@testset "Species-block expressions with BNG builtins" begin
    # Test that species ICs using BNG builtins (if, parameter refs) parse correctly
    net_str = """
    begin parameters
        1 k1 1.0
        2 N0 50
    end parameters
    begin species
        1 A() if((1==1),N0,0)
        2 B() N0
    end species
    begin reactions
        1 1 2 k1 #Rule1
    end reactions
    begin groups
        1 A_tot 1
        2 B_tot 2
    end groups
    """
    tmpfile = tempname() * ".net"
    write(tmpfile, net_str)
    rn = loadrxnetwork(BNGNetwork(), tmpfile; verbose = false)
    u0map = Catalyst.get_u0_map(rn)
    # if((1==1), N0, 0) evaluates to the symbolic parameter N0 (which has value 50)
    u0_vals = Dict(nameof(Symbolics.operation(Symbolics.unwrap(k))) => v
                   for (k, v) in u0map)
    @test isequal(u0_vals[:A], u0_vals[:B])  # both should be N0

    # verify the ODE builds and solves correctly (parameter substitution works)
    rn_c = complete(rn)
    prob = ODEProblem(rn_c, Float64[], (0.0, 1.0), Float64[])
    sol = solve(prob, Tsit5(); abstol = 1e-12, reltol = 1e-12)
    @test sol.retcode == ReturnCode.Success
    # A(0) should be 50 (= N0)
    @test sol[species(rn_c)[1]][1] ≈ 50.0
    rm(tmpfile)
end

@testset "Typed rate errors and edge cases" begin
    # MM error
    mm_net = """
    begin parameters
        1 kcat 1.0
        2 Km 0.5
    end parameters
    begin species
        1 S() 100
        2 E() 10
    end species
    begin reactions
        1 1,2 2 MM kcat Km #Rule1
    end reactions
    begin groups
        1 S_tot 1
    end groups
    """
    tmpfile = tempname() * ".net"
    write(tmpfile, mm_net)
    @test_throws ErrorException loadrxnetwork(BNGNetwork(), tmpfile; verbose = false)
    rm(tmpfile)

    # Hill error
    hill_net = """
    begin parameters
        1 Vmax 1.0
        2 Kh 0.5
        3 h 2
    end parameters
    begin species
        1 S() 100
        2 E() 10
    end species
    begin reactions
        1 1,2 2 Hill Vmax Kh h #Rule1
    end reactions
    begin groups
        1 S_tot 1
    end groups
    """
    tmpfile = tempname() * ".net"
    write(tmpfile, hill_net)
    @test_throws ErrorException loadrxnetwork(BNGNetwork(), tmpfile; verbose = false)
    rm(tmpfile)

    # Sat with stoich > 1 error
    sat_high_net = """
    begin parameters
        1 kcat 1.0
        2 Km 0.5
    end parameters
    begin species
        1 A() 100
        2 B() 0
    end species
    begin reactions
        1 1,1 2 Sat kcat Km #Rule1
    end reactions
    begin groups
        1 A_tot 1
    end groups
    """
    tmpfile = tempname() * ".net"
    write(tmpfile, sat_high_net)
    @test_throws ErrorException loadrxnetwork(BNGNetwork(), tmpfile; verbose = false)
    rm(tmpfile)

    # Ele typed rate (should work, just skip the keyword)
    ele_net = """
    begin parameters
        1 k1 0.5
    end parameters
    begin species
        1 A() 100
        2 B() 0
    end species
    begin reactions
        1 1 2 Ele k1 #Rule1
    end reactions
    begin groups
        1 A_tot 1
    end groups
    """
    tmpfile = tempname() * ".net"
    write(tmpfile, ele_net)
    rn = loadrxnetwork(BNGNetwork(), tmpfile; verbose = false)
    @test length(reactions(rn)) == 1
    # rate should be the parameter k1, not a new parameter Ele
    rx_rate_name = nameof(Symbolics.unwrap(reactions(rn)[1].rate))
    @test rx_rate_name == :k1
    rm(tmpfile)
end

# ═════════════════════════════════════════════════════════════════════════════
# Test: .net files without groups block
# ═════════════════════════════════════════════════════════════════════════════
@testset "No groups block (blbr.net)" begin
    datadir = joinpath(@__DIR__, "../data/blbr")
    fname = joinpath(datadir, "blbr.net")

    rn = loadrxnetwork(BNGNetwork(), fname; verbose = false)
    @test length(species(rn)) > 0
    @test length(reactions(rn)) > 0
    # no groups → no observables, but should still load
    @test isempty(observed(rn))

    # verify ODE can be built
    rn_c = complete(rn)
    prob = ODEProblem(rn_c, Float64[], (0.0, 1.0), Float64[])
    sol = solve(prob, Tsit5(); abstol = 1e-12, reltol = 1e-12)
    @test sol.retcode == ReturnCode.Success
end

# ═════════════════════════════════════════════════════════════════════════════
# Test: group name colliding with parameter name (toy-jim.net)
# ═════════════════════════════════════════════════════════════════════════════
@testset "Group-parameter name collision (toy-jim.net)" begin
    datadir = joinpath(@__DIR__, "../data/toy-jim")
    fname = joinpath(datadir, "toy-jim.net")

    rn = loadrxnetwork(BNGNetwork(), fname; verbose = false)
    @test length(species(rn)) > 0
    @test length(reactions(rn)) > 0

    # verify the model constructs without "duplicate names" error
    rn_c = complete(rn)
    prob = ODEProblem(rn_c, Float64[], (0.0, 1.0), Float64[])
    sol = solve(prob, Tsit5(); abstol = 1e-12, reltol = 1e-12)
    @test sol.retcode == ReturnCode.Success
end
