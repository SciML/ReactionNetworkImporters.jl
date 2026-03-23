using Catalyst, OrdinaryDiffEqTsit5, DataFrames, CSV
using ReactionNetworkImporters

# ─────────────────────────────────────────────────────────────────────────────
# Helper: load BNG .gdat file (handles # prefix on header)
# ─────────────────────────────────────────────────────────────────────────────
function load_gdat(path)
    lines = readlines(path)
    if startswith(lines[1], "#")
        lines[1] = lstrip(lines[1], ['#', ' '])
    end
    CSV.read(IOBuffer(join(lines, '\n')), DataFrame; delim = ' ', ignorerepeated = true)
end

"""
Load a .net file, solve the ODE, and compare against .gdat reference data
for a given observable column. Asserts relative error < atol.
"""
function test_ode_match(netfile, gdatfile, obs_col::Symbol;
                        tf, nsteps, atol = 1e-4, maxiters = 100_000)
    rn = loadrxnetwork(BNGNetwork(), netfile; verbose = false)
    rn_c = complete(rn)
    gdatdf = load_gdat(gdatfile)

    prob = ODEProblem(rn_c, Float64[], (0.0, tf), Float64[])
    sol = solve(prob, Tsit5(); abstol = 1e-12, reltol = 1e-12,
                saveat = tf / nsteps, maxiters)

    @test sol.retcode == ReturnCode.Success
    @test length(sol.t) == nrow(gdatdf)

    bng_vals = gdatdf[!, obs_col]
    cat_vals = sol[getproperty(rn_c, obs_col)]
    for (i, (bng, cat)) in enumerate(zip(bng_vals, cat_vals))
        if abs(bng) > 1e-10
            @test abs(bng - cat) < atol * abs(bng)
        else
            @test abs(bng - cat) < 1e-10
        end
    end
end

# ═════════════════════════════════════════════════════════════════════════════
# toy-jim.net — group-parameter name collision regression + higher-order
# ═════════════════════════════════════════════════════════════════════════════
@testset "toy-jim.net (group-param collision + ODE)" begin
    datadir = joinpath(@__DIR__, "../data/toy-jim")
    test_ode_match(
        joinpath(datadir, "toy-jim.net"),
        joinpath(datadir, "toy-jim.gdat"),
        :RecDim; tf = 120.0, nsteps = 120, atol = 1e-3
    )
end

# ═════════════════════════════════════════════════════════════════════════════
# localfunc.net — local function-style rates + fixed species
# ═════════════════════════════════════════════════════════════════════════════
@testset "localfunc.net (local functions + fixed species)" begin
    datadir = joinpath(@__DIR__, "../data/localfunc")
    test_ode_match(
        joinpath(datadir, "localfunc.net"),
        joinpath(datadir, "localfunc.gdat"),
        :Atot; tf = 10.0, nsteps = 40, atol = 1e-4
    )
end

# ═════════════════════════════════════════════════════════════════════════════
# test_synthesis_simple.net — null synthesis (0 -> X)
# ═════════════════════════════════════════════════════════════════════════════
@testset "test_synthesis_simple.net (null synthesis)" begin
    datadir = joinpath(@__DIR__, "../data/test_synthesis_simple")
    test_ode_match(
        joinpath(datadir, "test_synthesis_simple.net"),
        joinpath(datadir, "test_synthesis_simple.gdat"),
        :add_molecule; tf = 40.0, nsteps = 100, atol = 1e-6
    )
end

# ═════════════════════════════════════════════════════════════════════════════
# test_synthesis_cBNGL_simple.net — compartment + constant source
# ═════════════════════════════════════════════════════════════════════════════
@testset "test_synthesis_cBNGL_simple.net (compartment + source)" begin
    datadir = joinpath(@__DIR__, "../data/test_synthesis_cBNGL_simple")
    test_ode_match(
        joinpath(datadir, "test_synthesis_cBNGL_simple.net"),
        joinpath(datadir, "test_synthesis_cBNGL_simple.gdat"),
        :compartment_suffix; tf = 40.0, nsteps = 100, atol = 1e-6
    )
end

# ═════════════════════════════════════════════════════════════════════════════
# test_synthesis_complex_source_cBNGL.net — compartment + complex source
# ═════════════════════════════════════════════════════════════════════════════
@testset "test_synthesis_complex_source_cBNGL.net" begin
    datadir = joinpath(@__DIR__, "../data/test_synthesis_complex_source_cBNGL")
    test_ode_match(
        joinpath(datadir, "test_synthesis_complex_source_cBNGL.net"),
        joinpath(datadir, "test_synthesis_complex_source_cBNGL.gdat"),
        :vs_suffix; tf = 40.0, nsteps = 100, atol = 1e-6
    )
end

# ═════════════════════════════════════════════════════════════════════════════
# test_synthesis_complex.net — complex synthesis patterns
# ═════════════════════════════════════════════════════════════════════════════
@testset "test_synthesis_complex.net" begin
    datadir = joinpath(@__DIR__, "../data/test_synthesis_complex")
    test_ode_match(
        joinpath(datadir, "test_synthesis_complex.net"),
        joinpath(datadir, "test_synthesis_complex.gdat"),
        :Receptor; tf = 40.0, nsteps = 100, atol = 1e-6
    )
end

# ═════════════════════════════════════════════════════════════════════════════
# test_ANG_synthesis_simple.net — ANG synthesis
# ═════════════════════════════════════════════════════════════════════════════
@testset "test_ANG_synthesis_simple.net" begin
    datadir = joinpath(@__DIR__, "../data/test_ANG_synthesis_simple")
    test_ode_match(
        joinpath(datadir, "test_ANG_synthesis_simple.net"),
        joinpath(datadir, "test_ANG_synthesis_simple.gdat"),
        :add_molecule; tf = 40.0, nsteps = 100, atol = 1e-6
    )
end

# ═════════════════════════════════════════════════════════════════════════════
# univ_synth.net — compartment-aware null synthesis/degradation
# ═════════════════════════════════════════════════════════════════════════════
@testset "univ_synth.net (compartment synthesis/degradation)" begin
    datadir = joinpath(@__DIR__, "../data/univ_synth")
    test_ode_match(
        joinpath(datadir, "univ_synth.net"),
        joinpath(datadir, "univ_synth.gdat"),
        :B_CP; tf = 10.0, nsteps = 20, atol = 1e-4
    )
end

# ═════════════════════════════════════════════════════════════════════════════
# test_sbml_structured.net — compartment + functions
# ═════════════════════════════════════════════════════════════════════════════
@testset "test_sbml_structured.net (compartment + functions)" begin
    datadir = joinpath(@__DIR__, "../data/test_sbml_structured")
    test_ode_match(
        joinpath(datadir, "test_sbml_structured.net"),
        joinpath(datadir, "test_sbml_structured.gdat"),
        :MolA_cell; tf = 100.0, nsteps = 100, atol = 1e-4
    )
end

# ═════════════════════════════════════════════════════════════════════════════
# tlbr.net — ligand-receptor binding network
# ═════════════════════════════════════════════════════════════════════════════
@testset "tlbr.net (ligand-receptor network)" begin
    datadir = joinpath(@__DIR__, "../data/tlbr")
    test_ode_match(
        joinpath(datadir, "tlbr.net"),
        joinpath(datadir, "tlbr.gdat"),
        :Ltot; tf = 2.0, nsteps = 50, atol = 1e-4
    )
end

# ═════════════════════════════════════════════════════════════════════════════
# Haugh2b.net — medium-sized signaling network
# ═════════════════════════════════════════════════════════════════════════════
@testset "Haugh2b.net (signaling network)" begin
    datadir = joinpath(@__DIR__, "../data/Haugh2b")
    test_ode_match(
        joinpath(datadir, "Haugh2b.net"),
        joinpath(datadir, "Haugh2b.gdat"),
        :S2_P_tot; tf = 50.0, nsteps = 5, atol = 1e-6
    )
end
