using ReactionNetworkImporters, JET, Test
using Catalyst, SparseArrays

@testset "JET static analysis" begin
    # Test type stability of internal parsing functions
    @testset "Internal parsing functions" begin
        # Test stripchars returns String
        @test ReactionNetworkImporters.stripchars("ABC()") isa String
        @test ReactionNetworkImporters.stripchars("ABC()") == "ABC"

        # Test hasstripchars returns Bool
        @test ReactionNetworkImporters.hasstripchars("ABC()") isa Bool
        @test ReactionNetworkImporters.hasstripchars("ABC()") == true
        @test ReactionNetworkImporters.hasstripchars("ABC") == false

        # Test make_shortsym returns Symbol
        @test ReactionNetworkImporters.make_shortsym("ABC", 1) isa Symbol
        @test ReactionNetworkImporters.make_shortsym("VeryLongSpeciesName", 1) == :S1

        # Test seek_to_block returns Int
        lines = ["begin parameters", "end parameters"]
        @test ReactionNetworkImporters.seek_to_block(lines, 1, "begin parameters") isa Int
    end

    @testset "JET report_call on key functions" begin
        # Test that key internal functions don't have obvious dispatch issues
        # Note: We use target_modules to focus on ReactionNetworkImporters
        # and ignore issues from upstream dependencies (SymbolicUtils, etc.)

        # Test stripchars
        rep = JET.report_call(
            ReactionNetworkImporters.stripchars,
            (String,);
            target_modules = (ReactionNetworkImporters,)
        )
        @test length(JET.get_reports(rep)) == 0

        # Test hasstripchars
        rep = JET.report_call(
            ReactionNetworkImporters.hasstripchars,
            (String,);
            target_modules = (ReactionNetworkImporters,)
        )
        @test length(JET.get_reports(rep)) == 0

        # Test make_shortsym
        rep = JET.report_call(
            ReactionNetworkImporters.make_shortsym,
            (String, Int);
            target_modules = (ReactionNetworkImporters,)
        )
        @test length(JET.get_reports(rep)) == 0

        # Test seek_to_block
        rep = JET.report_call(
            ReactionNetworkImporters.seek_to_block,
            (Vector{String}, Int, String);
            target_modules = (ReactionNetworkImporters,)
        )
        @test length(JET.get_reports(rep)) == 0
    end

    @testset "MatrixNetwork construction" begin
        # Test that MatrixNetwork can be constructed without JET errors
        t = Catalyst.default_t()
        @species S1(t) S2(t)
        @parameters k1 k2

        rateexprs = [k1, k2]
        substoich = [1 0; 0 1]
        prodstoich = [0 1; 1 0]

        mn = MatrixNetwork(rateexprs, substoich, prodstoich;
            species = [S1, S2], params = [k1, k2], t = t)

        # Verify the construction produces the expected type
        @test mn isa MatrixNetwork
        @test length(mn.rateexprs) == length(rateexprs)
        @test mn.substoich == substoich
        @test mn.prodstoich == prodstoich
    end

    @testset "ComplexMatrixNetwork construction" begin
        # Test that ComplexMatrixNetwork can be constructed without JET errors
        t = Catalyst.default_t()
        @species S1(t) S2(t)
        @parameters k1

        rateexprs = [k1]
        stoichmat = [1 0; 0 1]
        incidencemat = [-1; 1;;]

        cmn = ComplexMatrixNetwork(rateexprs, stoichmat, incidencemat;
            species = [S1, S2], params = [k1], t = t)

        # Verify the construction produces the expected type
        @test cmn isa ComplexMatrixNetwork
        @test length(cmn.rateexprs) == length(rateexprs)
        @test cmn.stoichmat == stoichmat
        @test cmn.incidencemat == incidencemat
    end
end
