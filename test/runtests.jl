using SafeTestsets, Test

@time begin
    @time @safetestset "Quality Assurance" begin
        include("qa.jl")
    end
    @time @safetestset "BNG Birth-Death Test" begin
        include("test_nullrxs_odes.jl")
    end
    @time @safetestset "BNG Repressilator Test" begin
        include("test_repressilator_odes.jl")
    end
    @time @safetestset "BNG Higher Order Test" begin
        include("test_higherorder_odes.jl")
    end
    @time @safetestset "Matrix Input Test" begin
        include("test_mats.jl")
    end
end
