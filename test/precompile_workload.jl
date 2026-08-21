using ReactionNetworkImporters, Test

@testset "Representative public workflows" begin
    matrix_network = MatrixNetwork([1.0], reshape([1], 1, 1), reshape([0], 1, 1))
    matrix_system = loadrxnetwork(matrix_network; name = :precompile_matrix)
    @test nameof(matrix_system) == :precompile_matrix

    complex_network = ComplexMatrixNetwork(
        [1.0], [1 0], reshape([-1, 1], 2, 1)
    )
    complex_system = loadrxnetwork(complex_network; name = :precompile_complex_matrix)
    @test nameof(complex_system) == :precompile_complex_matrix
end
