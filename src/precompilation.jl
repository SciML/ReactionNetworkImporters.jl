using SparseArrays: sparse

@setup_workload begin
    @compile_workload begin
        rates = [1.0]
        substoich = reshape([1], 1, 1)
        prodstoich = reshape([0], 1, 1)

        matrix_network = MatrixNetwork(rates, substoich, prodstoich)
        loadrxnetwork(matrix_network; name = :precompile_matrix)

        sparse_network = MatrixNetwork(
            rates, sparse(substoich), sparse(prodstoich)
        )
        loadrxnetwork(sparse_network; name = :precompile_sparse_matrix)

        complex_network = ComplexMatrixNetwork(
            rates, [1 0], reshape([-1, 1], 2, 1)
        )
        loadrxnetwork(complex_network; name = :precompile_complex_matrix)
    end
end
