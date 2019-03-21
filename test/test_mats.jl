using DiffEqBiological, ReactionNetworkImporters, SparseArrays

# version giving parameters and rates
rs = @reaction_network begin
    k1, 2S1 --> S2
    k2, S2 --> 2S1
    k3, S1 + S2 --> S3
    k4, S3 --> S1 + S2
    k5, 3S3 --> 3S1
end k1 k2 k3 k4 k5

pars = [:k1, :k2, :k3, :k4, :k5]
substoich =[2 0 0;
            0 1 0;
            1 1 0;
            0 0 1;
            0 0 3]'
prodstoich = [0 1 0;
              2 0 0;
              0 0 1;
              1 1 0;
              3 0 0]'
prn = loadrxnetwork(MatrixNetwork(), "testnet", pars, substoich, prodstoich; params=pars)
@test rs == prn.rn

# sparse version
prn = loadrxnetwork(MatrixNetwork(), "testnet", pars, sparse(substoich), sparse(prodstoich); params=pars)
@test rs == prn.rn

# version with hard coded rates (no parameter symbols)
rs = @reaction_network begin
    1., 2S1 --> S2
    1., S2 --> 2S1
    1., S1 + S2 --> S3
    1., S3 --> S1 + S2
    1., 3S3 --> 3S1
end 
prn = loadrxnetwork(MatrixNetwork(), "testnet", ones(Float64,5), substoich, prodstoich)
@test rs == prn.rn

# sparse version
prn = loadrxnetwork(MatrixNetwork(), "testnet", ones(Float64,5), sparse(substoich), sparse(prodstoich))
@test rs == prn.rn



# version with species names and parameters
rs = @reaction_network begin
    k1, 2A --> B
    k2, B --> 2A
    k3, A + B --> C
    k4, C --> A + B
    k5, 3C --> 3A
end k1 k2 k3 k4 k5
species = [:A, :B, :C]
prn = loadrxnetwork(MatrixNetwork(), "testnet2", pars, substoich, prodstoich; species=species, params=pars)
@test rs == prn.rn

# sparse version
prn = loadrxnetwork(MatrixNetwork(), "testnet2", pars, sparse(substoich), sparse(prodstoich); species=species, params=pars)
@test rs == prn.rn
