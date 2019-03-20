using DiffEqBiological, ReactionNetworkImporters

# using DiffEqBiological
rs = @reaction_network begin
    k1, 2S1 --> S2
    k2, S2 --> 2S1
    k3, S1 + S2 --> S3
    k4, S3 --> S1 + S2
    k5, 3S3 --> 3S1
end k1 k2 k3 k4 k5

# model using matrices
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

prn = loadrxnetwork(MatrixNetwork(), "testnet", pars, pars, substoich, prodstoich)
@test rs == prn.rn


# using DiffEqBiological
rs = @reaction_network begin
    k1, 2A --> B
    k2, B --> 2A
    k3, A + B --> C
    k4, C --> A + B
    k5, 3C --> 3A
end k1 k2 k3 k4 k5
prn = loadrxnetwork(MatrixNetwork(), "testnet2", [:A,:B,:C], pars, pars, substoich, prodstoich)
@test rs == prn.rn
