using Catalyst, ReactionNetworkImporters, SparseArrays

# version giving parameters and rates
rs = @reaction_network begin
    k1, 2S1 --> S2
    k2, S2 --> 2S1
    k3, S1 + S2 --> S3
    k4, S3 --> S1 + S2
    k5, 3S3 --> 3S1
end k1 k2 k3 k4 k5

@parameters k1 k2 k3 k4 k5
pars = [k1,k2,k3,k4,k5]
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
prn = loadrxnetwork(MatrixNetwork(), pars, substoich, prodstoich; params=pars)
@test rs == prn.rn

# sparse version
prn = loadrxnetwork(MatrixNetwork(), pars, sparse(substoich), sparse(prodstoich); params=pars)
@test rs == prn.rn

# version with hard coded rates (no parameter symbols)
rs = @reaction_network begin
    1., 2S1 --> S2
    2., S2 --> 2S1
    3., S1 + S2 --> S3
    4., S3 --> S1 + S2
    5., 3S3 --> 3S1
end 
prn = loadrxnetwork(MatrixNetwork(), convert.(Float64,1:5), substoich, prodstoich)
@test rs == prn.rn

# sparse version
prn = loadrxnetwork(MatrixNetwork(), convert.(Float64,1:5), sparse(substoich), sparse(prodstoich))
@test rs == prn.rn



# version with species names, rate functions and symbolic parameters
rs = @reaction_network begin
    k1*A, 2A --> B
    k2, B --> 2A
    k3, A + B --> C
    k4, C --> A + B
    k5, 3C --> 3A
end k1 k2 k3 k4 k5
@parameters t
@variables A(t) B(t) C(t)
species = [A,B,C]
rates = [k1*A,k2,k3,k4,k5]
prn = loadrxnetwork(MatrixNetwork(), rates, substoich, prodstoich; species=species, params=pars)
@test rs == prn.rn

# sparse version
prn = loadrxnetwork(MatrixNetwork(), rates, sparse(substoich), sparse(prodstoich); species=species, params=pars)
@test rs == prn.rn
