using Catalyst, ReactionNetworkImporters, SparseArrays

# version giving parameters and rates
rs = @reaction_network begin
    k1, 2S₁ --> S₂
    k2, S₂ --> 2S₁
    k3, S₁ + S₂ --> S₃
    k4, S₃ --> S₁ + S₂
    k5, 3S₃ --> 3S₁
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
prn = loadrxnetwork(MatrixNetwork(), "testnet", pars, substoich, prodstoich; params=pars)
@test rs == prn.rn

# sparse version
prn = loadrxnetwork(MatrixNetwork(), "testnet", pars, sparse(substoich), sparse(prodstoich); params=pars)
@test rs == prn.rn

# version with hard coded rates (no parameter symbols)
rs = @reaction_network begin
    1., 2S₁ --> S₂
    2., S₂ --> 2S₁
    3., S₁ + S₂ --> S₃
    4., S₃ --> S₁ + S₂
    5., 3S₃ --> 3S₁
end 
prn = loadrxnetwork(MatrixNetwork(), "testnet", convert.(Float64,1:5), substoich, prodstoich)
@test rs == prn.rn

# sparse version
prn = loadrxnetwork(MatrixNetwork(), "testnet", convert.(Float64,1:5), sparse(substoich), sparse(prodstoich))
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
prn = loadrxnetwork(MatrixNetwork(), "testnet2", rates, substoich, prodstoich; species=species, params=pars)
@test rs == prn.rn

# sparse version
prn = loadrxnetwork(MatrixNetwork(), "testnet2", rates, sparse(substoich), sparse(prodstoich); species=species, params=pars)
@test rs == prn.rn
