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
compstoichmat =[2  0  1  0  0  3;
                 0  1  1  0  0  0;
                 0  0  0  1  3  0]
incidencemat  = [-1   1   0   0   0;
                 1  -1   0   0   0;
                 0   0  -1   1   0;
                 0   0   1  -1   0;
                 0   0   0   0  -1;
                 0   0   0   0   1]
prn = loadrxnetwork(MatrixNetwork(), pars, substoich, prodstoich; params=pars) # dense version
@test rs == prn.rn
prn = loadrxnetwork(MatrixNetwork(), pars, sparse(substoich), sparse(prodstoich); params=pars) # sparse version
@test rs == prn.rn

prn = LoadReacCompNetwork( pars, compstoichmat, incidencemat; params=pars) # dense version
@test rs == prn
prn = LoadReacCompNetwork( pars, sparse(compstoichmat), sparse(incidencemat); params=pars) # sparse version
@test rs == prn


# version with hard coded rates (no parameter symbols)
rs = @reaction_network begin
    1., 2S1 --> S2
    2., S2 --> 2S1
    3., S1 + S2 --> S3
    4., S3 --> S1 + S2
    5., 3S3 --> 3S1
end 
prn = loadrxnetwork(MatrixNetwork(), convert.(Float64,1:5), substoich, prodstoich)  # dense version
@test rs == prn.rn
prn = loadrxnetwork(MatrixNetwork(), convert.(Float64,1:5), sparse(substoich), sparse(prodstoich))  # sparse version 
@test rs == prn.rn

prn = LoadReacCompNetwork(convert.(Float64,1:5), compstoichmat, incidencemat)  # dense version
@test rs == prn
prn = LoadReacCompNetwork(convert.(Float64,1:5), sparse(compstoichmat), sparse(incidencemat)) # sparse version
@test rs == prn



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
prn = loadrxnetwork(MatrixNetwork(), rates, substoich, prodstoich; species=species, params=pars)  # dense version
@test rs == prn.rn
prn = loadrxnetwork(MatrixNetwork(), rates, sparse(substoich), sparse(prodstoich); species=species, params=pars) # sparse version
@test rs == prn.rn

prn = LoadReacCompNetwork(rates, compstoichmat, incidencemat; species = species, params=pars)  # dense version
@test rs == prn
prn = LoadReacCompNetwork(rates, sparse(compstoichmat), sparse(incidencemat); species = species, params=pars)  # sparse version
@test rs == prn


