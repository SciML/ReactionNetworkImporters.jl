using Catalyst, ReactionNetworkImporters, SparseArrays
using ModelingToolkit: nameof

# version giving parameters and rates
rs = @reaction_network rs1 begin
    k1, 2S1 --> S2
    k2, S2 --> 2S1
    k3, S1 + S2 --> S3
    k4, S3 --> S1 + S2
    k5, 3S3 --> 3S1
end

@parameters k1 k2 k3 k4 k5
pars = [k1, k2, k3, k4, k5]
substoich = [2 0 1 0 0;
             0 1 1 0 0;
             0 0 0 1 3]
prodstoich = [0 2 0 1 3;
              1 0 0 1 0;
              0 0 1 0 0]
compstoichmat = [2 0 1 0 0 3;
                 0 1 1 0 0 0;
                 0 0 0 1 3 0]
incidencemat = [-1 1 0 0 0;
                1 -1 0 0 0;
                0 0 -1 1 0;
                0 0 1 -1 0;
                0 0 0 0 -1;
                0 0 0 0 1]
mn1 = MatrixNetwork(pars, substoich, prodstoich; params = pars)
prn = loadrxnetwork(mn1; name = nameof(rs)) # dense version
@test rs == complete(prn.rn)
mn2 = MatrixNetwork(pars, sparse(substoich), sparse(prodstoich); params = pars)
prn = loadrxnetwork(mn2; name = nameof(rs)) # sparse version
@test rs == complete(prn.rn)

cmn1 = ComplexMatrixNetwork(pars, compstoichmat, incidencemat; params = pars)
prn = loadrxnetwork(cmn1; name = nameof(rs))
@test rs == complete(prn.rn)
cmn2 = ComplexMatrixNetwork(pars, sparse(compstoichmat), sparse(incidencemat);
    params = pars)
prn = loadrxnetwork(cmn2; name = nameof(rs))
@test rs == complete(prn.rn)

# version with hard coded rates (no parameter symbols)
rs = @reaction_network rs2 begin
    1.0, 2S1 --> S2
    2.0, S2 --> 2S1
    3.0, S1 + S2 --> S3
    4.0, S3 --> S1 + S2
    5.0, 3S3 --> 3S1
end
mn1 = MatrixNetwork(convert.(Float64, 1:5), substoich, prodstoich)
prn = loadrxnetwork(mn1; name = nameof(rs))   # dense version
@test rs == complete(prn.rn)
mn2 = MatrixNetwork(convert.(Float64, 1:5), sparse(substoich), sparse(prodstoich))
prn = loadrxnetwork(mn2; name = nameof(rs)) # sparse version
@test rs == complete(prn.rn)

cmn1 = ComplexMatrixNetwork(convert.(Float64, 1:5), compstoichmat, incidencemat)
prn = loadrxnetwork(cmn1; name = nameof(rs))
@test rs == complete(prn.rn)
cmn2 = ComplexMatrixNetwork(convert.(Float64, 1:5), sparse(compstoichmat),
    sparse(incidencemat))
prn = loadrxnetwork(cmn2; name = nameof(rs))
@test rs == complete(prn.rn)

# version with species names, rate functions and symbolic parameters
rs = @reaction_network rs3 begin
    k1 * A, 2A --> B
    k2, B --> 2A
    k3, A + B --> C
    k4, C --> A + B
    k5, 3C --> 3A
end
@variables t
@species A(t) B(t) C(t)
species = [A, B, C]
rates = [k1 * A, k2, k3, k4, k5]
mn1 = MatrixNetwork(rates, substoich, prodstoich; species = species, params = pars)
prn = loadrxnetwork(mn1; name = nameof(rs)) # dense version
@test rs == complete(prn.rn)
mn2 = MatrixNetwork(rates, sparse(substoich), sparse(prodstoich); species = species,
    params = pars)
prn = loadrxnetwork(mn2; name = nameof(rs)) # sparse version
@test rs == complete(prn.rn)

cmn1 = ComplexMatrixNetwork(rates, compstoichmat, incidencemat; species = species,
    params = pars)
prn = loadrxnetwork(cmn1; name = nameof(rs))
@test rs == complete(prn.rn)
cmn2 = ComplexMatrixNetwork(rates, sparse(compstoichmat), sparse(incidencemat);
    species = species, params = pars)
prn = loadrxnetwork(cmn2; name = nameof(rs))
@test rs == complete(prn.rn)
