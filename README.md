# ReactionNetworkImporters.jl

[![Build Status](https://travis-ci.org/SciML/ReactionNetworkImporters.jl.svg?branch=master)](https://travis-ci.org/SciML/ReactionNetworkImporters.jl)
<!-- [![Build status](https://ci.appveyor.com/api/projects/status/wqq5flk2w8asad78/branch/master?svg=true)](https://ci.appveyor.com/project/SciML/reactionnetworkimporters-jl/branch/master) -->

This package provides importers to load reaction networks into
[Catalyst.jl](https://github.com/SciML/Catalyst.jl)
[`ReactionSystem`s](https://catalyst.sciml.ai/dev/api/catalyst_api/#ModelingToolkit.ReactionSystem)
from several file formats. Currently it supports loading networks in the
following formats:
1. A *subset* of the BioNetGen .net file format. Alternatively, add `writeSBML()` at the bottom of a BioNetGen model and import SBML with [SBMLToolkit.jl](https://github.com/SciML/SBMLToolkit.jl).
2. Networks represented by dense or sparse substrate and product stoichiometric
   matrices.
3. Networks represented by dense or sparse complex stoichiometric and incidence matrices.   
<!-- 3. The basic format used by the [RSSA](https://www.cosbi.eu/research/prototypes/rssa) group at COSBI in their [model collection](https://www.cosbi.eu/prototypes/jLiexDeBIgFV4zxwnKiW97oc4BjTtIoRGajqdUz4.zip). -->

----
## Examples

### Loading a BioNetGen .net file
A simple network from the builtin BioNetGen bngl examples is the
[repressilator](data/repressilator/Repressilator.bngl). The `generate_network`
command in the bngl file outputs a reduced network description, i.e. a
[.net](data/repressilator/Repressilator.net) file, which can be loaded into a
Catalyst `ReactionSystem` as:
```julia
using ReactionNetworkImporters
fname = "PATH/TO/Repressilator.net"
prnbng = loadrxnetwork(BNGNetwork(), fname)
```
Here `BNGNetwork` is a type specifying the file format that is being loaded.
`prnbng` is a `ParsedReactionNetwork` structure with the following fields:
- `rn`, a Catalyst `ReactionSystem`
- `u₀`, a `Dict` mapping initial condition symbolic variables to numeric values
  and/or symbolic expressions.
- `p`, a `Dict` mapping parameter symbolic variables to numeric values and/or
  symbolic expressions.
- `varstonames`, a `Dict` mapping the internal symbolic variable of a species
  used in the generated `ReactionSystem` to a `String` generated from the name
  in the .net file. This is necessary as BioNetGen can generate exceptionally
  long species names, involving characters that lead to malformed species names
  when used with `Catalyst`.
- `groupstosyms`, a `Dict` mapping the `String`s representing names for any
  groups defined in the BioNetGen file to the corresponding symbolic variable
  representing the `ModelingToolkit` symbolic observable associated with the
  group.

Given `prnbng`, we can construct and solve the corresponding ODE model for the
reaction system by
```julia
using OrdinaryDiffEq, Catalyst
rn = prnbng.rn
tf = 100000.0
oprob = ODEProblem(rn, Float64[], (0.,tf), Float64[])
sol = solve(oprob, Tsit5(), saveat=tf/1000.)
```
Note that we specify empty parameter and initial condition vectors as these are
already stored in the generated `ReactionSystem`, `rn`. A `Dict` mapping each
symbolic species and parameter to its initial value or symbolic expression can
be obtained using `ModelingToolkit.defaults(rn)`. 

See the [Catalyst documentation](https://catalyst.sciml.ai/dev/) for how to
generate ODE, SDE, jump and other types of models.

### Loading a matrix representation
Catalyst `ReactionSystem`s can also be constructed from 
- substrate and product stoichiometric matrices. 
- complex stoichiometric and incidence matrices.

For example, here we both directly build a Catalyst
network using the `@reaction_network` macro, and then show how to build the same
network from these matrices using `ReactionNetworkImporters`:
```julia
# Catalyst network from the macro:
rs = @reaction_network begin
    k1, 2A --> B
    k2, B --> 2A
    k3, A + B --> C
    k4, C --> A + B
    k5, 3C --> 3A
end k1 k2 k3 k4 k5

# network from basic stoichiometry using ReactionNetworkImporters
@parameters k1 k2 k3 k4 k5
@variables t A(t) B(t) C(t)
species = [A,B,C]
pars = [k1,k2,k3,k4,k5]
substoich =[ 2  0  1  0  0;
            0  1  1  0  0;
            0  0  0  1  3]
prodstoich =  [0  2  0  1  3;
                1  0  0  1  0;
                0  0  1  0  0]
mn= MatrixNetwork(pars, substoich, prodstoich; species=species, 
                  params=pars) # a matrix network
prn = loadrxnetwork(mn) # dense version

# test the two networks are the same
@assert rs == prn.rn

# network from reaction complex stoichiometry 
stoichmat =[2  0  1  0  0  3;
                 0  1  1  0  0  0;
                 0  0  0  1  3  0]
incidencemat  = [-1   1   0   0   0;
                 1  -1   0   0   0;
                 0   0  -1   1   0;
                 0   0   1  -1   0;
                 0   0   0   0  -1;
                 0   0   0   0   1]
cmn= ComplexMatrixNetwork(pars, stoichmat, incidencemat; species=species, 
                          params=pars)  # a complex matrix network
prn = loadrxnetwork(cmn)

# test the two networks are the same
@assert rs == prn.rn
```

The basic usages are
```julia
mn = MatrixNetwork(rateexprs, substoich, prodstoich; species=Any[], 
                   params=Any[], t=nothing)
prn = loadrxnetwork(mn::MatrixNetwork) 

cmn = ComplexMatrixNetwork(rateexprs, stoichmat, incidencemat; species=Any[], 
                           params=Any[], t=nothing)
prn = loadrxnetwork(cmn::ComplexMatrixNetwork)
```
Here `MatrixNetwork` and `ComplexMatrixNetwork` are the types, which select that
we are constructing a substrate/product stoichiometric matrix-based or a
reaction complex matrix-based stoichiometric representation as input. See the
[Catalyst.jl API](https://catalyst.sciml.ai/dev/api/catalyst_api/) for more
discussion on these matrix representations, and how Catalyst handles symbolic
reaction rate expressions. These two types have the following fields:
- `rateexprs`, any valid
  [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) expression for
  the rates, or any basic number type. This can be a hardcoded rate constant
  like `1.0`, a parameter like `k1` above, or an general Symbolics expression
  involving parameters and species like `k*A`. 
- matrix inputs 
  - For `MatrixNetwork`
    - `substoich`, a number of species by number of reactions matrix with entry
      `(i,j)` giving the stoichiometric coefficient of species `i` as a
      substrate in reaction `j`. 
    - `prodstoich`, a number of species by number of reactions matrix with entry
      `(i,j)` giving the stoichiometric coefficient of species `i` as a product
      in reaction `j`.
  - For `ComplexMatrixNetwork`  
    - `stoichmat`, the complex stoichiometry matrix [defined
      here](https://catalyst.sciml.ai/dev/api/catalyst_api/#Catalyst.complexstoichmat).
    - `incidencemat`, the complex incidence matrix [defined
      here](https://catalyst.sciml.ai/dev/api/catalyst_api/#Catalyst.reactioncomplexes).      
- `species`, an optional vector of symbolic variables representing each species
  in the network. Can be constructed using the Symbolics.jl `@variables` macro.
  Each species should be dependent on the same time variable (`t` in the example
  above). 
- `parameters`, a vector of symbolic variables representing each parameter in
  the network. Can be constructed with the
  [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl)
  `@parameters` macro. If no parameters are used it is an optional keyword.
- `t`, an optional Symbolics.jl variable representing time as the independent
  variable of the reaction network. If not provided `Catalyst.DEFAULT_IV` is
  used to determine the default time variable.

For both input types, `loadrxnetwork` returns a `ParsedReactionNetwork`, `prn`,
with only the field, `prn.rn`, filled in. `prn.rn` corresponds to the generated
[Catalyst.jl
`ReactionSystem`](https://catalyst.sciml.ai/dev/api/catalyst_api/#Catalyst.ReactionSystem)
that represents the network.

Dispatches are added if `substoich` and `prodstoich` both have the type
`SparseMatrixCSC`in case of `MatrixNetwork` (or `stoichmat` and `incidencemat`
both have the type `SparseMatrixCSC` in case of `ComplexMatrixNetwork`), in
which case they are efficiently iterated through using the `SparseArrays`
interface.

If the keyword argument `species` is not set, the resulting reaction network
will simply name the species `S1`, `S2`,..., `SN` for a system with `N` total
species. `params` defaults to an empty vector, so that it does not need to be
set for systems with no parameters.

<!-- ### Loading a RSSA format network file
As the licensing is unclear we can not redistribute any example RSSA formatted networks. They can be downloaded from the model collection link listed above. Assuming you've saved both a reaction network file and corresponding initial condition file, they can be loaded as
```julia
initialconditionf = "PATH/TO/FILE"
networkf = "PATH/TO/FILE"
rssarn = loadrxnetwork(RSSANetwork(), "RSSARxSys", initialconditionf, networkf)
```
Here `RSSANetwork` specifies the type of the file to parse, and `RSSARxSys` gives the type of the generated `reaction_network`. `rssarn` is again a `ParsedReactionNetwork`, but only the `rn` and `u₀` fields will now be relevant (the remaining fields will be set to `nothing`). -->
