# ReactionNetworkImporters.jl

This package provides importers to load reaction networks into
[Catalyst.jl](https://docs.sciml.ai/Catalyst/stable/)
[`ReactionSystem`s](https://docs.sciml.ai/Catalyst/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem)
from several file formats. Currently, it supports loading networks in the
following formats:

 1. A *subset* of the BioNetGen .net file format.
 2. Networks represented by dense or sparse substrate and product stoichiometric
    matrices.
 3. Networks represented by dense or sparse complex stoichiometric and incidence matrices.

[SBMLToolkit.jl](https://docs.sciml.ai/SBMLToolkit/stable/) provides an
alternative for loading SBML files into Catalyst models, offering a much broader
set of supported features. It allows the import of models that include features
such as constant species, boundary condition species, events, constraint
equations, and more. SBML files can be generated from many standard modeling
tools, including BioNetGen, COPASI, and Virtual Cell.

## Installation

To install ReactionNetworkImporters.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("ReactionNetworkImporters")
```

## Contributing

  - Please refer to the
    [SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://github.com/SciML/ColPrac/blob/master/README.md)
    for guidance on PRs, issues, and other matters relating to contributing to SciML.

  - See the [SciML Style Guide](https://github.com/SciML/SciMLStyle) for common coding practices and other style decisions.
  - There are a few community forums:

      + The #diffeq-bridged and #sciml-bridged channels in the
        [Julia Slack](https://julialang.org/slack/)
      + The #diffeq-bridged and #sciml-bridged channels in the
        [Julia Zulip](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
      + On the [Julia Discourse forums](https://discourse.julialang.org)
      + See also [SciML Community page](https://sciml.ai/community/)

## Examples

### Loading a BioNetGen .net file

A simple network from the built-in BioNetGen bngl examples is the
[repressilator](https://github.com/SciML/ReactionNetworkImporters.jl/blob/65f5f23bc52a7eefe32ca2ccabef0298e8e5772d/data/repressilator/Repressilator.bngl).
The `generate_network` command in the bngl file outputs a reduced network
description, i.e. a
[.net](https://github.com/SciML/ReactionNetworkImporters.jl/blob/65f5f23bc52a7eefe32ca2ccabef0298e8e5772d/data/repressilator/Repressilator.net)
file, which can be loaded into a Catalyst `ReactionSystem` as:

```julia
using ReactionNetworkImporters
fname = "PATH/TO/Repressilator.net"
prnbng = loadrxnetwork(BNGNetwork(), fname)
```

Here `BNGNetwork` is a type specifying the file format that is being loaded.
`prnbng` is a `ParsedReactionNetwork` structure with the following fields:

  - `rn`, a Catalyst `ReactionSystem`
  - `u0`, a `Dict` mapping initial condition symbolic variables to numeric values
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
oprob = ODEProblem(rn, Float64[], (0.0, tf), Float64[])
sol = solve(oprob, Tsit5(), saveat = tf / 1000.0)
```

Note that we specify empty parameter and initial condition vectors, as these are
already stored in the generated `ReactionSystem`, `rn`. A `Dict` mapping each
symbolic species and parameter to its initial value or symbolic expression can
be obtained using `ModelingToolkit.defaults(rn)`.

See the [Catalyst documentation](https://docs.sciml.ai/Catalyst/stable/) for how to
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
rs = @reaction_network testnetwork begin
    k1, 2A --> B
    k2, B --> 2A
    k3, A + B --> C
    k4, C --> A + B
    k5, 3C --> 3A
end

# network from basic stoichiometry using ReactionNetworkImporters
@parameters k1 k2 k3 k4 k5
@variables t
@species A(t) B(t) C(t)
species = [A, B, C]
pars = [k1, k2, k3, k4, k5]
substoich = [2 0 1 0 0;
             0 1 1 0 0;
             0 0 0 1 3]
prodstoich = [0 2 0 1 3;
              1 0 0 1 0;
              0 0 1 0 0]
mn = MatrixNetwork(pars, substoich, prodstoich; species = species,
                   params = pars) # a matrix network
prn = loadrxnetwork(mn) # dense version

# test the two networks are the same
@assert rs == prn.rn

# network from reaction complex stoichiometry
stoichmat = [2 0 1 0 0 3;
             0 1 1 0 0 0;
             0 0 0 1 3 0]
incidencemat = [-1 1 0 0 0;
                1 -1 0 0 0;
                0 0 -1 1 0;
                0 0 1 -1 0;
                0 0 0 0 -1;
                0 0 0 0 1]
cmn = ComplexMatrixNetwork(pars, stoichmat, incidencemat; species = species,
                           params = pars)  # a complex matrix network
prn = loadrxnetwork(cmn; name = :testnetwork)

# test the two networks are the same
@assert rs == prn.rn
```

The basic usages are

```julia
mn = MatrixNetwork(rateexprs, substoich, prodstoich; species = Any[],
                   params = Any[], t = nothing)
prn = loadrxnetwork(mn::MatrixNetwork)

cmn = ComplexMatrixNetwork(rateexprs, stoichmat, incidencemat; species = Any[],
                           params = Any[], t = nothing)
prn = loadrxnetwork(cmn::ComplexMatrixNetwork)
```

Here `MatrixNetwork` and `ComplexMatrixNetwork` are the types, which select that
we are constructing a substrate/product stoichiometric matrix-based or a
reaction complex matrix-based stoichiometric representation as input. See the
[Catalyst.jl API](https://docs.sciml.ai/Catalyst/stable/api/catalyst_api/) for more
discussion on these matrix representations, and how Catalyst handles symbolic
reaction rate expressions. These two types have the following fields:

  - `rateexprs`, any valid
    [Symbolics.jl](https://docs.sciml.ai/Symbolics/stable/) expression for
    the rates, or any basic number type. This can be a hard-coded rate constant
    like `1.0`, a parameter like `k1` above, or a general Symbolics expression
    involving parameters and species like `k*A`.

  - matrix inputs

      + For `MatrixNetwork`

          * `substoich`, a number of species by number of reactions matrix, with entry
            `(i,j)` giving the stoichiometric coefficient of species `i` as a
            substrate in reaction `j`.
          * `prodstoich`, a number of species by number of reactions matrix, with entry
            `(i,j)` giving the stoichiometric coefficient of species `i` as a product
            in reaction `j`.

      + For `ComplexMatrixNetwork`

          * `stoichmat`, the complex stoichiometry matrix [defined
            here](https://docs.sciml.ai/Catalyst/stable/api/catalyst_api/#Catalyst.complexstoichmat).
          * `incidencemat`, the complex incidence matrix [defined
            here](https://docs.sciml.ai/Catalyst/stable/api/catalyst_api/#Catalyst.reactioncomplexes).
  - `species`, an optional vector of symbolic variables representing each species
    in the network. Can be constructed using the Catalyst.jl `@species` macro.
    Each species should be dependent on the same time variable (`t` in the example
    above).
  - `parameters`, a vector of symbolic variables representing each parameter in
    the network. Can be constructed with the
    [ModelingToolkit.jl](https://docs.sciml.ai/ModelingToolkit/stable/)
    `@parameters` macro. If no parameters are used, it is an optional keyword.
  - `t`, an optional Symbolics.jl variable representing time as the independent
    variable of the reaction network. If not provided, `Catalyst.default_t()` is
    used to determine the default time variable.

For both input types, `loadrxnetwork` returns a `ParsedReactionNetwork`, `prn`,
with only the field, `prn.rn`, filled in. `prn.rn` corresponds to the generated
[Catalyst.jl
`ReactionSystem`](https://docs.sciml.ai/Catalyst/stable/api/catalyst_api/#Catalyst.ReactionSystem)
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

## Reproducibility

```@raw html
<details><summary>The documentation of this SciML package was built using these direct dependencies,</summary>
```

```@example
using Pkg # hide
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>and using this machine and Julia version.</summary>
```

```@example
using InteractiveUtils # hide
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>
```

```@example
using Pkg # hide
Pkg.status(; mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```

```@eval
using TOML
using Markdown
version = TOML.parse(read("../../Project.toml", String))["version"]
name = TOML.parse(read("../../Project.toml", String))["name"]
link_manifest = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
                "/assets/Manifest.toml"
link_project = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
               "/assets/Project.toml"
Markdown.parse("""You can also download the
[manifest]($link_manifest)
file and the
[project]($link_project)
file.
""")
```
