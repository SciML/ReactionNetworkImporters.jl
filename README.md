# ReactionNetworkImporters.jl

This package provides importers to load reaction networks from several file formats. Currently it supports loading networks in the following formats:
1. The basic format used by the [RSSA](https://www.cosbi.eu/research/prototypes/rssa) group at COSBI in their [model collection](https://www.cosbi.eu/prototypes/jLiexDeBIgFV4zxwnKiW97oc4BjTtIoRGajqdUz4.zip).
2. A *subset* of the BioNetGen .net file format.

Imported networks can then be output to
1. [DiffEqBiological](https://github.com/JuliaDiffEq/DiffEqBiological.jl/), as a `min_reaction_network`.
