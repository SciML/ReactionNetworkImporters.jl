# History

## v1.0.0

Breaking release for compatibility with Catalyst.jl v16.

### Breaking changes

- `loadrxnetwork` now returns a `ReactionSystem` directly instead of a
  `ParsedReactionNetwork`. The `ParsedReactionNetwork` type has been removed.
- Initial conditions and parameter values from BNG files are stored as
  `initial_conditions` on the `ReactionSystem`, and are also available via the
  Catalyst metadata accessors `Catalyst.get_u0_map(rn)` and
  `Catalyst.get_parameter_map(rn)`.
- BNG species name mappings and group symbols are stored as system metadata,
  accessible via the new exported accessors `get_varstonames(rn)` and
  `get_groupstosyms(rn)`.
- `ModelingToolkit` is no longer a dependency; replaced by `ModelingToolkitBase`.
- Test dependency `CSVFiles` replaced by `CSV.jl`.
- Requires Catalyst.jl v16+.

### Migration guide

| Before (v0.x)                              | After (v1.0)                        |
| ------------------------------------------ | ----------------------------------- |
| `prn = loadrxnetwork(BNGNetwork(), file)`  | `rn = loadrxnetwork(BNGNetwork(), file)` |
| `prn.rn`                                   | `rn` (the return value itself)      |
| `prn.u0`                                   | `Catalyst.get_u0_map(rn)`          |
| `prn.p`                                    | `Catalyst.get_parameter_map(rn)`   |
| `prn.varstonames`                          | `get_varstonames(rn)`              |
| `prn.groupstosyms`                         | `get_groupstosyms(rn)`             |
| `rs == prn.rn`                             | `Catalyst.isequivalent(rs, rn)`    |
| `prn = loadrxnetwork(mn)`; `prn.rn`       | `rn = loadrxnetwork(mn)`           |
