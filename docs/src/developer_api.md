# Developer API

`NetworkFileFormat` is a developer extension interface. It is stable for
packages that add an importer, but applications should normally use one of the
concrete formats documented on the [API](@ref) page.

## Extending an importer

To support another serialized reaction-network format:

1. Define a concrete subtype of `NetworkFileFormat`.
2. Implement `ReactionNetworkImporters.loadrxnetwork(::YourFormat, input; kwargs...)`.
3. Return a Catalyst `ReactionSystem`; forward `name` and other supported
   keywords to the `ReactionSystem` constructor where applicable.
4. Document the accepted `input`, the file-format rules, and any metadata the
   importer adds to the returned system.

The generic `loadrxnetwork` function is the extension point. Extensions must
not overload a concrete built-in importer type or rely on parser internals.

```@docs
ReactionNetworkImporters.NetworkFileFormat
```
