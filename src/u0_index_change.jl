# Ensures that `prnbng.u₀` works.
function Base.getproperty(prnbng::ParsedReactionNetwork, name::Symbol)
    if name === :u₀
        return getfield(prnbng, :u0)
    else
        return getfield(prnbng, name)
    end
end

# Ensures that `prnbng.u₀ = ...` works.
function Base.setproperty!(prnbng::ParsedReactionNetwork, name::Symbol, x)
    if name === :u₀
        return setfield!(prnbng, :u0, x)
    else
        return setfield!(prnbng, name, x)
    end
end