import ReactionNetworkImporters: NetworkFileFormat, loadrxnetwork

struct GenericOnlyFormat <: NetworkFileFormat end

function loadrxnetwork(::GenericOnlyFormat, payload; name = :generic_only)
    return (name = name, payload = payload)
end

@testset "NetworkFileFormat generic extension interface" begin
    format = GenericOnlyFormat()

    @test format isa NetworkFileFormat
    @test loadrxnetwork(format, :network) == (name = :generic_only, payload = :network)
    @test loadrxnetwork(format, :network; name = :custom) ==
          (name = :custom, payload = :network)
end
