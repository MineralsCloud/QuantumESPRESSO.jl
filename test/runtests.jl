using QuantumESPRESSO
using Test

@testset "QuantumESPRESSO.jl" begin
    include("NamelistsTests.jl")
    include("FortranDataTypeTests.jl")
end
