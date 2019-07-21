#=
FortranDataTypeTests:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-20
=#
module FortranDataTypeTests

using Test

using QuantumESPRESSO.FortranDataType

@testset "Test `guesstype`" begin
    @test guesstype("1.0E-6") == Float32
    @test guesstype("3.2767e+2") == Float32
    @test guesstype("1.89e-14") == Float32
    @test guesstype("-0.65e-2") == Float32
    @test guesstype("+1e8") == Float32

    @test guesstype("1.0D-6") == Float64
    @test guesstype("1d8") == Float64
    @test guesstype("0.") == Float64
    @test guesstype("1.00") == Float64
    @test guesstype("+3.141593") == Float64
    @test guesstype("+3.1415926535d+0") == Float64
    @test guesstype("-4.78d+6") == Float64
    @test guesstype("1.0d+0") == Float64

    @test guesstype("124") == Integer
    @test guesstype("-448") == Integer
    @test guesstype("0") == Integer
    @test guesstype("32767") == Integer
    @test guesstype("2147483647") == Integer
    @test guesstype("-9874") == Integer

    @test guesstype(".true.") == Bool
    @test guesstype(".t.") == Bool
    @test guesstype(".false.") == Bool
    @test guesstype(".f.") == Bool

    @test guesstype("''") == String
    @test guesstype("\"./tmp234\"") == String
    @test guesstype("'david'") == String
end  # testset

end
