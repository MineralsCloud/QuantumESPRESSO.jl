#=
FortranDataTypeTests:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-20
=#
module FortranDataTypeTests

using Test

using QuantumESPRESSO.FortranDataType

@testset "Test `whattype`" begin
    @test whattype("1.0E-6") == Float32
    @test whattype("3.2767e+2") == Float32
    @test whattype("1.89e-14") == Float32
    @test whattype("-0.65e-2") == Float32
    @test whattype("+1e8") == Float32

    @test whattype("1.0D-6") == Float64
    @test whattype("1d8") == Float64
    @test whattype("0.") == Float64
    @test whattype("1.00") == Float64
    @test whattype("+3.141593") == Float64
    @test whattype("+3.1415926535d+0") == Float64
    @test whattype("-4.78d+6") == Float64
    @test whattype("1.0d+0") == Float64

    @test whattype("124") == Integer
    @test whattype("-448") == Integer
    @test whattype("0") == Integer
    @test whattype("32767") == Integer
    @test whattype("2147483647") == Integer
    @test whattype("-9874") == Integer

    @test whattype(".true.") == Bool
    @test whattype(".t.") == Bool
    @test whattype(".false.") == Bool
    @test whattype(".f.") == Bool

    @test whattype("''") == String
    @test whattype("\"./tmp234\"") == String
    @test whattype("'david'") == String
end  # testset

end
