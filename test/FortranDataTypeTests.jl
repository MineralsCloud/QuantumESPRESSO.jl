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

    @test guesstype("(5.229, -4.78)") == Complex{Float64}
    @test guesstype("(0.0,1.0)") == Complex{Float64}
    @test guesstype("(0.0,1)") == Complex{Real}
    @test guesstype("(3.2767e+2, -0.65e-2)") == Complex{Float32}

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

@testset "Test `parse*` functions" begin
    @test parsefloat(Float32, "1.0E-6") === 1.0f-6
    @test parsefloat(Float32, "3.2767e+2") === 3.2767f2
    @test parsefloat(Float32, "1.89e-14") === 1.89f-14
    @test parsefloat(Float32, "-0.65e-2") === -0.65f-2
    @test parsefloat(Float32, "+1e8") === 1f8

    @test parsefloat(Float64, "1.0D-6") === 1.0e-6
    @test parsefloat(Float64, "1d8") === 1e8
    @test parsefloat(Float64, "0.") === 0.
    @test parsefloat(Float64, "1.00") === 1.00
    @test parsefloat(Float64, "+3.141593") === 3.141593
    @test parsefloat(Float64, "+3.1415926535d+0") === +3.1415926535
    @test parsefloat(Float64, "-4.78d+6") === -4.78e6
    @test parsefloat(Float64, "1.0d+0") === 1.0

    @test parsecomplex(Complex{Float64}, "(5.229, -4.78)") === Complex(5.229, -4.78)
    @test parsecomplex(Complex{Float64}, "(0.0,1.0)") === Complex(0.0,1.0)
    @test parsecomplex(Complex{Float64}, "(0.0,1)") === Complex(0.0,1)
    @test parsecomplex(Complex{Float32}, "(3.2767e+2, -0.65e-2)") === Complex{Float32}(3.2767e2, -0.65e-2)

    @test parseint(Int, "124") === 124
    @test parseint(Int, "-448") === -448
    @test parseint(Int, "0") === 0
    @test parseint(Int, "32767") === 32767
    @test parseint(Int, "2147483647") === 2147483647
    @test parseint(Int, "-9874") === -9874

    @test parsebool(".true.") === true
    @test parsebool(".t.") === true
    @test parsebool(".false.") === false
    @test parsebool(".f.") === false

    @test parsestring("''") === ""
    @test parsestring("\"./tmp234\"") === "./tmp234"
    @test parsestring("'david'") === "david"
end  # testset

end
