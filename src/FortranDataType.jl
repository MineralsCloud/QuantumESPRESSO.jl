"""
# module FortranDataType



# Examples

```jldoctest
julia>
```
"""
module FortranDataType

using Compat: isnothing
using Rematch: MatchFailure

export FortranCode,
    @f_str,
    guesstype,
    to_fortran

const FORTRAN_INT = r"(?<=\s|^)([-+]?\d+)(?=\s|$)"
const FORTRAN_FLOAT = r"[-+]?\d*\.?\d+((:?[ed])[-+]?\d+)?"i
const FORTRAN_BOOL = r"\.(true|false|t|f)\."i
const FORTRAN_STRING = r"[\'\"](.*)[\'\"]"
const FORTRAN_COMPLEX = r"\([-+]?\d*\.?\d+((:?[ed])[-+]?\d+)?,\s*[-+]?\d*\.?\d+((:?[ed])[-+]?\d+)?\)"i

struct FortranCode{T <: AbstractString}
    data::T
end  # struct FortranCode

macro f_str(str)
    # Must escape the variable to return it to the calling context and executed there.
    return :(FortranCode($(esc(str))))
end  # macro f_str

"""
    guesstype(s::FortranCode)

# Examples
```jldoctest
julia> guesstype("1.0E-6")
Float32

julia> guesstype("1.0D-6")
Float64

julia> guesstype("1.039624557")
Float64
```
"""
function guesstype(s::FortranCode)
    str = s.data
    # Must put `Complex` in front of `AbstractFloat`, or the complex number will be matched by 2
    # floats. Must put `String` in front of `Number`s, or strings containing numbers will be matched.
    for (regex, type) in zip((FORTRAN_STRING, FORTRAN_INT, FORTRAN_COMPLEX, FORTRAN_FLOAT, FORTRAN_BOOL), (String, Integer, Complex, AbstractFloat, Bool))
        if !isnothing(match(regex, str))
            if regex == FORTRAN_FLOAT
                # 'd' and no-'d' returns double precision, 'e' returns single precision
                occursin(r"e"i, str) ? (return Float32) : return Float64
            elseif regex == FORTRAN_COMPLEX
                r1, r2 = split(str, ",")
                T1, T2 = guesstype(@f_str(r1[2:end])), guesstype(@f_str(r2[1:end - 1]))
                return Complex{Base.promote_type(T1, T2)}
            else
                return type
            end  # if-else
        end  # if
    end
    throw(MatchFailure("No type could be suggested for '$(str)'!"))
end  # function guesstype

function captured(regex, str)
    m = match(regex, str)
    isnothing(m) && throw(MatchFailure("Cannot match Fortran data $(str)!"))
    return m.captures
end  # function captured

function Base.parse(::Type{T}, s::FortranCode) where {T <: Integer}
    str = s.data
    captures = captured(FORTRAN_INT, str)
    isabstracttype(T) && return parse(Int, captures[1])
    return parse(T, captures[1])
end
function Base.parse(::Type{T}, s::FortranCode) where {T <: AbstractFloat}
    str = s.data
    parse(T, replace(str, r"d"i => "e"))
end
function Base.parse(::Type{Complex{T}}, s::FortranCode) where {T <: AbstractFloat}
    str = s.data
    r1, r2 = split(str, ",")
    a, b = parse(T, @f_str(r1[2:end])), parse(T, @f_str(r2[1:end - 1]))
    return Complex(a, b)
end
function Base.parse(::Type{Bool}, s::FortranCode) where {T <: AbstractFloat}
    str = s.data
    captures = captured(FORTRAN_BOOL, str)
    captures[1] in ("true", "t") && return true
    captures[1] in ("false", "f") && return false
end
function Base.parse(::Type{T}, s::FortranCode) where {T <: AbstractString}
    str = s.data
    captures = captured(FORTRAN_STRING, str)
    return "$(captures[1])"
end

function to_fortran(v::Int)
    FortranCode(string(v))
end  # function to_fortran
function to_fortran(v::Float32; scientific::Bool = false)
    str = string(v)
    scientific && return FortranCode(replace(str, r"f"i => "E"))
    return FortranCode(string(v))
end  # function to_fortran
function to_fortran(v::Float64; scientific::Bool = false)
    str = string(v)
    scientific && return FortranCode(replace(str, r"e"i => "D"))
    return FortranCode(string(v))
end  # function to_fortran
function to_fortran(v::Bool)
    v ? (return FortranCode(".true.")) : return FortranCode(".false.")
end  # function to_fortran
function to_fortran(v::AbstractString)
    return FortranCode("'$(v)'")
end  # function to_fortran

function Base.string(s::FortranCode)
    return string(s.data)
end  # function Base.string

end
