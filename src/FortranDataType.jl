"""
# module FortranDataType



# Examples

```jldoctest
julia>
```
"""
module FortranDataType

using Rematch: MatchFailure

export guesstype,
    parseint,
    parsefloat,
    parsebool,
    parsestring,
    @f_str

const FORTRAN_INT = r"(?<=\s|^)([-+]?\d+)(?=\s|$)"
const FORTRAN_FLOAT = r"[-+]?\d*\.?\d+((:?[ed])[-+]?\d+)?"i
const FORTRAN_BOOL = r"(\.(true|false|t|f)\.)"i
const FORTRAN_STRING = r"[\'\"](.*)[\'\"]"
const FORTRAN_COMPLEX = r"\([-+]?\d*\.?\d+((:?[ed])[-+]?\d+)?,\s*[-+]?\d*\.?\d+((:?[ed])[-+]?\d+)?\)"i

macro f_str(s)
    return :($s)
end  # macro f_str

"""

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
function guesstype(str::AbstractString)
    # Must put `Complex` in front of `AbstractFloat`, or the complex number will be matched by 2
    # floats. Must put `String` in front of `Number`s, or strings containing numbers will be matched.
    for (regex, type) in zip((FORTRAN_STRING, FORTRAN_INT, FORTRAN_COMPLEX, FORTRAN_FLOAT, FORTRAN_BOOL), (String, Integer, Complex, AbstractFloat, Bool))
        if !isnothing(match(regex, str))
            if regex == FORTRAN_FLOAT
                # 'd' and no-'d' returns double precision, 'e' returns single precision
                occursin(r"e"i, str) ? (return Float32) : return Float64
            elseif regex == FORTRAN_COMPLEX
                r1, r2 = split(str, ",")
                T1, T2 = guesstype(r1[2:end]), guesstype(r2[1:end - 1])
                return Complex{Base.promote_type(T1, T2)}
            else
                return type
            end  # if-else
        end  # if
    end
    throw(MatchFailure("No type could be suggested for '$(str)'!"))
end  # function guesstype

function captured(regex, s)
    m = match(regex, s)
    isnothing(m) && throw(ParseError("Cannot parse Fortran data $(s)!"))
    return m.captures
end  # function captured

function parseint(::Type{T}, s) where {T <: Integer}
    captures = captured(FORTRAN_INT, s)
    return parse(T, captures[1])
end  # function parseint

function parsefloat(::Type{T}, s) where {T <: AbstractFloat}
    captures = captured(FORTRAN_FLOAT, s)
    length(captures) == 4 && return parse(T, string(captures[1], ".", captures[3], "e", captures[4]))
    return parse(T, string(captures[1], ".", captures[3]))
end  # function parsefloat

function parsecomplex(::Type{Tuple{T1, T2}}, s) where {T1 <: AbstractFloat, T2 <: AbstractFloat}
    r1, r2 = split(s, ",")
    a, b = parsefloat(T1, r1[2:end]), parsefloat(T2, r2[1:end - 1])
    return Complex(a, b)
end  # function parsecomplex

function parsebool(s)
    captures = captured(FORTRAN_BOOL, s)
    captures[1] in ("true", "t") && return true
    captures[1] in ("false", "f") && return false
end  # function parsebool

function parsestring(s)
    captures = captured(FORTRAN_STRING, s)
    return "$(captures[1])"
end  # function parsestring

end
