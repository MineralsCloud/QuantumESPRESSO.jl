"""
# module FortranDataType



# Examples

```jldoctest
julia>
```
"""
module FortranDataType

using Rematch: MatchFailure

export whattype,
    parseint,
    parsefloat,
    parsebool,
    parsestring

const FORTRAN_INT = r"(?<=\s|^)([-+]?\d+)(?=\s|$)"
const FORTRAN_FLOAT = r"[-+]?\d*\.?\d+((:?[ed])[-+]?\d+)?"i
const FORTRAN_BOOL = r"(\.(true|false|t|f)\.)"i
const FORTRAN_STRING = r"[\'\"](.*)[\'\"]"
const FORTRAN_COMPLEX = r"\([-+]?\d*\.?\d+((:?[ed])[-+]?\d+)?,\s*[-+]?\d*\.?\d+((:?[ed])[-+]?\d+)?\)"i

"""

# Examples
```jldoctest
julia> whattype("1.0E-6")
Float32

julia> whattype("1.0D-6")
Float64

julia> whattype("1.039624557")
Float64
```
"""
function whattype(str::AbstractString)
    for (regex, type) in zip((FORTRAN_INT, FORTRAN_FLOAT, FORTRAN_COMPLEX, FORTRAN_BOOL, FORTRAN_STRING), (Integer, AbstractFloat, Complex, Bool, String))
        if !isnothing(match(regex, str))
            if regex == FORTRAN_FLOAT
                # 'd' and no-'d' returns double precision, 'e' returns single precision
                occursin(r"e"i, str) ? (return Float32) : return Float64
            else
                return type
            end  # if-else
        end  # if
    end
    throw(MatchFailure("No type could be suggested for '$(str)'!"))
end

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
    captures = captured(FORTRAN_INT, s)
    length(captures) == 4 && return parse(T, string(captures[1], ".", captures[3], "e", captures[4]))
    return parse(T, string(captures[1], ".", captures[3]))
end  # function parsefloat

function parsebool(s)
    captures = captured(FORTRAN_INT, s)
    captures[1] in ("true", "t") && return true
    captures[1] in ("false", "f") && return false
end  # function parsebool

function parsestring(s)
    captures = captured(FORTRAN_STRING, s)
    return "$(captures[1])"
end  # function parsestring

end
