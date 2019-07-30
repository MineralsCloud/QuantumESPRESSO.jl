"""
# module Namelists



# Examples

```jldoctest
julia>
```
"""
module Namelists

using Compat: isnothing
using Fortran90Namelists.FortranToJulia: FortranData

using QuantumESPRESSOBase.Namelists
using QuantumESPRESSOBase.Namelists.PW

export read_namelist

function read_title_line(title_line, regex)
    m = match(regex, title_line)
    if isnothing(m)
        # The first line should be '&<NAMELIST>', if it is not, either the regular expression
        # wrong or something worse happened.
        error("No match found in $(title_line)!")
    else
        namelist_name = m.captures[1]  # The first parenthesized subgroup will be `namelist_name`.
    end
    return lowercase(namelist_name)
end  # function read_title_line

function read_namelist(lines)
    result = Dict()
    namelist_name = read_title_line(first(lines), r"&(\w+)"i)
    regex = r"\s*([\w\d]+)(?:\((.*)\))?\s*=\s*([^,\n]*)"i
    T = Dict(
        "control" => ControlNamelist,
        "system" => SystemNamelist,
        "electrons" => ElectronsNamelist,
        "ions" => IonsNamelist,
        "cell" => CellNamelist
    )[namelist_name]
    for line in Iterators.drop(lines, 1)  # Drop the title line
        str = strip(line)
        # Use '=' as the delimiter, split the stripped line into a key and a value.
        # Skip this line if a line starts with '!' (comment) or this line is empty ('').
        isempty(str) || any(startswith(str, x) for x in ('!', '/')) && continue

        m = match(regex, str)
        isnothing(m) && error("Matching not found!")
        captures = m.captures
        k = Symbol(string(captures[1]))
        if !isnothing(captures[2])  # An entry with multiple values, e.g., `celldm[2] = 3.0`.
            val = parse(Float64, FortranData(string(captures[3])))
            # If `celldm` occurs before, push the new value, else create a vector of pairs.
            index = parse(Int, captures[2])
            haskey(result, k) ? result[k] = fillbyindex!(result[k], Pair(index, val)) : result[k] = fillbyindex!([], Pair(index, val))
        else
            v = FortranData(string(captures[3]))
            # `result` is a `Dict{Symbol,Any}`, we need to parse `FortranData` from QuantumESPRESSO's input
            # as type of the field of the namelist.
            result[k] = parse(fieldtype(T, k), v)
        end
    end
    dict = merge(to_dict(T()), result)
    return T((dict[f] for f in fieldnames(T))...)
end  # function read_namelist

function fillbyindex!(x::AbstractVector, p::Pair{Int, T}) where {T}
    index, value = p
    if isempty(x)
        x = Vector{Union{Missing, T}}(missing, index)
    else
        if index > length(x)
            append!(x, Vector{Union{Missing, T}}(missing, index - length(x)))
        else
        end
    end
    x[index] = value
    return x
end

end
