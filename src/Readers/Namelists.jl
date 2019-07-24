"""
# module Namelists



# Examples

```jldoctest
julia>
```
"""
module Namelists

using Compat: isnothing

using QuantumESPRESSO.FortranDataType
using QuantumESPRESSO.Namelists
using QuantumESPRESSO.Namelists.PW

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
        k, v = split(str, '=', limit = 2)
        k = strip(k)
        v = first(split(strip(rstrip(strip(v), ',')), '!'))  # Ignore trailing comma of the line
        # `result` is a `Dict{Symbol,Any}`, we need to parse `FortranCode` from QuantumESPRESSO's input
        # as type of the field of the namelist.
        result[Symbol(k)] = parse(fieldtype(T, Symbol(k)), FortranCode(v))
    end
    dict = merge(to_dict(T()), result)
    return T((dict[f] for f in fieldnames(T))...)
end  # function read_namelist

end
