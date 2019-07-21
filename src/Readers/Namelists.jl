"""
# module Namelists



# Examples

```jldoctest
julia>
```
"""
module Namelists

using Compat: isnothing

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
    return titlecase(namelist_name)
end  # function read_title_line

function read_namelist(lines)
    result = Dict()
    namelist_name = read_title_line(first(lines), r"&(\w+)"i)
    for line in Iterators.drop(lines, 1)  # Drop the title line
        str = strip(line)
        # Use '=' as the delimiter, split the stripped line into a key and a value.
        # Skip this line if a line starts with '!' (comment) or this line is empty ('').
        isempty(str) || any(startswith(str, x) for x in ('!', '/')) && continue
        k, v = split(str, '=', limit=2)
        k = strip(k)
        v = first(split(strip(rstrip(strip(v), ',')), '!'))  # Ignore trailing comma of the line
        result[k] = v
    end
    namelist = Symbol("$(namelist_name)Namelist")
    parameters = Symbol(join(["$k=$v" for (k, v) in result], ","))
    eval(quote
        $(namelist)($(parameters))
    end)
end  # function read_namelist

end
