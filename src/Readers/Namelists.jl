"""
# module Namelists



# Examples

```jldoctest
julia>
```
"""
module Namelists

export read_namelist

function read_namelist(lines)
    result = Dict()
    for line in lines
        str = strip(line)
        # Use '=' as the delimiter, split the stripped line into a key and a value.
        # Skip this line if a line starts with '&' (namelist caption) or '!' (comment) or
        # this line is empty ('').
        isempty(str) || any(startswith(str, x) for x in ('!', '/', '&')) && continue
        k, v = split(str, '=', limit=2)
        k = strip(k)
        v = first(split(strip(rstrip(strip(v), ',')), '!'))  # Ignore trailing comma of the line
        result[k] = v
    end
    return result
end  # function read_namelist

end
