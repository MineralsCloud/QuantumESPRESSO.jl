#=
namelists:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-17
=#
export read_namelist

function read_namelist(lines)
    result = Dict()
    for line in lines
        s = strip(line)
        # Use '=' as the delimiter, split the stripped line into a key and a value.
        # Skip this line if a line starts with '&' (namelist caption) or '!' (comment) or
        # this line is empty ('').
        (startswith(s, '&') || startswith(s, '!') || startswith(s, '/') || isempty(s)) && continue
        k, v = split(s, '=', limit=2)
        k = strip(k)
        v = first(split(strip(rstrip(strip(v), ',')), '!'))  # Ignore trailing comma of the line
        result[k] = v
    end
    return result
end  # function read_namelist
