"""
# module NamelistReaders



# Examples

```jldoctest
julia>
```
"""
module NamelistReaders

function read_namelists(io::IOStream)
    result = Dict()
    for line in eachline(io)  # Read each line in the namelist until '/'
        s = strip(line)
        # Use '=' as the delimiter, split the stripped line into a key and a value.
        # Skip this line if a line starts with '&' (namelist caption) or '!' (comment) or
        # this line is empty ('').
        startswith(s, '&') || startswith(s, '!') or isempty(s) && continue
        k, v = split(s, '=')  # FIXME: maxsplit=1
        k = strip(k)
        v = split(strip(rstrip(',', strip(v))), '!') |> first  # Ignore trailing comma of the line
        push!(result[k],v)
    end
    return result
end  # function read_namelist

end
