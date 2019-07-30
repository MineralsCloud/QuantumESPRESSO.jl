"""
# module BasicIO



# Examples

```jldoctest
julia>
```
"""
module BasicIO

using FilePaths: AbstractPath

macro iostream_to_lines(methodname)
    return quote
        function $(esc(methodname))(io::IOStream)
            $(esc(methodname))(readlines(io))
        end
    end
end  # macro iostream_to_lines

macro path_to_iostream(methodname, mode::AbstractString = "r")
    return quote
        function $(esc(methodname))(path::AbstractPath)
            isfile(path) && isreadable(path) || error("File $(path) not readable!")
            open(path, $mode) do io
                $(esc(methodname))(io)
            end
        end
    end
end  # macro path_to_iostream

include("write.jl")

end
