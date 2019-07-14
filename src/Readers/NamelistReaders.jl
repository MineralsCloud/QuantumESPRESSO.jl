"""
# module NamelistReaders



# Examples

```jldoctest
julia>
```
"""
module NamelistReaders

ENV["PYTHON"] = joinpath(@__FILE__, "../../deps/f90nml/bin/python")
ENV["PYCALL_JL_RUNTIME_PYTHON"] = joinpath(@__FILE__, "../../deps/f90nml/bin/python")

println(Sys.which("python"))

using PyCall

pyimport("f90nml")

end
