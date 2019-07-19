"""
# module Namelists



# Examples

```jldoctest
julia>
```
"""
module Namelists

export Namelist

abstract type Namelist end

include("PW.jl")
include("PH.jl")

end
