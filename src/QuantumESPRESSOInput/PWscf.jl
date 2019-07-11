"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using Parameters: @with_kw

using QuantumESPRESSO.QuantumESPRESSOInput

export PWscfInput

@with_kw struct PWscfInput <: QuantumESPRESSOInput
    namelists::Dict
    cards::Dict
end  # struct PWscfInput

end
