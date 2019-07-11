"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using Parameters: @with_kw

using QuantumESPRESSO.Namelists.PWscf
using QuantumESPRESSO.Cards.PWscf
using QuantumESPRESSO.QuantumESPRESSOInput

export PWscfInput

@with_kw struct PWscfInput <: Input
    namelists::Dict; @assert all([
        isempty(setdiff([:CONTROL, :SYSTEM, :ELECTRONS], keys(namelists))),
        isa(namelists[:CONTROL], ControlNamelist),
        isa(namelists[:SYSTEM], SystemNamelist),
        isa(namelists[:ELECTRONS], ElectronsNamelist)
    ])
    cards::Dict; @assert all([
        isempty(setdiff([:ATOMIC_SPECIES, :ATOMIC_POSITIONS, :K_POINTS], keys(cards))),
        isa(cards[:ATOMIC_SPECIES], AtomicSpeciesCard),
        isa(cards[:ATOMIC_POSITIONS], AtomicPositionCard),
        isa(cards[:K_POINTS], KPointsCard)
    ])
end  # struct PWscfInput

end
