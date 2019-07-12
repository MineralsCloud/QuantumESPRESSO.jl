"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

using Parameters: @with_kw

using QuantumESPRESSO.Namelists
using QuantumESPRESSO.Namelists.PWscf
using QuantumESPRESSO.Cards
using QuantumESPRESSO.Cards.PWscf
using QuantumESPRESSO.QuantumESPRESSOInput

export PWscfInput,
    namelists,
    cards

@with_kw struct PWscfInput <: Input
    control::ControlNamelist = ControlNamelist()
    system::SystemNamelist
    electrons::ElectronsNamelist = ElectronsNamelist()
    ions::IonsNamelist = IonsNamelist()
    cell::CellNamelist = CellNamelist()
    atomicspecies::AtomicSpeciesCard
    atomicpositions::AtomicPositionCard
    kpoints::KPointsCard
    cellparameters::CellParametersCard
end  # struct PWscfInput

filter_field_by_supertype(obj, ::Type{T}) where {T} = filter(x->isa(x, T), map(x->getfild(ob, x), fieldnames(typeof(obj))) |> collect)

namelists(input::PWscfInput) = filter_field_by_supertype(input, Namelist)

cards(input::PWscfInput) = filter_field_by_supertype(input, Card)

end
