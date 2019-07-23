"""
# module PW



# Examples

```jldoctest
julia>
```
"""
module PW

using IterTools: fieldvalues
using Parameters: @with_kw

using QuantumESPRESSO
using QuantumESPRESSO.Namelists
using QuantumESPRESSO.Namelists.PW
using QuantumESPRESSO.Cards
using QuantumESPRESSO.Cards.PW
using QuantumESPRESSO.QuantumESPRESSOInput

export PWInput,
    namelists,
    cards

@with_kw struct PWInput <: Input
    control::ControlNamelist = ControlNamelist()
    system::SystemNamelist
    electrons::ElectronsNamelist = ElectronsNamelist()
    ions::IonsNamelist = IonsNamelist()
    cell::CellNamelist = CellNamelist()
    atomicspecies::AtomicSpeciesCard
    atomicpositions::AtomicPositionsCard
    kpoints::KPointsCard
    cellparameters::CellParametersCard
end  # struct PWInput

filter_field_by_supertype(obj, ::Type{T}) where {T} = filter(x->isa(x, T), map(x->getfild(ob, x), fieldnames(typeof(obj))) |> collect)

namelists(input::PWInput) = filter_field_by_supertype(input, Namelist)

cards(input::PWInput) = filter_field_by_supertype(input, Card)

function QuantumESPRESSO.to_qe(input::PWInput; indent::AbstractString = "    ", sep::AbstractString = " ")::String
    return join(map(to_qe, fieldvalues(input)), "\n")
end  # function to_qe

end
