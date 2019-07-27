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

filter_field_by_supertype(obj, ::Type{T}) where {T} = filter(x->isa(x, T), map(x->getfield(obj, x), fieldnames(typeof(obj))) |> collect)

namelists(input::PWInput) = filter_field_by_supertype(input, Namelist)

cards(input::PWInput) = filter_field_by_supertype(input, Card)

function QuantumESPRESSO.to_qe(input::PWInput; indent::AbstractString = "    ", sep::AbstractString = " ", debug::Bool = true)::String
    if debug
        return join(map(to_qe, fieldvalues(input)), "\n")
    else
        str = ""
        for namelist in namelists(input)
            str *= to_qe(to_dict(namelist))
        end
        for card in cards(input)
            str *= to_qe(card)
        end
        return str
    end # if
end  # function to_qe

end
