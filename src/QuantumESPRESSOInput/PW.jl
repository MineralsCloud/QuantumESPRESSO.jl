"""
# module PW



# Examples

```jldoctest
julia>
```
"""
module PW

using Parameters: @with_kw
using Rematch: @match

using QuantumESPRESSO.Namelists
using QuantumESPRESSO.Namelists.PW
using QuantumESPRESSO.Cards
using QuantumESPRESSO.Cards.PW
using QuantumESPRESSO.QuantumESPRESSOInput

export PWInput,
    typefield,
    namelists,
    cards

@with_kw struct PWInput <: Input
    control::ControlNamelist = ControlNamelist()
    system::SystemNamelist
    electrons::ElectronsNamelist = ElectronsNamelist()
    ions::IonsNamelist = IonsNamelist()
    cell::CellNamelist = CellNamelist()
    atomicspecies::AtomicSpeciesCard
    atomicpositions::AtomicPositionCard
    kpoints::KPointsCard
    cellparameters::CellParametersCard
end  # struct PWInput

function typefield(type::T) where {T <: Union{Namelist, Card}}
    if T <: Namelist
        @match T begin
            ControlNamelist => :control
            SystemNamelist => :system
            ElectronsNamelist => :electrons
            IonsNamelist => :ions
            CellNamelist => :cell
        end
    else  # T <: Card
        @match T begin
            AtomicSpeciesCard => :atomicspecies
            AtomicPositionCard => :atomicpositions
            KPointsCard => :kpoints
            CellParametersCard => :cellparameters
        end
    end  # if-else
end  # function typefield

filter_field_by_supertype(obj, ::Type{T}) where {T} = filter(x->isa(x, T), map(x->getfild(ob, x), fieldnames(typeof(obj))) |> collect)

namelists(input::PWInput) = filter_field_by_supertype(input, Namelist)

cards(input::PWInput) = filter_field_by_supertype(input, Card)

end
