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
using QuantumESPRESSO.Cards

export AtomicSpecies,
    AtomicSpeciesCard,
    AtomicPosition,
    AtomicPositionsCard,
    CellParametersCard,
    KPoint,
    MonkhorstPackGrid,
    GammaPoint,
    SpecialKPoint,
    KPointsCard

# =============================== AtomicSpecies ============================== #
struct AtomicSpecies{A <: AbstractString,B <: Real,C <: AbstractString}
    atom::A
    mass::B
    pseudopotential::C
end

function QuantumESPRESSO.to_qe(data::AtomicSpecies; sep::AbstractString = " ")::String
    return join(map(string, fieldvalues(data)), sep)
end  # function to_qe

struct AtomicSpeciesCard{T <: AbstractVector{<: AtomicSpecies}} <: Card
    data::T
end

function QuantumESPRESSO.to_qe(card::AtomicSpeciesCard; indent::AbstractString = "    ", sep::AbstractString = " ")::String
    """
    ATOMIC_SPECIES
    $(join(["$(indent)$(to_qe(x; sep = sep))" for x in card.data], "\n"))
    """
end  # function to_qe
# ============================================================================ #

# ============================== AtomicPosition ============================== #
@with_kw struct AtomicPosition{A <: AbstractString,B <: AbstractVector{<: Real},C <: AbstractVector{Int}}
    atom::A
    pos::B; @assert length(pos) == 3
    if_pos::C = [1, 1, 1]; @assert length(if_pos) == 3
end

@with_kw struct AtomicPositionsCard{A <: AbstractString,B <: AbstractVector{<: AtomicPosition}} <: Card
    option::A = "alat"; @assert option in allowed_options(AtomicPositionsCard)
    data::B
end
# ============================================================================ #

# ============================== CellParameters ============================== #
@with_kw struct CellParametersCard{A <: AbstractString,B <: AbstractMatrix} <: Card
    option::A = "alat"; @assert option in allowed_options(CellParametersCard)
    data::B; @assert size(data) == (3, 3)
end
# ============================================================================ #

# ================================== KPoint ================================== #
abstract type KPoint end

@with_kw struct MonkhorstPackGrid{A <: AbstractVector{Int},B <: AbstractVector{Int}} <: KPoint
    grid::A; @assert length(grid) == 3
    offsets::B; @assert length(offsets) == 3
end

struct GammaPoint <: KPoint end

@with_kw struct SpecialKPoint{A <: AbstractVector{Float64},B <: Real} <: KPoint
    coordinates::A; @assert length(coordinates) == 3
    weight::B
end

@with_kw struct KPointsCard{A <: AbstractString,B <: AbstractVector{<: KPoint}} <: Card
    option::A = "tpiba"; @assert option in allowed_options(KPointsCard)
    data::B
    @assert begin
        if option == "automatic"
            eltype(data) <: MonkhorstPackGrid
        elseif option == "gamma"
            eltype(data) <: GammaPoint
        else  # option in ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
            eltype(data) <: SpecialKPoint
        end
    end
end
# ============================================================================ #

# ================================== Methods ================================= #
Cards.option(card::AtomicSpeciesCard) = nothing

Cards.allowed_options(::Type{<: AtomicPositionsCard}) = ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
Cards.allowed_options(::Type{<: CellParametersCard}) = ("alat", "bohr", "angstrom")
Cards.allowed_options(::Type{<: KPointsCard}) = ("tpiba", "automatic", "crystal", "gamma", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")

QuantumESPRESSO.name(::Type{<: AtomicSpeciesCard}) = :atomicspecies
QuantumESPRESSO.name(::Type{<: AtomicPositionsCard}) = :atomicpositions
QuantumESPRESSO.name(::Type{<: KPointsCard}) = :kpoints
QuantumESPRESSO.name(::Type{<: CellParametersCard}) = :cellparameters
# ============================================================================ #

end
