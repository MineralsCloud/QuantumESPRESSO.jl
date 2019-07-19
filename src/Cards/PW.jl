"""
# module PW



# Examples

```jldoctest
julia>
```
"""
module PW

using Crystals
using Parameters: @with_kw

using QuantumESPRESSO.Cards

export AtomicSpecies,
    AtomicSpeciesCard,
    AtomicPosition,
    AtomicPositionCard,
    CellParametersCard,
    KPoint,
    MonkhorstPackGrid,
    GammaPoint,
    SpecialKPoint,
    KPointsCard,
    option,
    allowed_options

# =============================== AtomicSpecies ============================== #
struct AtomicSpecies{A <: AbstractString, B <: Real, C <: AbstractString}
    atom::A
    mass::B
    pseudopotential::C
end  # struct AtomicSpecies

@with_kw struct AtomicSpeciesCard{A <: AbstractVector{AtomicSpecies}} <: Card
    option = nothing
    data::A
end  # struct AtomicSpeciesCard
# ============================================================================ #

# ============================== AtomicPosition ============================== #
@with_kw struct AtomicPosition{A <: AbstractString, B <: AbstractVector{<: Real}, C <: AbstractVector{Int}}
    atom::A
    pos::B; @assert length(pos) == 3
    if_pos::C = [1, 1, 1]; @assert length(if_pos) == 3
end  # struct AtomicPosition

@with_kw struct AtomicPositionCard{A <: AbstractString, B <: AbstractVector{AtomicPosition}} <: Card
    option::A = "alat"; @assert option in allowed_options(AtomicPositionCard)
    data::B
end  # struct AtomicPositionCard
# ============================================================================ #

# ============================== CellParameters ============================== #
@with_kw struct CellParametersCard{A <: AbstractString, B <: AbstractMatrix} <: Card
    option::A = "alat"; @assert option in allowed_options(CellParametersCard)
    data::B; @assert size(lattice) == (3, 3)
end  # struct CellParametersCard
# ============================================================================ #

# ================================== KPoint ================================== #
abstract type KPoint end

@with_kw struct MonkhorstPackGrid{A <: AbstractVector{Int}, B <: AbstractVector{Int}} <: KPoint
    grid::A; @assert length(grid) == 3
    offsets::B; @assert length(offsets) == 3
end  # struct MonkhorstPackGrid

struct GammaPoint <: KPoint end

@with_kw struct SpecialKPoint{A <: AbstractVector{Float64}, B <: Real} <: KPoint
    coordinates::A; @assert length(coordinates) == 3
    weight::B
end  # struct SpecialKPoint

@with_kw struct KPointsCard{A <: AbstractString, B <: AbstractVector{KPoint}} <: Card
    option::A = "tpiba"; @assert option in allowed_options(KPointsCard)
    points::B
end  # struct KPointsCard
# ============================================================================ #

# ================================== Methods ================================= #
option(card::Card) = getfield(card, :option)

allowed_options(::Type{<: Card}) = nothing
allowed_options(::Type{AtomicPositionCard}) = ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
allowed_options(::Type{CellParametersCard}) = ("alat", "bohr", "angstrom")
allowed_options(::Type{KPointsCard}) = ("tpiba", "automatic", "crystal", "gamma", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
# ============================================================================ #

end
