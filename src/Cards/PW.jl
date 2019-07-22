"""
# module PW



# Examples

```jldoctest
julia>
```
"""
module PW

using Parameters: @with_kw
using Setfield: @lens, set

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
    KPointsCard,
    allowed_options,
    name,
    evolve

# =============================== AtomicSpecies ============================== #
struct AtomicSpecies{A <: AbstractString, B <: Real, C <: AbstractString}
    atom::A
    mass::B
    pseudopotential::C
end  # struct AtomicSpecies

function evolve(data::AtomicSpecies, dict::Dict{Symbol, T}) where {T}
    for (k, v) in dict
        lens = @macroexpand(quote
            @lens _.$(esc(k))
        end)
        data = set(data, lens, v)
    end
    return data
end  # function evolve

@with_kw struct AtomicSpeciesCard{A <: AbstractVector{<: AtomicSpecies}} <: Card
    option::Nothing = nothing
    data::A
end  # struct AtomicSpeciesCard
# ============================================================================ #

# ============================== AtomicPosition ============================== #
@with_kw struct AtomicPosition{A <: AbstractString, B <: AbstractVector{<: Real}, C <: AbstractVector{Int}}
    atom::A
    pos::B; @assert length(pos) == 3
    if_pos::C = [1, 1, 1]; @assert length(if_pos) == 3
end  # struct AtomicPosition

@with_kw struct AtomicPositionsCard{A <: AbstractString, B <: AbstractVector{<: AtomicPosition}} <: Card
    option::A = "alat"; @assert option in allowed_options(AtomicPositionsCard)
    data::B
end  # struct AtomicPositionsCard
# ============================================================================ #

# ============================== CellParameters ============================== #
@with_kw struct CellParametersCard{A <: AbstractString, B <: AbstractMatrix} <: Card
    option::A = "alat"; @assert option in allowed_options(CellParametersCard)
    data::B; @assert size(data) == (3, 3)
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

@with_kw struct KPointsCard{A <: AbstractString, B <: AbstractVector{<: KPoint}} <: Card
    option::A = "tpiba"; @assert option in allowed_options(KPointsCard)
    data::B
end  # struct KPointsCard
# ============================================================================ #

# ================================== Methods ================================= #
allowed_options(::Type{<: AtomicPositionsCard}) = ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
allowed_options(::Type{<: CellParametersCard}) = ("alat", "bohr", "angstrom")
allowed_options(::Type{<: KPointsCard}) = ("tpiba", "automatic", "crystal", "gamma", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")

name(::Type{<: AtomicSpeciesCard}) = :atomicspecies
name(::Type{<: AtomicPositionsCard}) = :atomicpositions
name(::Type{<: KPointsCard}) = :kpoints
name(::Type{<: CellParametersCard}) = :cellparameters
# ============================================================================ #

end
