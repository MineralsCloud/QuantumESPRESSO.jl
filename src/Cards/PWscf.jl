"""
# module PWscf



# Examples

```jldoctest
julia>
```
"""
module PWscf

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
    allowed_options

struct AtomicSpecies
    atom::String
    mass::Float64
    pseudopotential::String
end  # struct AtomicSpecies

@with_kw struct AtomicSpeciesCard <: Card
    option = nothing
    data::Vector{AtomicSpecies}
end  # struct AtomicSpeciesCard

struct AtomicPosition
    atom::String
    position::Vector{Float64}
end  # struct AtomicPosition

@with_kw struct AtomicPositionCard <: Card
    option::String = "alat"; @assert option in allowed_options(AtomicPositionCard)
    data::Vector{AtomicPosition}
end  # struct AtomicPositionCard

@with_kw struct CellParametersCard <: Card
    option::String = "alat"
    lattice::Crystal = Crystal(eye(3)u"nm")
end  # struct CellParametersCard

abstract type KPoint end

@with_kw struct MonkhorstPackGrid <: KPoint
    grid::Vector{Int}; @assert length(grid) == 3
    offsets::Vector{Int}; @assert length(offsets) == 3
end  # struct MonkhorstPackGrid

struct GammaPoint <: KPoint end

@with_kw struct SpecialKPoint <: KPoint
    coordinates::Vector{Float64}; @assert length(coordinates) == 3
    weight::Real
end  # struct SpecialKPoint

@with_kw struct KPointsCard <: Card
    option::String = "tpiba"
    points::Vector{KPoint}
end  # struct KPointsCard

allowed_options(::Card) = nothing
allowed_options(::AtomicPositionCard) = ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
allowed_options(::CellParametersCard) = ("alat", "bohr", "angstrom")
allowed_options(::KPointsCard) = ("tpiba", "automatic", "crystal", "gamma", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")

end
