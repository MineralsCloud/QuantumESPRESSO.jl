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
    atom::AbstractString
    mass::Float64
    pseudopotential::AbstractString
end  # struct AtomicSpecies

@with_kw struct AtomicSpeciesCard <: Card
    option = nothing
    data::AbstractVector{AtomicSpecies}
end  # struct AtomicSpeciesCard

struct AtomicPosition
    atom::AbstractString
    position::AbstractVector{Float64}
end  # struct AtomicPosition

@with_kw struct AtomicPositionCard <: Card
    option::AbstractString = "alat"; @assert option in allowed_options(AtomicPositionCard)
    data::AbstractVector{AtomicPosition}
end  # struct AtomicPositionCard

@with_kw struct CellParametersCard <: Card
    option::AbstractString = "alat"; @assert option in allowed_options(CellParametersCard)
    lattice::AbstractMatrix
end  # struct CellParametersCard

abstract type KPoint end

@with_kw struct MonkhorstPackGrid <: KPoint
    grid::AbstractVector{Int}; @assert length(grid) == 3
    offsets::AbstractVector{Int}; @assert length(offsets) == 3
end  # struct MonkhorstPackGrid

struct GammaPoint <: KPoint end

@with_kw struct SpecialKPoint <: KPoint
    coordinates::AbstractVector{Float64}; @assert length(coordinates) == 3
    weight::Real
end  # struct SpecialKPoint

@with_kw struct KPointsCard <: Card
    option::AbstractString = "tpiba"; @assert option in allowed_options(KPointsCard)
    points::AbstractVector{KPoint}
end  # struct KPointsCard

allowed_options(::Type{<: Card}) = nothing
allowed_options(::Type{AtomicPositionCard}) = ("alat", "bohr", "angstrom", "crystal", "crystal_sg")
allowed_options(::Type{CellParametersCard}) = ("alat", "bohr", "angstrom")
allowed_options(::Type{KPointsCard}) = ("tpiba", "automatic", "crystal", "gamma", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")

end
