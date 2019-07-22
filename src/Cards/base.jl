#=
base:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-22
=#
using Parameters
using SimpleTraits

export Card,
    HasOption,
    name,
    option,
    allowed_options

abstract type Card end

@traitdef HasOption{T}

name(::Type{<: Card}) = error("Undefined name!")

@traitfn option(card::T) where {T; HasOption{T}} = getfield(card, :option)
@traitfn option(card::T) where {T; !HasOption{T}} = nothing

allowed_options(::Type{<: Card}) = error("No allowed options defined!")

function Parameters.reconstruct(card::Card, newdict::AbstractDict)
    :option in keys && error("If you want to change the option of a card, reconstruct a new one!")
    return reconstruct(card, newdict)
end  # function reconstruct
