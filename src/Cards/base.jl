#=
base:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-22
=#
export Card,
    name,
    option,
    allowed_options

abstract type Card end

name(::Type{<: Card}) = error("Undefined name!")

option(card::Card) = getfield(card, :option)

allowed_options(::Type{<: Card}) = nothing
