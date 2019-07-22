#=
base:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-22
=#
export Card,
    name,
    option,
    allowed_options,
    evolve

abstract type Card end

name(::Type{<: Card}) = error("Undefined name!")

option(card::Card) = getfield(card, :option)

allowed_options(::Type{<: Card}) = error("No allowed options defined!")

function evolve(card::Card, newdict::AbstractDict)
    :option in keys && error("If you want to change the option of a card, reconstruct a new one!")
    return reconstruct(card, newdict)
end  # function evolve
function evolve(card::Card; kwargs...)
    return evolve(card, kwargs)
end  # function evolve
