#=
base:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-22
=#
export InputEntry, name, to_qe

abstract type InputEntry end

name(::Type{<:InputEntry}) = error("Undefined name!")

function to_qe end
