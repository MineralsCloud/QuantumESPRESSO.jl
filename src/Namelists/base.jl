#=
base:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-21
=#
using Parameters: type2dict, reconstruct

export Namelist

abstract type Namelist end

function to_dict(nml::Namelist)::Dict{Symbol,Any}
    return type2dict(nml)
end  # function to_dict

function evolve(nml::T, kwargs...)::T where {T <: Namelist}
    return reconstruct(nml, kwargs...)
end  # function evolve
function evolve(nml::T, dict::AbstractDict{Symbol,Any})::T where {T <: Namelist}
    return reconstruct(nml, dict)
end  # function evolve
