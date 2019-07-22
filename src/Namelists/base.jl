#=
base:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-21
=#
using FilePaths: AbstractPath, extension
import JSON
using Parameters: type2dict, reconstruct
import YAML

export Namelist,
    to_dict,
    evolve

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

function to_qe(nml::Namelist, indent::AbstractString = "    ")::String
    entries = Dict(key => to_fortran(value) for (key, value) in to_dict(nml))
    """
    &$(name(nml))
    $(join(["$(indent)$(key) = $(value)" for (key, value) in entries], "\n"))
    /
    """
end  # function to_qe

function dump(path::AbstractPath, nml::Namelist)
    entries = Dict(key => to_fortran(value) for (key, value) in to_dict(nml))
    iswritable(path) || error("File $(path) not writable!")
    open(path, "r+") do io
        if extension(path) == "json"
            JSON.print(io, entries)
        elseif extension(path) == "yaml" || extension(path) == "yml"
            YAML.dump(io, entries)
        else
            error("Unknown extension type given!")
        end  # if-elseif-else
    end
end  # function dump
