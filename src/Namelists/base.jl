#=
base:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-21
=#
using FilePaths: AbstractPath, extension, exists
import JSON
using Parameters: type2dict

using QuantumESPRESSO.FortranDataType
using QuantumESPRESSO.Yaml

export Namelist,
    name,
    to_dict

abstract type Namelist end

name(::Type{<: Namelist}) = error("Undefined name!")

function to_dict(nml::Namelist)::Dict{Symbol,Any}
    return type2dict(nml)
end  # function to_dict

function to_qe(nml::Namelist, indent::AbstractString = "    ")::String
    entries = Dict(key => to_fortran(value) for (key, value) in to_dict(nml))
    """
    &$(name(nml))
    $(join(["$(indent)$(key) = $(value)" for (key, value) in entries], "\n"))
    /
    """
end  # function to_qe

function Base.dump(path::AbstractPath, nml::Namelist)
    exists(path) || touch(path)
    entries = Dict(key => to_fortran(value) for (key, value) in to_dict(nml))
    iswritable(path) || error("File $(path) not writable!")
    open(path, "r+") do io
        if extension(path) == "json"
            JSON.print(io, entries)
        elseif extension(path) == "yaml" || extension(path) == "yml"
            for p in pairs(entries)
                Yaml.write(io, p)
            end
        else
            error("Unknown extension type given!")
        end  # if-elseif-else
    end
end  # function Base.dump
