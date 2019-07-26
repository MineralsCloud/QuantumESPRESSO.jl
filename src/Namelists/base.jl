#=
base:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-21
=#
using FilePaths: AbstractPath, extension, exists
import JSON
using Parameters: type2dict

using QuantumESPRESSO
using QuantumESPRESSO.FortranDataType

export Namelist, to_dict, dropdefault

abstract type Namelist <: InputEntry end

QuantumESPRESSO.name(::Type{<:Namelist}) = error("Undefined name!")

function to_dict(nml::Namelist)::Dict{Symbol,Any}
    return type2dict(nml)
end # function to_dict

function dropdefault(nml::Namelist)
    default = typeof(nml)()
    result = Dict{Symbol,Any}()
    for (k, v) in to_dict(nml)
        if v != getfield(default, :k)
            result[k] = v
        end
    end
    return result
end

function QuantumESPRESSO.to_qe(dict::AbstractDict; indent::AbstractString = "    ")::String
    for (key, value) in dict
        if value isa Vector{<:Pair}
            for x in value
                content *= "$(indent)$(key)($(x.first)) = $(string(to_fortran(x.second)))\n"
            end
        else
            content *= "$(indent)$(key) = $(string(to_fortran(value)))\n"
        end
    end
    return content * "/\n"
end

function QuantumESPRESSO.to_qe(nml::Namelist; indent::AbstractString = "    ")::String
    return "&$(name(typeof(nml)))\n" * to_qe(to_dict(nml); indent = indent)
end # function to_qe

function Base.dump(path::AbstractPath, nml::Namelist)
    exists(path) || touch(path)
    entries = Dict(key => to_fortran(value) for (key, value) in to_dict(nml))
    iswritable(path) || error("File $(path) not writable!")
    open(path, "r+") do io
        if extension(path) == "json"
            JSON.print(io, entries)
        elseif extension(path) == "yaml" || extension(path) == "yml"
            @warn "Currently not supported!"
        else
            error("Unknown extension type given!")
        end # if-elseif-else
    end
end # function Base.dump
