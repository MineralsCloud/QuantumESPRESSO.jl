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

export Namelist, to_dict

abstract type Namelist <: InputEntry end

QuantumESPRESSO.name(::Type{<:Namelist}) = error("Undefined name!")

function to_dict(nml::Namelist)::Dict{Symbol,Any}
    return type2dict(nml)
end # function to_dict

function QuantumESPRESSO.to_qe(nml::Namelist; indent::AbstractString = "    ")::String
    entries = to_dict(nml)
    content = "&$(name(typeof(nml)))\n"
    for (key, value) in entries
        if value isa Vector{<:Pair}
            for x in value
                content *= "$(indent)$(key)($(x.first)) = $(string(to_fortran(x.second)))\n"
            end
        else
            content *= "$(indent)$(key) = $(string(to_fortran(value)))\n"
        end
    end
    return content * "/\n"
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
