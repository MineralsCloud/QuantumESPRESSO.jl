#=
base:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-21
=#
using Parameters: type2dict

export Namelist

abstract type Namelist end

function to_dict(nml::Namelist)::Dict{Symbol,Any}
    return type2dict(nml)
end  # function to_dict
