#
# ComfyCommons.jl
# Copyright 2017, 2018 Mirko Bunse
#
#
# Comfortable utility functions.
#
#
# ComfyCommons.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ComfyCommons.jl.  If not, see <http://www.gnu.org/licenses/>.
#
module Yaml # sub-module of ComfyCommons


import YAML
export load_file, write_file, expand, interpolate, interpolate!


"""
    load_file(filename, more_constrictors=nothing; kwargs...)

See `YAML.load_file`, which is extended here by keyword arguments.

You can specify additional configurations with the keyword arguments. For example, if you
run `load_file("bla.yml", more = "stuff")`, you read the `bla.yml` file and set the `more`
property to `stuff`.
"""
function load_file(filename::AbstractString, more_constructors::YAML._constructor=nothing;
                   kwargs...)
    config = YAML.load_file(filename, more_constructors)
    for (k, v) in kwargs
        config[string(k)] = v
    end
    return config
end


#
# YAML file writing (may contribute to https://github.com/dcjones/YAML.jl/issues/29)
#
"""
    write_file(filename, config)

Write the given configuration to a YAML file.
"""
function write_file(filename::AbstractString, config::Dict{Any,Any}, prefix::AbstractString="")
    if (!endswith(filename, ".yml"))
        warn("The provided filename $filename does not end on '.yml'. Still writing...")
    end
    file = open(filename, "w")
    print(file, prefix)
    write(file, config)
    close(file)
end

function write(io::IO, config::Dict{Any,Any}, level::Int=0, ignorelevel::Bool=false)
    for (i, tup) in enumerate(config)
        write(io, tup, level, ignorelevel ? i == 1 : false)
    end
end

function write(io::IO, config::Pair{Any,Any}, level::Int=0, ignorelevel::Bool=false)
    print(io, indent(string(config[1]) * ":", level, ignorelevel)) # print key
    if (typeof(config[2]) <: Dict || typeof(config[2]) <: AbstractArray)
        print(io, "\n")
    else
        print(io, " ")
    end
    write(io, config[2], level + 1) # print value
end

function write(io::IO, config::AbstractArray, level::Int=0, ignorelevel::Bool=false)
    for (i, elem) in enumerate(config)
        print(io, indent("- ", level))   # print sequence element character '-'
        write(io, elem, level + 1, true) # print value, ignore first indent
    end
end

write(io::IO, config::Any, level::Int=0, ignorelevel::Bool=false) =
    print(io, string(config) * "\n") # no indent required

indent(str::AbstractString, level::Int, ignorelevel::Bool=false) =
    repeat("  ", ignorelevel ? 0 : level) * str
#
# end of file writing
#



_INTERPOLATE_DOC = """
    interpolate(configuration, property; kwargs...)
    interpolate!(configuration, property; kwargs...)

Interpolate the value of a `property` in some `configuration` with the values referenced therein.

For example, if the `configuration` specifies two properties `prop: "\$(foo)bar"` and
`foo: "bar"`, interpolating the value of `foo` into `prop` yields `prop: "barbar"`. Namely,
`\$(foo)` has been replaced by the value of `foo`. If only `prop: "\$(foo)bar"` is given but
no property `foo` is configured, you can specify it as a keyword argument, calling
`interpolate(conf, "prop", foo="bar")`, which yields the same result. In any case, `prop`
has to refer to `foo` with a dollar sign and enclosing brackets (`\$(foo)`).

The `property` may be a vector specifying a root-to-leaf path in the configuration tree.
"""
@doc _INTERPOLATE_DOC interpolate
@doc _INTERPOLATE_DOC interpolate!

function interpolate(conf::Dict{Any,Any}, property::AbstractArray; kwargs...)
    replace(_getindex(conf, property...), r"\$\([a-zA-Z_]+\)" => s -> begin
        s = s[3:end-1] # remove leading '$' and enclosing brackets
        if haskey(conf, s)
            conf[s]
        elseif haskey(kwargs, Symbol(s))
            kwargs[Symbol(s)]
        else
            error("Key $s not in config and not supplied as keyword argument.")
        end
    end)
end

interpolate(conf::Dict{Any,Any}, property::Any; kwargs...) =
    interpolate(conf, [property]; kwargs...)

interpolate!(conf::Dict{Any,Any}, property::AbstractArray; kwargs...) =
    _setindex!(conf, interpolate(conf, property; kwargs...), property...)

interpolate!(conf::Dict{Any,Any}, property::Any; kwargs...) =
    interpolate!(conf, [property]; kwargs...)


"""
    expand(configuration, property)

Expand the values of a `configuration` with the content of a `property` specified therein.

Expansion means that if the `property` is a vector with multiple elements, each element
defines a copy of the configuration where only this element is stored in the `property`.
Thus, a vector of configurations is returned, with all elements similar to the original
`configuration` but each containing another value of the `property`. The `property` may be
a vector specifying a root-to-leaf path in the configuration tree.
"""
expand(config::Dict{Any,Any}, property::Any) =
    [ Dict(config..., property => v) for v in _getindex(config, property) ]

# deeper in the configuration tree, it gets slightly more complex
expand(config::Dict{Any,Any}, property::AbstractArray) =
    [ begin
        c = deepcopy(config)
        _setindex!(c, v, property...)
        c
    end for v in vcat(_getindex(config, property...)...) ]

# multiple expansions
expand(config::Dict{Any,Any}, properties::Any...) = # Any also matches AbstractArrays
    vcat([ expand(expansion, properties[2:end]...) for expansion in expand(config, properties[1]) ]...)

# functions that complement the usual indexing with varargs
# (overriding Base.getindex and Base.setindex would screw up these methods)
@inbounds _getindex(val::Dict, keys::Any...) =
    if length(keys) > 1
        _getindex(val[keys[1]], keys[2:end]...) # descend into the val[keys[1]] sub-tree
    elseif haskey(val, keys[1])
        val[keys[1]]
    else
        Any[nothing] # return a dummy value - this will never be set (see _setindex! below)
    end

@inbounds _getindex(arr::AbstractArray, keys::Any...) =
    [ try _getindex(val, keys...) catch; end for val in arr ]

@inbounds _setindex!(val::Dict, value::Any, keys::Any...) =
    if length(keys) > 1
        _setindex!(val[keys[1]], value, keys[2:end]...)
    elseif haskey(val, keys[1]) # only update existing mappings because expansion never adds new mappings
        val[keys[1]] = value
    end

@inbounds _setindex!(arr::AbstractArray, value::Any, keys::Any...) =
    [ try _setindex!(val, value, keys...) catch; end for val in arr ]


end
