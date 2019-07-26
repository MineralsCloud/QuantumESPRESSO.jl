#=
base:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-22
=#
export InputEntry, name, to_qe

abstract type InputEntry end

name(::Type{<:InputEntry}) = error("Undefined name!")

function to_qe(dict::AbstractDict; indent::AbstractString = "    ")::String
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
