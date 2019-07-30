#=
write:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-29
=#
using QuantumESPRESSOBase
using QuantumESPRESSOBase.QuantumESPRESSOInput.PW

function Base.write(io::IO, pw::PWInput, debug::Bool = true)
    write(io, to_qe(pw, debug = debug))
end
function Base.write(path::AbstractPath, pw::PWInput, debug::Bool = true)
    open(path, "r+") do io
        write(io, pw, debug)
    end
end
