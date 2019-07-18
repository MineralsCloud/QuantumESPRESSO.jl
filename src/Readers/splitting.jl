#=
splitting:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-17
=#
using DataStructures: SortedDict

export get_namelist_identifier_indices

const NAMELIST_END = raw"/\s*[\r\n]"
const NAMELIST_STARTS = "&CONTROL", "&SYSTEM", "&ELECTRONS", "&IONS", "&CELL"
const CARD_STARTS = "ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS", "CELL_PARAMETERS", "OCCUPATIONS", "CONSTRAINTS", "ATOMIC_FORCES"

function get_namelist_identifier_indices(io::IOStream)
    match_records = Dict()
    str = read(io, String)
    for pattern in NAMELIST_STARTS
        m = match(Regex("($(pattern))(.*)($(NAMELIST_END))", "i"), str)
        m === nothing ? continue : match_records[pattern] = range(m.offsets[1]; stop=m.offsets[3])
    end
    return SortedDict(sort(collect(match_records), by=x->x[2]))
end  # function get_namelist_identifier_indices
