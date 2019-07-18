#=
splitting:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-17
=#
using DataStructures: SortedDict

const NAMELIST_SEP = r"/\s*[\r\n]"
const NAMELIST_IDENTIFIERS = r"&CONTROL"i, r"&SYSTEM"i, r"&ELECTRONS"i, r"&IONS"i, r"&CELL"i
const CARD_IDENTIFIERS = r"ATOMIC_SPECIES"i, r"ATOMIC_POSITIONS"i, r"K_POINTS"i, r"CELL_PARAMETERS"i, r"OCCUPATIONS"i, r"CONSTRAINTS"i, r"ATOMIC_FORCES"i

function get_namelist_identifier_indices(io::IOStream)
    match_records = Dict()
    for pattern in NAMELIST_IDENTIFIERS
        m0 = findall(str -> occursin(pattern, str), io)
        m0 === nothing && continue
        m1 = findnext(str -> occursin(NAMELIST_SEP, str), io, m0)
        match_records[pattern] = range(m0, m1)
    end
    return SortedDict(sort(collect(match_records), by=x->x[2]))
end  # function get_namelist_identifier_indices
