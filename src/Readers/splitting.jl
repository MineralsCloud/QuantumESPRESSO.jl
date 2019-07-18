#=
splitting:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-17
=#
const NAMELIST_SEP = r"/\s*[\r\n]"
const NAMELIST_IDENTIFIERS = r"&CONTROL"i, r"&SYSTEM"i, r"&ELECTRONS"i, r"&IONS"i, r"&CELL"i
const CARD_IDENTIFIERS = r"ATOMIC_SPECIES"i, r"ATOMIC_POSITIONS"i, r"K_POINTS"i, r"CELL_PARAMETERS"i, r"OCCUPATIONS"i, r"CONSTRAINTS"i, r"ATOMIC_FORCES"i

struct RangeIndices
    startindex
    endindex
end  # struct RangeIndices

function get_namelist_identifier_positions(io::IOStream; include_heading::Bool=true, include_ending::Bool=false)
    match_records = Dict()
    for pattern in NAMELIST_IDENTIFIERS
        m0 = findfirst(pattern, io)
        m0 === nothing && continue
        m1 = findnext(NAMELIST_SEP, io)
        match_records[pattern] = RangeIndices(startindex=(include_heading ? m0.start() : m0.end()), endindex=(include_ending ? m1.end() : m1.start()))
    end
    m = sorted(match_records.items(), key=operator.itemgetter(1))
    return OrderedDict(m)
end  # function get_namelist_identifier_positions
