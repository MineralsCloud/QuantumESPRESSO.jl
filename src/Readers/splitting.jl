#=
splitting:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-17
=#
using DataStructures: SortedDict, OrderedDict
using IterUtils: throw_which_occursin

export get_namelist_identifier_indices

const NAMELIST_END = raw"/\s*[\r\n]"
const NAMELIST_STARTS = "&CONTROL", "&SYSTEM", "&ELECTRONS", "&IONS", "&CELL"
const CARD_STARTS = "ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS", "CELL_PARAMETERS", "OCCUPATIONS", "CONSTRAINTS", "ATOMIC_FORCES"

function get_card_identifier_indices(io::IOStream)
    records = OrderedDict()
    for (i, line) in enumerate(eachline(io))
        str = strip(line)
        isempty(str) || startswith(str, r"[!#]") && continue
        cardname = throw_which_occursin(CARD_STARTS, str)
        cardname === nothing ? continue : records[cardname] = i
    end  # for
    return records
end  # function get_namelist_identifier_indices

function get_namelist_identifier_indices(io::IOStream)
    records = OrderedDict()
    for (i, line) in enumerate(eachline(io))
        str = strip(line)
        isempty(str) || startswith(str, r"[!#]") && continue
        namelistname = throw_which_occursin(NAMELIST_STARTS, str)
        namelistname === nothing ? continue : records[namelistname] = i
    end  # for
    return records
end  # function get_namelist_identifier_indices
