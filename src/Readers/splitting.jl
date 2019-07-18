#=
splitting:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-17
=#
using DataStructures: SortedDict, OrderedDict
using FilePaths: AbstractPath
using IterUtils: throw_which_occursin

export get_namelist_identifier_indices,
    get_card_identifier_indices

const NAMELIST_END = r"/\s*[\r\n]"
const NAMELIST_STARTS = r"&CONTROL"i, r"&SYSTEM"i, r"&ELECTRONS"i, r"&IONS"i, r"&CELL"i  # regex: "&(.[^,]*)"
const CARD_STARTS = r"ATOMIC_SPECIES"i, r"ATOMIC_POSITIONS"i, r"K_POINTS"i, r"CELL_PARAMETERS"i, r"OCCUPATIONS"i, r"CONSTRAINTS"i, r"ATOMIC_FORCES"i

function get_card_identifier_indices(io::IOStream)
    records = OrderedDict()
    for (i, line) in enumerate(eachline(io))
        str = strip(line)
        isempty(str) || startswith(str, '!') || startswith(str, '#') && continue
        cardname = throw_which_occursin(CARD_STARTS, str)
        cardname === nothing ? continue : records[cardname] = i
    end  # for
    return records
end  # function get_namelist_identifier_indices
function get_card_identifier_indices(path::AbstractPath)
    isfile(path) && isreadable(path) || error("File $(path) not readable!")
    open(path, "r") do io
        get_card_identifier_indices(io)
    end
end  # function get_card_identifier_indices

function get_namelist_identifier_indices(io::IOStream)
    records = OrderedDict()
    for (i, line) in enumerate(eachline(io))
        str = strip(line)
        isempty(str) || startswith(str, '!') || startswith(str, '#') && continue
        namelistname = throw_which_occursin(NAMELIST_STARTS, str)
        namelistname === nothing ? continue : records[namelistname] = i
    end  # for
    return records
end  # function get_namelist_identifier_indices
function get_namelist_identifier_indices(path::AbstractPath)
    isfile(path) && isreadable(path) || error("File $(path) not readable!")
    open(path, "r") do io
        get_namelist_identifier_indices(io)
    end
end  # function get_card_identifier_indices
