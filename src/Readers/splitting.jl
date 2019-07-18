#=
splitting:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-17
=#
using DataStructures
using FilePaths: AbstractPath

export get_namelist_identifier_indices,
    get_card_identifier_indices

const NAMELIST_END = r"/\s*[\r\n]"
const NAMELIST_STARTS = "&CONTROL", "&SYSTEM", "&ELECTRONS", "&IONS", "&CELL"  # regex: "&(.[^,]*)"
const CARD_STARTS = "ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS", "CELL_PARAMETERS", "OCCUPATIONS", "CONSTRAINTS", "ATOMIC_FORCES"

function get_card_identifier_indices(io::IOStream)
    records = OrderedDict()
    for (i, line) in enumerate(eachline(io))
        str = strip(line)
        isempty(str) || startswith(str, '!') || startswith(str, '#') && continue
        for cardname in CARD_STARTS
            if occursin(Regex("$cardname", "i"), str)
                records[cardname] = i
            else
                continue
            end  # if-else
        end  # for
    end  # for
    if haskey(records, "OCCUPATIONS")
        index = last(collect(values(get_namelist_identifier_indices(seek(io, 1)))))
        records["OCCUPATIONS"] < index && pop!(records, "OCCUPATIONS")
    end  # if
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
        for namelistname in NAMELIST_STARTS
            if occursin(Regex("$namelistname", "i"), str)
                records[namelistname] = i
            else
                continue
            end  # if-else
        end  # for
    end  # for
    return records
end  # function get_namelist_identifier_indices
function get_namelist_identifier_indices(path::AbstractPath)
    isfile(path) && isreadable(path) || error("File $(path) not readable!")
    open(path, "r") do io
        get_namelist_identifier_indices(io)
    end
end  # function get_card_identifier_indices
