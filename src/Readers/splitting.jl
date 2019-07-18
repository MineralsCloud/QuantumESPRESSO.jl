#=
splitting:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-17
=#
using DataStructures
using FilePaths: AbstractPath

export namelist_identifier_linenumbers,
    card_identifier_linenumbers

const NAMELIST_END = r"/\s*[\r\n]"
const NAMELIST_STARTS = "&CONTROL", "&SYSTEM", "&ELECTRONS", "&IONS", "&CELL"  # regex: "&(.[^,]*)"
const CARD_STARTS = "ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS", "CELL_PARAMETERS", "OCCUPATIONS", "CONSTRAINTS", "ATOMIC_FORCES"

function card_identifier_linenumbers(io::IOStream)
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
        linenumber = last(collect(values(namelist_identifier_linenumbers(seek(io, 1)))))
        records["OCCUPATIONS"] < linenumber && pop!(records, "OCCUPATIONS")
    end  # if
    return records
end  # function card_identifier_linenumbers
function card_identifier_linenumbers(path::AbstractPath)
    isfile(path) && isreadable(path) || error("File $(path) not readable!")
    open(path, "r") do io
        get_card_identifier_indices(io)
    end
end  # function card_identifier_linenumbers

function namelist_identifier_linenumbers(io::IOStream)
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
end  # function namelist_identifier_linenumbers
function namelist_identifier_linenumbers(path::AbstractPath)
    isfile(path) && isreadable(path) || error("File $(path) not readable!")
    open(path, "r") do io
        namelist_identifier_linenumbers(io)
    end
end  # function namelist_identifier_linenumbers
