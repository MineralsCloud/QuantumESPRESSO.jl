#=
splitting:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-17
=#
using DataStructures
using FilePaths: AbstractPath
using ResumableFunctions

export namelist_identifier_linenumbers,
    namelist_lineranges,
    card_identifier_linenumbers,
    input_identifier_linenumbers,
    dispatch_readers

const NAMELIST_END = '/'  # Not a regex anymore, since I strip everyline
const NAMELIST_STARTS = "&CONTROL", "&SYSTEM", "&ELECTRONS", "&IONS", "&CELL"  # regex: "&(.[^,]*)"
const CARD_STARTS = "ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS", "CELL_PARAMETERS", "OCCUPATIONS", "CONSTRAINTS", "ATOMIC_FORCES"

function namelist_identifier_linenumbers(io::IOStream)
    records = OrderedDict()
    for (i, line) in enumerate(eachline(io))
        str = strip(line)
        isempty(str) || startswith(str, '!') || startswith(str, '#') && continue
        for namelistname in NAMELIST_STARTS
            occursin(Regex("$namelistname", "i"), str) ? records[namelistname] = i : continue
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

function namelist_lineranges(io::IOStream)
    records = OrderedDict()
    for (i, line) in enumerate(eachline(io))
        str = strip(line)
        isempty(str) || startswith(str, '!') || startswith(str, '#') && continue
        for namelistname in NAMELIST_STARTS
            occursin(Regex("$namelistname", "i"), str) ? records[namelistname] = i : continue
        end  # for
        if startswith(str, NAMELIST_END)
            lastkey = last(collect(keys(records)))
            records[lastkey] = range(records[lastkey]; stop = i)
        end  # if
    end  # for
    return records
end  # function namelist_lineranges
function namelist_lineranges(path::AbstractPath)
    isfile(path) && isreadable(path) || error("File $(path) not readable!")
    open(path, "r") do io
        namelist_lineranges(io)
    end
end  # function namelist_lineranges

function card_identifier_linenumbers(io::IOStream)
    records = OrderedDict()
    for (i, line) in enumerate(eachline(io))
        str = strip(line)
        isempty(str) || startswith(str, '!') || startswith(str, '#') && continue
        for cardname in CARD_STARTS
            occursin(Regex("$cardname", "i"), str) ? records[cardname] = i : continue
        end  # for
    end  # for
    if haskey(records, "OCCUPATIONS")
        # Remember to rewind the `io`
        linenumber = last(collect(values(namelist_identifier_linenumbers(seekstart(io)))))
        records["OCCUPATIONS"] < linenumber && pop!(records, "OCCUPATIONS")
    end  # if
    return records
end  # function card_identifier_linenumbers
function card_identifier_linenumbers(path::AbstractPath)
    isfile(path) && isreadable(path) || error("File $(path) not readable!")
    open(path, "r") do io
        card_identifier_linenumbers(io)
    end
end  # function card_identifier_linenumbers

function input_identifier_linenumbers(io::IOStream)
    # Remember to rewind the `io`
    Dict("namelists" => namelist_identifier_linenumbers(io), "cards" => card_identifier_linenumbers(seekstart(io)))
end  # function input_identifier_linenumbers
function input_identifier_linenumbers(path::AbstractPath)
    isfile(path) && isreadable(path) || error("File $(path) not readable!")
    open(path, "r") do io
        input_identifier_linenumbers(io)
    end
end  # function input_identifier_linenumbers

function dispatch_readers(io::IOStream)
    linenumbers = input_identifier_linenumbers(io)
    namelist_linenumbers = linenumbers["namelists"]
    card_linenumbers = linenumbers["cards"]
    namelists = Dict()
    cards = Dict()
    for (k, v) in namelist_linenumbers
        namelists[k] = read_namelist(io[v])
    end  # for
    for (k, v) in card_linenumbers
        card[k] = begin
            k == "ATOMIC_SPECIES" && read_atomicspecies(io[v])
            k == "ATOMIC_POSITIONS" && read_atomicpositions(io[v])
            k == "K_POINTS" && read_kpoints(io[v])
            k == "CELL_PARAMETERS" && read_cellparameters(io[v])
            # TODO: Other cards
        end
    end  # for
    return Dict("namelists" => namelists, "cards" => cards)
end  # function dispatch_readers
function dispatch_readers(path::AbstractPath)
    isfile(path) && isreadable(path) || error("File $(path) not readable!")
    open(path, "r") do io
        dispatch_readers(io)
    end
end  # function dispatch_readers

@resumable function iterate_io_between(io::IOStream, start::Int, stop::Int)
    for i in eachindex(1:stop)
        if i == start
            @yield io
        else
            continue
        end  # if-else
    end  # for
end  # function iterate_io_between
