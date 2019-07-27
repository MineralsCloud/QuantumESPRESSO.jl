#=
splitting:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-17
=#
using DataStructures
using FilePaths: AbstractPath
using Parameters

using QuantumESPRESSO: name
using QuantumESPRESSO.Readers.Namelists
using QuantumESPRESSO.Readers.Cards.PW
using QuantumESPRESSO.QuantumESPRESSOInput.PW

export namelist_identifier_linenumbers,
    namelist_lineranges,
    card_identifier_linenumbers,
    card_lineranges,
    input_identifier_linenumbers,
    input_lineranges,
    dispatch_readers,
    form_input_object

const NAMELIST_END = '/'  # Not a regex anymore, since I strip everyline
const NAMELIST_STARTS = "&CONTROL", "&SYSTEM", "&ELECTRONS", "&IONS", "&CELL"  # regex: "&(.[^,]*)"
const CARD_STARTS = "ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS", "CELL_PARAMETERS", "OCCUPATIONS", "CONSTRAINTS", "ATOMIC_FORCES"

macro iostream_to_lines(methodname)
    return quote
        function $(esc(methodname))(io::IOStream)
            $(esc(methodname))(readlines(io))
        end
    end
end  # macro iostream_to_lines

macro path_to_iostream(methodname)
    return quote
        function $(esc(methodname))(path::AbstractPath)
            isfile(path) && isreadable(path) || error("File $(path) not readable!")
            open(path, "r") do io
                $(esc(methodname))(io)
            end
        end
    end
end  # macro path_to_iostream

function namelist_identifier_linenumbers(lines)
    records = OrderedDict()
    for (i, line) in enumerate(lines)
        str = strip(line)
        isempty(str) || startswith(str, '!') || startswith(str, '#') && continue
        for namelistname in NAMELIST_STARTS
            occursin(Regex("$namelistname", "i"), str) ? records[namelistname] = i : continue
        end  # for
    end  # for
    return records
end  # function namelist_identifier_linenumbers
@iostream_to_lines namelist_identifier_linenumbers
@path_to_iostream namelist_identifier_linenumbers

function namelist_lineranges(lines)
    records = OrderedDict()
    for (i, line) in enumerate(lines)
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
    if all(isincreasing(x) for x in values(records))
        return records
    else
        error("Something went wrong!")
    end  # if-else
end  # function namelist_lineranges
@iostream_to_lines namelist_lineranges
@path_to_iostream namelist_lineranges

function card_identifier_linenumbers(lines)
    records = OrderedDict()
    for (i, line) in enumerate(lines)
        str = strip(line)
        isempty(str) || startswith(str, '!') || startswith(str, '#') && continue
        for cardname in CARD_STARTS
            occursin(Regex("$cardname", "i"), str) ? records[cardname] = i : continue
        end  # for
    end  # for
    if haskey(records, "OCCUPATIONS")
        linenumber = last(collect(values(namelist_identifier_linenumbers(lines))))
        records["OCCUPATIONS"] < linenumber && pop!(records, "OCCUPATIONS")
    end  # if
    return records
end  # function card_identifier_linenumbers
@iostream_to_lines card_identifier_linenumbers
@path_to_iostream card_identifier_linenumbers

function card_lineranges(lines)
    records = OrderedDict()
    for (i, line) in enumerate(lines)
        str = strip(line)
        isempty(str) || startswith(str, '!') || startswith(str, '#') && continue
        for cardname in CARD_STARTS
            occursin(Regex("$cardname", "i"), str) ? records[cardname] = i : continue
        end  # for
    end  # for
    if haskey(records, "OCCUPATIONS")
        linenumber = last(collect(values(namelist_identifier_linenumbers(lines))))
        records["OCCUPATIONS"] < linenumber && pop!(records, "OCCUPATIONS")
    end  # if
    for (i, (k, v)) in enumerate(records)
        if i == length(values(records))
            records[k] = v:length(lines)
        else
            nextkey = collect(keys(records))[i + 1]
            records[k] = v:(records[nextkey] - 1)
        end  # if
    end  # for
    if all(isincreasing(x) for x in values(records))
        return records
    else
        error("Something went wrong!")
    end  # if-else
end  # function card_lineranges
@iostream_to_lines card_lineranges
@path_to_iostream card_identifier_linenumbers

function input_identifier_linenumbers(lines)
    Dict("namelists" => namelist_identifier_linenumbers(lines), "cards" => card_identifier_linenumbers(lines))
end  # function input_identifier_linenumbers
@iostream_to_lines input_identifier_linenumbers
@path_to_iostream input_identifier_linenumbers

function input_lineranges(lines)
    Dict("namelists" => namelist_lineranges(lines), "cards" => card_lineranges(lines))
end  # function input_lineranges
@iostream_to_lines input_lineranges
@path_to_iostream input_lineranges

function dispatch_readers(lines)
    lineranges = input_lineranges(lines)
    namelist_lineranges = lineranges["namelists"]
    card_lineranges = lineranges["cards"]
    namelists = Dict()
    cards = Dict()
    for (k, v) in namelist_lineranges
        namelists[k] = read_namelist(iterate_lines_between(lines, v))
    end  # for
    for (k, v) in card_lineranges
        cards[k] = begin
            if k == "ATOMIC_SPECIES"
                read_atomicspecies(iterate_lines_between(lines, v))
            elseif k == "ATOMIC_POSITIONS"
                read_atomicpositions(iterate_lines_between(lines, v))
            elseif k == "K_POINTS"
                read_kpoints(iterate_lines_between(lines, v))
            elseif k == "CELL_PARAMETERS"
                read_cellparameters(iterate_lines_between(lines, v))
            else
                error("Unrecognized card name $(k)!")
            end  # if-elseif-else
            # TODO: Other cards
        end
    end  # for
    return Dict("namelists" => namelists, "cards" => cards)
end  # function dispatch_readers
@iostream_to_lines dispatch_readers
@path_to_iostream dispatch_readers

function form_input_object(lines)
    dict = dispatch_readers(lines)
    d = Dict()
    for v in values(dict["namelists"])
        d[name(typeof(v))] = v
    end  # for
    for v in values(dict["cards"])
        d[name(typeof(v))] = v
    end  #
    return PWInput(; d...)
end  # function form_input_object

isincreasing(r::UnitRange) = r.stop > r.start ? true : false

function iterate_lines_between(lines, start::Int, stop::Int)
    Iterators.take(Iterators.drop(lines, start - 1), stop - start + 1)
end  # function iterate_lines_between
iterate_lines_between(lines, r::UnitRange) = iterate_lines_between(lines, r.start, r.stop)
