#=
cards:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-17
=#
using Crystals

using QuantumESPRESSO.Cards.PWscf

export read_atomicspecies,
    read_atomicpositions,
    read_kpoints,
    read_cellparameters

function read_title_line(title_line, regex, default_option)
    m = match(regex, title_line)
    if m === nothing
        # The first line should be '<CARD> {<option>}', if it is not, either the regular expression
        # wrong or something worse happened.
        error("No match found in (title_line)!")
    else
        option = m.captures[1]  # The first parenthesized subgroup will be `option`.
    end
    if isempty(option)
        @warn "No option is found, default option '(default_option)' will be set!"
        option = default_option
    end
    return option
end  # function read_title_line

function preprocess_line(line)
    str = strip(line)
    # If this line is an empty line or a line of comment.
    # Comments lines in cards can be introduced by either a "!" or a "#" character in the first position of a line.
    isempty(str) || any(startswith(str, x) for x in ('!', '#')) && return nothing
    # Do not start any line in cards with a "/" character.
    str == '/' && error("Do not start any line in cards with a '/' character!")
    return str
end  # function preprocess_line

function read_atomicspecies(lines)
    atomic_species = []
    for line in lines
        str = strip(line)
        # Skip the title line, any empty line, or a line of comment.
        (isempty(str) || startswith(strip(str), '!') || occursin(r"ATOMIC_SPECIES"i, str)) && continue
        m = match(r"(\S+)\s*(-?\d*\.?\d*)\s*(\S+)\s*", str)
        if m === nothing
            @warn "No match found in the line $(line)!"
        else
            name, mass, pseudopotential = m.captures
            push!(atomic_species, AtomicSpecies(name, parse(Float64, mass), pseudopotential))
        end
    end
    return AtomicSpeciesCard(data=atomic_species)
end  # function read_atomicspecies

function read_atomicpositions(lines)
    atomic_positions = []
    option = read_title_line(first(lines), r"ATOMIC_POSITIONS\s*(?:[({])?\s*(\w*)\s*(?:[)}])?"i, "alat")
    for line in Iterators.drop(lines, 1)  # Drop the title line
        str = preprocess_line(line)
        str === nothing && continue

        if match(r"\{.*\}", str) !== nothing
            m = match(r"(\w+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*\{\s*([01])?\s*([01])?\s*([01])?\s*\}", str)
            name, x, y, z, if_pos1, if_pos2, if_pos3 = m.captures
            push!(atomic_positions, AtomicPosition(name, map(x->parse(Float64, x), [x, y, z])))
        else
            m = match(r"(\w+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)", str)
            if m === nothing
                @warn "No match found in the line $(line)!"
            else
                name, x, y, z = m.captures
                push!(atomic_positions, AtomicPosition(name, map(x->parse(Float64, x), [x, y, z])))
            end
        end
    end
    return AtomicPositionCard(option=option, data=atomic_positions)
end  # function read_atomicpositions

function read_kpoints(lines)
    option = read_title_line(first(lines), r"K_POINTS\s*(?:[({])?\s*(\w*)\s*(?:[)}])?"i, "tbipa")

    option == "gamma" && return KPointsCard(option=option, points=GammaPoint())

    if option == "automatic"
        for line in Iterators.drop(lines, 1)  # Drop the title line
            str = preprocess_line(line)

            sp = split(str)
            grid, offsets = map(x -> parse(Int, x), sp[1:3]), map(x -> parse(Int, x), sp[4:6])
            return KPointsCard(option=option, points=[MonkhorstPackGrid(grid=grid, offsets=offsets)])
        end
    end

    option in ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c") && return nothing

    error("Unknown option '$option' given!")
end  # function read_kpoints

function read_cellparameters(lines)
    cell_params = []
    option = read_title_line(first(lines), r"CELL_PARAMETERS\s*[\{\(]?\s*(\w*)\s*[\}\)]?"i, "bohr")

    for line in Iterators.drop(lines, 1)  # Drop the title line
        str = preprocess_line(line)
        str === nothing && continue

        m = match(r"(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*", str)
        if m !== nothing
            v1, v2, v3 = m.captures
            cell_params = vcat(cell_params, map(x->parse(Float64, x), [v1, v2, v3]))
        end
    end
    return CellParametersCard(option, reshape(cell_params, (3, 3)))
end  # function read_cellparameters
