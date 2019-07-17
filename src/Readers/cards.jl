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

function read_atomicspecies(io::IOStream)
    atomic_species = []
    for line in io[2:end]
        # Skip the title line, any empty line, or a line of comment.
        isempty(line) || startswith(strip(line), '!') && continue
        m = match(r"(\S+)\s*(-?\d*\.?\d*)\s*(\S+)\s*", strip(line))
        if m === nothing
            @warn "No match found in the line $(line)!"
        else
            name, mass, pseudopotential = m.captures
            push!(atomic_species, AtomicSpecies(name, mass, pseudopotential))
        end
    end
    return AtomicSpeciesCard(data=atomic_species)
end  # function read_atomicspecies

function read_atomicpositions(io::IOStream)
    atomic_positions = []
    title_line = first(io)
    m = match(r"ATOMIC_POSITIONS\s*(?:[({])?\s*(\w*)\s*(?:[)}])?"i, title_line)
    m === nothing && error("No match found in the line '$(title_line)'! Something went wrong!")
    option = m.captures[2]
    if isempty(option)
        @warn "No option is found, default option 'alat' will be set! " *
              "Not specifying units is DEPRECATED and will no longer be allowed in the future"
        option = "alat"
    end
    for line in io[2:end]
        # If this line is an empty line or a line of comment.
        isempty(line) || startswith(strip(line), '!') && continue
        strip(line) == '/' && error("Do not start any line in cards with a '/' character!")
        if match(r"\\{.*\\}", line)
            m = match(r"(\w+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*\\{\s*([01])?\s*([01])?\s*([01])?\s*\\}",
                strip(line))
            name, x, y, z, if_pos1, if_pos2, if_pos3 = m.captures
            push!(atomic_positions, AtomicPosition(atom=name, position=[x, y, z]))
        else
            m = match(r"(\w+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)", strip(line))
            if m === nothing
                @warn "No match found in the line $(line)!"
            else
                name, x, y, z = m.captures
                push!(atomic_positions, AtomicPosition(atom=name, position=[x, y, z]))
            end
        end
    end
    return AtomicPositionCard(option=option, data=atomic_positions)
end  # function read_atomicpositions

function read_kpoints(io::IOStream)
    title_line = first(io)
    m = match(r"K_POINTS\s*(?:[({])?\s*(\w*)\s*(?:[)}])?"i, title_line)
    m === nothing && error("Match not found! Check your option!")
    option = match.captures[2]  # The second parenthesized subgroup will be `option`.

    isempty(option) && error("Option is not given! you must give one!")

    option == "gamma" && return KPointsCard(option=option, points=GammaPoint())

    if option == "automatic"
        for line in io[2:end]
            # If this line is an empty line or a line of comment.
            isempty(line) || startswith(strip(line), '!') && continue
            strip(line) == '/' && error("Do not start any line in cards with a '/' character!")
            line = split(line)
            grid, offsets = map(x -> parse(Int, x), line[1:3]), map(x -> parse(Int, x), line[4:7])
            return KPointsCard(option=option, points=MonkhorstPackGrid(grid=grid, offsets=offsets))
        end
    end

    option in ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c") && return nothing

    error("Unknown option '$option' given!")
end  # function read_kpoints

function read_cellparameters(io::IOStream)
    cell_params = []
    title_line = first(io)
    m = match(r"CELL_PARAMETERS\s*\\{?\s*(\w*)\s*\\}?"i, title_line)
    if match === nothing
        # The first line should be 'CELL_PARAMETERS blahblahblah', if it is not, either the regular expression
        # wrong or something worse happened.
        error("No match found! Check you 'CELL_PARAMETERS' line!")
    end
    option = match.captures[2]  # The second parenthesized subgroup will be `option`.
    if isempty(option)
        @warn "Not specifying unit is DEPRECATED and will no longer be allowed in the future!"
        option = "bohr"
    end
    for line in io[2:end]
        strip(line) == '/' && error("Do not start any line in cards with a '/' character!")
        if match(r"(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*", strip(line))
            v1, v2, v3 = match(r"(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*", strip(line)).captures
            push!(cell_params, [v1, v2, v3])
        end
    end
    return CellParametersCard(option, Crystal(cell_params))
end  # function read_cellparameters
