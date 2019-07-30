"""
# module PW



# Examples

```jldoctest
julia>
```
"""
module PW

using Compat: isnothing
using Fortran90Namelists.FortranToJulia: FortranData

using QuantumESPRESSOBase.Cards.PW

export read_atomicspecies,
    read_atomicpositions,
    read_kpoints,
    read_cellparameters

function read_title_line(title_line, regex, default_option)
    m = match(regex, title_line)
    if isnothing(m)
        # The first line should be '<CARD> {<option>}', if it is not, either the regular expression
        # wrong or something worse happened.
        error("No match found in $(title_line)!")
    else
        option = m.captures[1]  # The first parenthesized subgroup will be `option`.
    end
    if isempty(option)
        @warn "No option is found, default option '$(default_option)' will be set!"
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
    atomic_species = AtomicSpecies[]
    for line in Iterators.drop(lines, 1)  # Drop the title line
        str = preprocess_line(line)
        isnothing(str) && continue

        m = match(r"(\S+)\s*(-?\d*\.?\d*)\s*(\S+)\s*", str)
        if isnothing(m)
            @warn "No match found in the line $(line)!"
        else
            atom, mass, pseudopotential = m.captures
            push!(atomic_species, AtomicSpecies(string(atom), parse(Float64, FortranData(mass)), string(pseudopotential)))
        end
    end
    return AtomicSpeciesCard(atomic_species)
end  # function read_atomicspecies

function read_atomicpositions(lines)
    atomic_positions = AtomicPosition[]
    option = read_title_line(first(lines), r"ATOMIC_POSITIONS\s*(?:[({])?\s*(\w*)\s*(?:[)}])?"i, "alat")
    for line in Iterators.drop(lines, 1)  # Drop the title line
        str = preprocess_line(line)
        isnothing(str) && continue

        if !isnothing(match(r"\{.*\}", str))
            m = match(r"(\w+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*\{\s*([01])?\s*([01])?\s*([01])?\s*\}", str)
            atom, x, y, z, if_pos1, if_pos2, if_pos3 = m.captures
            push!(atomic_positions, AtomicPosition(atom = string(atom),
                                                   pos = [parse(Float64, FortranData(p)) for p in (x, y, z)],
                                                   if_pos = [parse(Float64, FortranData(x)) for x in (if_pos1, if_pos2, if_pos3)]))
        else
            m = match(r"(\w+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)", str)
            if isnothing(m)
                @warn "No match found in the line $(line)!"
            else
                atom, x, y, z = m.captures
                push!(atomic_positions, AtomicPosition(atom = string(atom),
                pos = [parse(Float64, FortranData(p)) for p in (x, y, z)]))
            end
        end
    end
    return AtomicPositionsCard(option = string(option), data = atomic_positions)
end  # function read_atomicpositions

function read_kpoints(lines)
    option = read_title_line(first(lines), r"K_POINTS\s*(?:[({])?\s*(\w*)\s*(?:[)}])?"i, "tbipa")

    option == "gamma" && return KPointsCard(option = string(option), points = GammaPoint())

    if option == "automatic"
        for line in Iterators.drop(lines, 1)  # Drop the title line
            str = preprocess_line(line)
            isnothing(str) && continue

            sp = split(str)
            grid, offsets = [parse(Int, FortranData(x)) for x in sp[1:3]], [parse(Int, FortranData(x)) for x in sp[4:6]]
            return KPointsCard(option = string(option), data = [MonkhorstPackGrid(grid = grid, offsets = offsets)])
        end
    end

    if option in ("tpiba", "crystal", "tpiba_b", "crystal_b", "tpiba_c", "crystal_c")
        kpoints = SpecialKPoint[]
        for line in Iterators.drop(lines, 1)  # Drop the title line
            str = preprocess_line(line)
            isnothing(str) && continue

            sp = split(str)
            length(sp) == 1 && (nks = parse(Int, FortranData(first(sp))))
            (@isdefined nks) && break
        end
        for line in lines
            str = preprocess_line(line)
            isnothing(str) && continue

            sp = split(str)
            length ≠ 4 && error("Unknown input given!")
            push!(kpoints, SpecialKPoint(collect(parse(Float64, FortranData(x)) for x in sp[1:3]), parse(Float64, FortranData(sp[4]))))
        end
        length(kpoints) ≠ nks && throw(DimensionMismatch("The length of k-points $(length(kpoints)) is not equal to $(nks)!"))
        return KPointsCard(option = string(option), data = kpoints)
    end

    error("Unknown option '$option' given!")
end  # function read_kpoints

function read_cellparameters(lines)
    cell_params = []
    option = read_title_line(first(lines), r"CELL_PARAMETERS\s*[\{\(]?\s*(\w*)\s*[\}\)]?"i, "bohr")

    for line in Iterators.drop(lines, 1)  # Drop the title line
        str = preprocess_line(line)
        isnothing(str) && continue

        m = match(r"(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*(-?\d*\.\d*)\s*", str)
        if !isnothing(m)
            v1, v2, v3 = m.captures
            cell_params = vcat(cell_params, [parse(Float64, FortranData(x)) for x in (v1, v2, v3)])
        end
    end
    return CellParametersCard(string(option), reshape(cell_params, (3, 3)))
end  # function read_cellparameters

end
