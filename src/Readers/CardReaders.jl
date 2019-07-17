"""
# module CardReaders



# Examples

```jldoctest
julia>
```
"""
module CardReaders

using QuantumESPRESSO.Cards.PWscf

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
    return AtomicSpeciesCard(data = atomic_species)
end  # function read_atomicspecies

function read_atomicpositions(io::IOStream)
    atomic_positions = []
    title_line = io[1]
    m = match(r"ATOMIC_POSITIONS\s*(?:[({])?\s*(\w*)\s*(?:[)}])?", title_line, flags = re.IGNORECASE)
    m === nothing && error("No match found in the line '$(title_line)'! Something went wrong!")
    option = m.captures[2]
    if isempty(option)
        @warn "No option is found, default option 'alat' will be set! " *
              "Not specifying units is DEPRECATED and will no longer be allowed in the future"
        option = "alat"
    end
    for line in s[2:end]
        # If this line is an empty line or a line of comment.
        isempty(line) || startswith(strip(line), '!') && continue
        if strip(line) == '/'
            error("Do not start any line in cards with a '/' character!")
        end
        if match(r"\\{.*\\}", line):
            m = match(r"(\w+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*\\{\s*([01])?\s*([01])?\s*([01])?\s*\\}",
                strip(line))
            name, x, y, z, if_pos1, if_pos2, if_pos3 = m.captures
            push!(atomic_positions, AtomicPosition(atom = name, position = [x, y, z]))
        else
            m = match("(\w+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)", strip(line))
            if m === nothing
                @warn "No match found in the line $(line)!"
            else
                name, x, y, z = m.captures
                push!(atomic_positions, AtomicPosition(atom = name, position = [x, y, z]))
            end
        end
    end
    return AtomicPositionCard(option = option, data = atomic_positions)
end  # function read_atomicpositions

end
