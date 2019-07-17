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
        match = match(r"(\S+)\s*(-?\d*\.?\d*)\s*(\S+)\s*", strip(line))
        if match === nothing
            @warn "No match found in the line $(line)!"
        else
            name, mass, pseudopotential = match.captures
            push!(atomic_species, AtomicSpecies(name, mass, pseudopotential))
        end
    end
    return AtomicSpeciesCard(data = atomic_species)
end  # function read_atomicspecies

end
