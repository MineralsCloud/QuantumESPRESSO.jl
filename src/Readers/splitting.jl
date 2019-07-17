#=
splitting:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-07-17
=#
const NAMELIST_IDENTIFIERS = "&CONTROL", "&SYSTEM", "&ELECTRONS", "&IONS", "&CELL"
const CARD_IDENTIFIERS = "ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS", "CELL_PARAMETERS", "OCCUPATIONS", "CONSTRAINTS", "ATOMIC_FORCES"

struct RangeIndices
    startindex
    endindex
end  # struct RangeIndices

function get_namelist_identifier_positions(io::IOStream)
    match_records = Dict()
    identifiers = NAMELIST_IDENTIFIERS
    for pattern in identifiers
        # ``re.compile`` will produce a regular expression object, on which we can use its ``search`` method.
        m0 = re.compile(pattern, flags=re.IGNORECASE).search(io, match_records.get(pattern, pos))
        m0 === nothing && continue
        m1 = re.compile(self.namelist_sep, flags=re.IGNORECASE).search(io, m0.end())
        match_records[pattern] = RangeIndices(startindex={True: m0.start(), False: m0.end()}[include_heading], endindex={True: m1.end(), False: m1.start()}[include_ending])
    end
    # The following one-line code first sort the ``dict`` *positions* by its values, i.e., a ``RangeIndices`` tuple.
    # Then we get a list of tuples, with first entries to be the identifiers and second indices to be
    # the indices.
    # For example, if ``x = {'a': 2, 'b': 3, 'c': 1}``, then
    # >>> sorted(x.items(), key=operator.itemgetter(1))
    # [('c', 1), ('a', 2), ('b', 3)]
    # Tuples are compared lexicographically using comparison of corresponding elements, thus compared with their
    # *begin* entry.
    m = sorted(match_records.items(), key=operator.itemgetter(1))
    return OrderedDict(m)
end  # function get_namelist_identifier_positions
