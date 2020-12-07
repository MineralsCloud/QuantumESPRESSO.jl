module Outputs

using QuantumESPRESSOParsers.Outputs: SubroutineError

export SubroutineError

module PWscf

using QuantumESPRESSOParsers.Outputs.PWscf:
    Diagonalization,
    Preamble,
    Diagonalization,
    Davidson,
    ConjugateGradient,
    ProjectedPreconditionedConjugateGradient,
    parse_fft_base_info,
    parse_symmetries,
    parse_ibz,
    parse_stress,
    parse_iteration_time,
    parse_bands,
    parse_all_electron_energy,
    parse_energy_decomposition,
    parse_paw_contribution,
    parse_smearing_energy,
    parse_version,
    parse_parallel_info,
    parse_fft_dimensions,
    parse_iteration_head,
    parse_electrons_energies,
    parse_clock,
    # whatinput,
    isoptimized,
    isjobdone,
    tryparsefirst,
    parsefirst,
    tryparseall,
    parseall,
    tryparselast,
    parselast,
    tryparsenext,
    parsenext,
    tryparsefinal,
    parsefinal

export Diagonalization,
    Preamble,
    Davidson,
    ConjugateGradient,
    ProjectedPreconditionedConjugateGradient,
    parse_fft_base_info,
    parse_symmetries,
    parse_ibz,
    parse_stress,
    parse_iteration_time,
    parse_bands,
    parse_all_electron_energy,
    parse_energy_decomposition,
    parse_paw_contribution,
    parse_smearing_energy,
    parse_version,
    parse_parallel_info,
    parse_fft_dimensions,
    parse_iteration_head,
    parse_electrons_energies,
    parse_clock,
    # whatinput,
    isoptimized,
    isjobdone,
    tryparsefirst,
    parsefirst,
    tryparseall,
    parseall,
    tryparselast,
    parselast,
    tryparsenext,
    parsenext,
    tryparsefinal,
    parsefinal

end # module PWscf

module PHonon

using QuantumESPRESSOParsers.Outputs.PHonon: parse_frequency, parse_dos

export parse_frequency, parse_dos

end # module PHonon

end # module Outputs
