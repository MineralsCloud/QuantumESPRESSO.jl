module Inputs

using QuantumESPRESSOBase.Inputs:
    Namelist,
    Card,
    InputEntry,
    getoption,
    allowed_options,
    titleof,
    inputstring,
    compulsory_namelists,
    optional_namelists,
    compulsory_cards,
    optional_cards,
    allnamelists,
    allcards
using QuantumESPRESSOParsers.Inputs: InvalidInput, InputFile

export InvalidInput,
    InputFile, getnamelists, getcards, getoption, allowed_options, titleof, inputstring

module PWscf

using QuantumESPRESSOBase.Inputs.PWscf:
    ControlNamelist,
    SystemNamelist,
    ElectronsNamelist,
    IonsNamelist,
    CellNamelist,
    DosNamelist,
    BandsNamelist,
    AtomicSpecies,
    AtomicSpeciesCard,
    AtomicPosition,
    AtomicPositionsCard,
    CellParametersCard,
    AtomicForce,
    AtomicForcesCard,
    MonkhorstPackGrid,
    GammaPoint,
    SpecialKPoint,
    KPointsCard,
    PWInput,
    optconvert,
    xmldir,
    wfcfiles,
    getoption,
    allowed_options,
    allnamelists,
    allcards,
    compulsory_namelists,
    optional_namelists,
    compulsory_cards,
    optional_cards
using QuantumESPRESSOParsers.Inputs.PWscf

export ControlNamelist,
    SystemNamelist,
    ElectronsNamelist,
    IonsNamelist,
    CellNamelist,
    DosNamelist,
    BandsNamelist,
    AtomicSpecies,
    AtomicSpeciesCard,
    AtomicPosition,
    AtomicPositionsCard,
    CellParametersCard,
    AtomicForce,
    AtomicForcesCard,
    MonkhorstPackGrid,
    GammaPoint,
    SpecialKPoint,
    KPointsCard,
    PWInput,
    optconvert,
    xmldir,
    wfcfiles,
    getoption,
    allowed_options,
    allnamelists,
    allcards,
    compulsory_namelists,
    optional_namelists,
    compulsory_cards,
    optional_cards

end # module PWscf

module PHonon

using QuantumESPRESSOBase.Inputs.PHonon:
    PhNamelist,
    Q2rNamelist,
    MatdynNamelist,
    DynmatNamelist,
    Q2rInput,
    DynmatInput,
    PhInput,
    MatdynInput

export PhNamelist,
    Q2rNamelist, MatdynNamelist, DynmatNamelist, Q2rInput, DynmatInput, PhInput, MatdynInput

end # module PHonon

end # module Inputs
