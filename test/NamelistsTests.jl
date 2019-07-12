#=
NamelistsTests.jl:
- Julia version: 1.0
- Author: qz
- Date: Jul 11, 2019
=#
module NamelistsTests

using Test

using QuantumESPRESSO.Namelists.PWscf
using QuantumESPRESSO.Cards.PWscf
using QuantumESPRESSO.QuantumESPRESSOInput.PWscf

as = AtomicSpeciesCard(data=[AtomicSpecies("Fe", 55.845, "Fe.pseudopotential")])
ap = AtomicPositionCard(data=[AtomicPosition("Fe", [0, 0, 0])])
cell = CellParametersCard()
k = KPointsCard(points=[GammaPoint()])

pw = PWscfInput(system=SystemNamelist(celldm=ones(6)), atomicspecies=as, atomicpositions=ap,
    kpoints=k, cellparameters=cell
)

end
