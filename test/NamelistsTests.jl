#=
NamelistsTests.jl:
- Julia version: 1.0
- Author: qz
- Date: Jul 11, 2019
=#
module NamelistsTests

using Test

using QuantumESPRESSO.Namelists.PW
using QuantumESPRESSO.Cards.PW
using QuantumESPRESSO.QuantumESPRESSOInput.PW

as = AtomicSpeciesCard(data=[AtomicSpecies("Fe", 55.845, "Fe.pseudopotential")])
ap = AtomicPositionCard(data=[AtomicPosition("Fe", [0, 0, 0])])
cell = CellParametersCard()
k = KPointsCard(points=[GammaPoint()])

pw = PWInput(system=SystemNamelist(celldm=ones(6)), atomicspecies=as, atomicpositions=ap,
    kpoints=k, cellparameters=cell
)

end
