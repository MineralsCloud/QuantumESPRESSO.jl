#=
NamelistsTests.jl:
- Julia version: 1.0
- Author: qz
- Date: Jul 11, 2019
=#
module NamelistsTests

using LinearAlgebra
using Test

using QuantumESPRESSO.Namelists.PW
using QuantumESPRESSO.Cards.PW
using QuantumESPRESSO.QuantumESPRESSOInput.PW

as = AtomicSpeciesCard(data=[AtomicSpecies("Fe", 55.845, "Fe.pseudopotential")])
ap = AtomicPositionCard(data=[AtomicPosition(atom="Fe", pos=[0, 0, 0])])
cell = CellParametersCard(data=diagm(0=>[1, 1, 1]))
k = KPointsCard(data=[GammaPoint()])

pw = PWInput(system=SystemNamelist(celldm=ones(6)), atomicspecies=as, atomicpositions=ap,
    kpoints=k, cellparameters=cell
)

end
