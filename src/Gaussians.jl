module Gaussians

using LinearAlgebra
using SpecialFunctions

using Formatting

include("sto_ng.jl")

export STO_NG, sto_g,
    KineticOperator, CoulombPotential,
    ElectronRepulsionIntegrals

end
