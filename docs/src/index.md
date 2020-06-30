```@meta
CurrentModule = Gaussians
```

# Gaussians

## Simple usage

This is the HeH⁺ example of Szabo & Ostlund (1996), where the two
constituent atoms are modelled by one STO-3G orbital each.

```julia-repl
julia> using Gaussians

julia> using LinearAlgebra

julia> R = [0.0, 1.4632] # Positions of the nuclei
2-element Array{Float64,1}:
 0.0
 1.4632

julia> ζ = [2.0925, 1.24] # Exponents
2-element Array{Float64,1}:
 2.0925
 1.24

julia> Z = [2.0, 1.0] # Nuclear charges
2-element Array{Float64,1}:
 2.0
 1.0

julia> Φ = [sto_g(3, R[i], ζ[i]) for i = 1:2]
2-element Array{STO_NG{Float64,Float64},1}:
 STO-3G basis function: [0.18298567211025116×1s(0.48084429026249986), 0.5871364000040699×1s(1.7766911481187495), 0.6070752318149509×1s(9.753934615874998)]
 STO-3G basis function: [0.08347403618820594×1s(0.16885615680000002), 0.2678387030861255×1s(0.6239134896), 0.27693435931394883×1s(3.4252500160000006)]
```

With the basis defined, we can compute the various integrals:

```julia-repl
julia> S = Matrix(Φ', I, Φ)
2×2 Array{Float64,2}:
 1.0      0.45077
 0.45077  1.0

julia> T = Matrix(Φ', KineticOperator(), Φ)
2×2 Array{Float64,2}:
 2.16431   0.167013
 0.167013  0.760033

julia> V₁ = Matrix(Φ', CoulombPotential(-Z[1], R[1]), Φ)
2×2 Array{Float64,2}:
 -4.13983  -1.10291
 -1.10291  -1.26525

julia> V₂ = Matrix(Φ', CoulombPotential(-Z[2], R[2]), Φ)
2×2 Array{Float64,2}:
 -0.67723   -0.411305
 -0.411305  -1.22662

julia> eris = ElectronRepulsionIntegrals(Φ)
2 ElectronRepulsionIntegrals stored in a Array{Float64,4} with 6 unique integrals:
#1     1   1   1   1: 1.307152e+00
#2     1   2   1   1: 4.372793e-01
#3     1   2   1   2: 1.772671e-01
#4     2   2   1   1: 6.057034e-01
#5     2   2   1   2: 3.117946e-01
#6     2   2   2   2: 7.746084e-01
```

## Reference

Equations tagged (SO.X.Y) refer to

    - Szabo, A., & Ostlund, N. S. (1996). Modern Quantum Chemistry:
      Introduction to Advanced Electronic Structure Theory (Dover Books on
      Chemistry). Dover Publications.

```@index
```

```@autodocs
Modules = [Gaussians]
```
