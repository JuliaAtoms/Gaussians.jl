var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = Gaussians","category":"page"},{"location":"#Gaussians","page":"Home","title":"Gaussians","text":"","category":"section"},{"location":"#Simple-usage","page":"Home","title":"Simple usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This is the HeH⁺ example of Szabo & Ostlund (1996), where the two constituent atoms are modelled by one STO-3G orbital each.","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Gaussians\n\njulia> using LinearAlgebra\n\njulia> R = [0.0, 1.4632] # Positions of the nuclei\n2-element Array{Float64,1}:\n 0.0\n 1.4632\n\njulia> ζ = [2.0925, 1.24] # Exponents\n2-element Array{Float64,1}:\n 2.0925\n 1.24\n\njulia> Z = [2.0, 1.0] # Nuclear charges\n2-element Array{Float64,1}:\n 2.0\n 1.0\n\njulia> Φ = [sto_g(3, R[i], ζ[i]) for i = 1:2]\n2-element Array{STO_NG{Float64,Float64},1}:\n STO-3G basis function: [0.18298567211025116×1s(0.48084429026249986), 0.5871364000040699×1s(1.7766911481187495), 0.6070752318149509×1s(9.753934615874998)]\n STO-3G basis function: [0.08347403618820594×1s(0.16885615680000002), 0.2678387030861255×1s(0.6239134896), 0.27693435931394883×1s(3.4252500160000006)]","category":"page"},{"location":"","page":"Home","title":"Home","text":"With the basis defined, we can compute the various integrals:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> S = Matrix(Φ', I, Φ)\n2×2 Array{Float64,2}:\n 1.0      0.45077\n 0.45077  1.0\n\njulia> T = Matrix(Φ', KineticOperator(), Φ)\n2×2 Array{Float64,2}:\n 2.16431   0.167013\n 0.167013  0.760033\n\njulia> V₁ = Matrix(Φ', CoulombPotential(-Z[1], R[1]), Φ)\n2×2 Array{Float64,2}:\n -4.13983  -1.10291\n -1.10291  -1.26525\n\njulia> V₂ = Matrix(Φ', CoulombPotential(-Z[2], R[2]), Φ)\n2×2 Array{Float64,2}:\n -0.67723   -0.411305\n -0.411305  -1.22662\n\njulia> eris = ElectronRepulsionIntegrals(Φ)\n2 ElectronRepulsionIntegrals stored in a Array{Float64,4} with 6 unique integrals:\n#1     1   1   1   1: 1.307152e+00\n#2     1   2   1   1: 4.372793e-01\n#3     1   2   1   2: 1.772671e-01\n#4     2   2   1   1: 6.057034e-01\n#5     2   2   1   2: 3.117946e-01\n#6     2   2   2   2: 7.746084e-01","category":"page"},{"location":"#Reference","page":"Home","title":"Reference","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Equations tagged (SO.X.Y) refer to","category":"page"},{"location":"","page":"Home","title":"Home","text":"- Szabo, A., & Ostlund, N. S. (1996). Modern Quantum Chemistry:\n  Introduction to Advanced Electronic Structure Theory (Dover Books on\n  Chemistry). Dover Publications.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Gaussians]","category":"page"},{"location":"#Gaussians.STO_NG","page":"Home","title":"Gaussians.STO_NG","text":"STO_NG(n, ℓ, α, c, R)\n\nSlater-type orbital constructed from N primitive Gaussian-type orbitals with principal quantum numbers n, angular momenta ℓ, exponents α, and contraction coefficients c, all of which are centred at R (which may be a scalar for linear molecules or an AbstractVector with 3 elements for arbitrary molecules).\n\nSo far, only 1s primitive Gaussians, which are written as\n\nbeginequation\ntagSOB1\ng_mathrm1s(alpha) = (2alphapi)^34exp(-alpha r^2)\nendequation\n\nare supported.\n\njulia> STO_NG([1], [0], [0.5], [1.0])\nSTO-1G basis function:\n    0.423777 × 1s(  0.500000)\n\n\njulia> STO_NG([1,1], [0,0], [0.5,0.25], [1,1]/√2)\nSTO-2G basis function:\n    0.299656 × 1s(  0.500000)\n    0.178176 × 1s(  0.250000)\n\n\n\n\n\n","category":"type"},{"location":"#Gaussians.sto_g","page":"Home","title":"Gaussians.sto_g","text":"sto_g(N, R[, ζ=1.0])\n\nStandard contractions for N=123, centred at R, and with exponent ζ.\n\n\n\n\n\n","category":"function"}]
}
