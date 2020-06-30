using Gaussians
using Test
using LinearAlgebra
using StaticArrays

@testset "Gaussians.jl" begin
    @testset "STO-NG" begin
        @testset "HeH⁺" begin
            ζ = [2.0925, 1.24]
            Z = [2.0, 1.0]
            @testset "$label" for (label,R) in [
                ("Scalar", [0.0, 1.4632]),
                ("Vector", [SA[0.0, 0.0, 0.0], SVector{3}(1.4632*normalize(rand(3)))])
            ]
                Φ = [sto_g(3, R[i], ζ[i]) for i = 1:2]

                S = Matrix(Φ', I, Φ)
                @test S ≈ [1 0.4507704116;0.4507704116 1]

                T = Matrix(Φ', KineticOperator(), Φ)
                @test T ≈ [2.164313 0.167013; 0.167013 0.760033] atol=1e-5

                V₁ = Matrix(Φ', CoulombPotential(-Z[1], R[1]), Φ)
                @test V₁ ≈ [-4.139827 -1.102912; -1.102912 -1.265246] atol=1e-5

                V₂ = Matrix(Φ', CoulombPotential(-Z[2], R[2]), Φ)
                @test V₂ ≈ [-0.677230 -0.411305; -0.411305 -1.226615] atol=1e-5

                eris = ElectronRepulsionIntegrals(Φ)
                display(eris)

                @test eris[1,1,1,1] ≈ 1.307152 atol=1e-5
                @test eris[2,1,1,1] ≈ 0.437279 atol=1e-5
                @test eris[2,1,2,1] ≈ 0.177267 atol=1e-5
                @test eris[2,2,1,1] ≈ 0.605703 atol=1e-5
                @test eris[2,2,2,1] ≈ 0.311795 atol=1e-5
                @test eris[2,2,2,2] ≈ 0.774608 atol=1e-5

                @test eris[2,1,1,1] == eris[1,2,1,1] == eris[1,1,2,1] == eris[1,1,1,2]
                @test eris[2,1,2,1] == eris[1,2,2,1] == eris[1,2,1,2] == eris[2,1,1,2]
                @test eris[2,2,1,1] == eris[1,1,2,2]
                @test eris[2,2,2,1] == eris[2,2,1,2] == eris[2,1,2,2] == eris[1,2,2,2]
            end
        end
    end
end
