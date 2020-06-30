# * STO-NG

@doc raw"""
    STO_NG(n, ℓ, α, c, R)

Slater-type orbital constructed from `N` primitive Gaussian-type
orbitals with principal quantum numbers `n`, angular momenta `ℓ`,
exponents `α`, and contraction coefficients `c`, all of which are
centred at `R` (which may be a scalar for linear molecules or an
`AbstractVector` with 3 elements for arbitrary molecules).

So far, only 1s primitive Gaussians, which are written as
```math
\begin{equation}
\tag{SO.B.1}
g_{\mathrm{1s}}(\alpha) = (2\alpha/\pi)^{3/4}\exp(-\alpha r^2)
\end{equation}
```
are supported.

```jldoctest
julia> STO_NG([1], [0], [0.5], [1.0])
STO-1G basis function:
    0.423777 × 1s(  0.500000)


julia> STO_NG([1,1], [0,0], [0.5,0.25], [1,1]/√2)
STO-2G basis function:
    0.299656 × 1s(  0.500000)
    0.178176 × 1s(  0.250000)
```
"""
struct STO_NG{T,C}
    n::Vector{Int}
    ℓ::Vector{Int}
    α::Vector{T}
    c::Vector{T}
    R::C # If C is a scalar type, then we're considering linear
    # molecules, else arbitrary
    function STO_NG(n, ℓ, α::Vector{T}, c::Vector{T}, R::C=zero(T), ζ=one(T)) where {T,C}
        all(isone, n) && all(iszero, ℓ) ||
            throw(ArgumentError("Only 1s Gaussians support so far"))
        N = length(α)
        N == length(n) == length(ℓ) == length(c) ||
            throw(DimensionMismatch())

        α = α .* ζ^2
        c = c .* (2α/π).^(3/4)

        new{T,C}(n, ℓ, α, c, R)
    end
end
Base.eltype(g::STO_NG{T}) where T = T
Base.adjoint(g::STO_NG) = Adjoint(g)
Base.zip(g::STO_NG) = zip(g.α,g.c)

const spectroscopic_labels = "spdfghiklmnoqrtuvwxyz"

function Base.show(io::IO, g::STO_NG)
    N = length(g.α)
    write(io, "STO-$(N)G basis function: [")
    write(io, join(["$(g.c[i])×$(g.n[i])$(spectroscopic_labels[g.ℓ[i]+1])($(g.α[i]))"
                    for i = 1:N], ", "))
    write(io, "]")
end

function Base.show(io::IO, ::MIME"text/plain", g::STO_NG)
    N = length(g.α)
    write(io, "STO-$(N)G basis function:\n")
    for i ∈ 1:N
        printfmtln(io, "{1:12.6f} ×{2:>2s}{3:s}({4:10.6f})",
                   g.c[i], g.n[i], spectroscopic_labels[g.ℓ[i]+1], g.α[i])
    end
end

function Base.show(io::IO, gc::Adjoint{<:Any,<:STO_NG})
    write(io, "Adjoint ")
    show(io, parent(gc))
end

function Base.show(io::IO, mime::MIME"text/plain", gc::Adjoint{<:Any,<:STO_NG})
    write(io, "Adjoint ")
    show(io, mime, parent(gc))
end

# # Docstring does not work with Documenter.jl due to https://github.com/JuliaDocs/Documenter.jl/issues/1192
# @doc raw"""
#     (g::STO_NG)(r)

# Evaluate the GTO `g` at ``r`` contracting over the primitive 1s Gaussians

# ```math
# \begin{equation}
# \tag{SO.B.1}
# g_{\mathrm{1s}}(\alpha) = (2\alpha/\pi)^{3/4}\exp(-\alpha r^2)
# \end{equation}
# ```
# """
function (g::STO_NG{T})(r) where T
    v = zero(T)
    rR² = sum(abs2, r-g.R)
    for (α,c) in zip(g.α,g.c)
        v += c*(2α/π)^(3/4)*exp(-α*rR²)
    end
    v
end

# ** STO-NG matrix elements

function LinearAlgebra.dot(a::STO_NG{T}, A, b::STO_NG{U}) where {T,U}
    s = zero(promote_type(T,U))
    RA = a.R
    RB = b.R
    R² = sum(abs2, RA-RB)

    for (α,c) in zip(a.α,a.c)
        for (β,d) in zip(b.α,b.c)
            s += c*d*primitive_dot(α, A, β, R², RA, RB)
        end
    end
    s
end

LinearAlgebra.dot(a::STO_NG, b::STO_NG) =
    LinearAlgebra.dot(a::STO_NG, I, b::STO_NG)

Base.:(*)(a::Adjoint{<:Any,<:STO_NG}, b::STO_NG) = dot(parent(a), b)
Base.:(*)(a::STO_NG, b::Adjoint{<:Any,<:STO_NG}) = conj(dot(parent(b), a))

function LinearAlgebra.Matrix(Ac::Adjoint{<:Any, <:AbstractVector{<:STO_NG{T}}},
                              O,
                              B::AbstractVector{<:STO_NG{U}}) where {T,U}
    m,n = length(Ac),length(B)
    M = zeros(promote_type(T,U), m, n)
    for i ∈ 1:m
        a = parent(Ac[i])
        for j ∈ 1:n
            b = B[j]
            M[i,j] = if O isa UniformScaling && i == j
                # This is a bit ugly, but the diagonal overlaps are
                # not numerically one.
                O.λ
            else
                dot(a, O, b)
            end
        end
    end
    M
end

common_center(α, RA, β, RB, αpβ⁻¹=inv(α+β)) =
    (α*RA + β*RB)*αpβ⁻¹

# *** Overlap integrals

function primitive_dot(α, I::UniformScaling, β, R², RA, RB)
    αpβ⁻¹ = inv(α+β)
    I.λ*(π*αpβ⁻¹)^(3/2)*exp(-α*β*R²*αpβ⁻¹)
end

# *** Kinetic operator

struct KineticOperator end
Base.show(io::IO, ::KineticOperator) = write(io, "-½∇²")

function primitive_dot(α, ::KineticOperator, β, R², RA, RB)
    αpβ⁻¹ = inv(α+β)
    μ = α*β*αpβ⁻¹
    μ*(3 - 2μ*R²)*(π*αpβ⁻¹)^(3/2)*exp(-μ*R²)
end

# *** Boys function

function boys0(t)
    if abs(t) < √(eps(t))
        one(t) - t/3
    else
        (√(π/t)*erf(√(t)))/2
    end
end

# *** Coulomb potential

struct CoulombPotential{T,C}
    Q::T
    R::C
end
Base.show(io::IO, VC::CoulombPotential) = write(io, "$(VC.Q)/|r-$(VC.R)|")

function primitive_dot(α, VC::CoulombPotential, β, RAB², RA, RB)
    αpβ⁻¹ = inv(α+β)
    RP = common_center(α, RA, β, RB, αpβ⁻¹)
    RPC² = sum(abs2, RP-VC.R)
    2π*αpβ⁻¹*VC.Q*exp(-α*β*αpβ⁻¹*RAB²)*boys0((α+β)*RPC²)
end

# *** Electron repulsion integrals

struct ElectronRepulsionIntegrals{Storage}
    storage::Storage
end

# Dumb, dense array storage

function Base.setindex!(V::ElectronRepulsionIntegrals{<:AbstractArray{<:Any,4}}, v, i, j, k, l)
    VV = V.storage
    VV[i,j,k,l] = VV[j,i,k,l] = VV[i,j,l,k] = VV[j,i,l,k] =
        VV[k,l,i,j] = VV[k,l,j,i] = VV[l,k,i,j] = VV[l,k,j,i] = v
end

Base.getindex(V::ElectronRepulsionIntegrals{<:AbstractArray{<:Any,4}}, i, j, k, l) =
    V.storage[i,j,k,l]

function unique_eris(fun::Function, n)
    for i ∈ 1:n
        ii = (i*(i-1))÷2
        for j ∈ i:n
            ij = ii + j
            for k ∈ 1:n
                kk = (k*(k-1))÷2
                for l ∈ k:n
                    kl = kk + l
                    kl > ij && continue
                    fun(i,j,k,l)
                end
            end
        end
    end
end

function num_unique_eris(n)
    N = 0
    unique_eris(n) do _,_,_,_
        N += 1
    end
    N
end

function Base.show(io::IO, ::MIME"text/plain", ElectronRepulsionIntegrals::ElectronRepulsionIntegrals)
    n = size(ElectronRepulsionIntegrals.storage,1)
    N = num_unique_eris(n)
    write(io, "$(n) ElectronRepulsionIntegrals stored in a $(typeof(ElectronRepulsionIntegrals.storage)) with $N unique integrals:\n")
    fmt = FormatExpr("#{1:<3d} {2:>3d} {3:>3d} {4:>3d} {5:>3d}: {6:<12.6e}")
    m = 1
    unique_eris(n) do i,j,k,l
        printfmtln(io, fmt, m, i, j, k, l, ElectronRepulsionIntegrals[i,j,k,l])
        m += 1
    end
end

function primitive_eri(α, β, γ, δ, RAB², RCD², RA, RB, RC, RD)
    # This is in charge cloud notation, i.e. RA and RB are in the same
    # coordinate r₁, with A being the ket and B the bra. Similarly for
    # RC and RD in the coordinate r₂.

    αpβ⁻¹ = inv(α+β)
    γpδ⁻¹ = inv(γ+δ)

    RP = common_center(α, RA, β, RB, αpβ⁻¹)
    RQ = common_center(γ, RC, δ, RD, γpδ⁻¹)
    RPQ² = sum(abs2, RP - RQ)

    2(π^(5/2))/((α+β)*(γ+δ)*√(α+β+γ+δ)) *
        exp(-α*β*αpβ⁻¹*RAB² - γ*δ*γpδ⁻¹*RCD²) *
        boys0((α+β)*(γ+δ)/(α+β+γ+δ)*RPQ²)
end

function ElectronRepulsionIntegrals(Φ::AbstractVector{<:STO_NG{T}}) where T
    # Find all electron repulsion between pairs A C of electrons from
    # the ket side and pairs B D of electrons from the bra
    # side. Formally, the pair A C needs to be conjugated, but since
    # we work with real Gaussians, and put any complexity into the
    # expansion coefficients, we can simply evaluate the integrals in
    # real arithmetic and employ symmetries such as V(2111) == V(1211)
    # == V(1121) == V(1112).

    n = length(Φ)
    V = ElectronRepulsionIntegrals(zeros(T, n, n, n, n))
    
    for i ∈ 1:n
        ii = (i*(i-1))÷2
        ϕi = Φ[i]
        RA = ϕi.R
        for j ∈ i:n
            ij = ii + j
            ϕj = Φ[j]
            RB = ϕj.R
            RAB² = sum(abs2, RA-RB)
            for k ∈ 1:n
                kk = (k*(k-1))÷2
                ϕk = Φ[k]
                RC = ϕk.R
                for l ∈ k:n
                    kl = kk + l
                    kl > ij && continue
                    ϕl = Φ[l]
                    RD = ϕl.R
                    RCD² = sum(abs2, RC-RD)
                    v = zero(T)
                    for (α,a) in zip(ϕi)
                        for (β,b) in zip(ϕj)
                            μ = a*b
                            for (γ,c) in zip(ϕk)
                                ν = μ*c
                                for (δ,d) in zip(ϕl)
                                    v += ν*d*primitive_eri(α, β, γ, δ, RAB², RCD², RA, RB, RC, RD)
                                end
                            end
                        end
                    end
                    V[i,j,k,l] = v
                end
            end
        end
    end
    
    V
end

# ** Standard contractions

"""
    sto_g(N, R[, ζ=1.0])

Standard contractions for ``N=1,2,3``, centred at ``R``, and with
exponent ``ζ``.
"""
function sto_g(N::Int, R, ζ=1.0)
    if N == 1
        STO_NG([1], [0], [0.270950], [1.0], R, ζ)
    elseif N == 2
        STO_NG([1,1], [0,0], [0.151623, 0.851819], [0.678914, 0.430129], R, ζ)
    elseif N == 3
        STO_NG([1,1,1], [0,0,0], [0.109818, 0.405771, 2.22766],
               [0.444635, 0.535328, 0.154329], R, ζ)
    else
        throw(ArgumentError("Invalid N = $N"))
    end
end
