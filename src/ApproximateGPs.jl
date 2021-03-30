module ApproximateGPs

using AbstractGPs
using GPLikelihoods

abstract type AbstractVariationalGP <: AbstractGPs.AbstractGP end

const JITT = 1e-8
struct VariationalGP{Tgp<:LatentFiniteGP,Tμ<:AbstractVector,TΣ} <: AbstractVariationalGP
    f::Tgp
    μ::Tμ
    Σ::TΣ
end

function VariationalGP(gp::LatentFiniteGP; μ::AbstractVector=randn(length(gp)), Σ::AbstractMatrix=Matrix{Float64}(I(length(gp))), kwargs...)
    length(μ) == length(gp) || error("Initial variational mean has wrong dimensions $(length(μ))")
    size(Σ, 1) == size(Σ, 2) == length(gp) || error("Initial variational covariance has wrong dimensions $(size(Σ))")
    return VariationalGP(gp, μ, Σ)
end

function VariationalGP(gp::FiniteGP, lik; kwargs...)
    return VariationalGP(LatentFiniteGP(gp, lik); kwargs...)
end

function VariationalGP(gp::GP, X::AbstractVector, lik; kwargs...)
    jitt = haskey(kwargs, :jitt) ? kwargs[:jitt] : JITT
    return VariationalGP(LatentFiniteGP(gp(X, jitt)), lik; kwargs...)
end

function fit!(gp::VariationalGP, y::AbstractVector, optimiser::GPOptimiser; n_iter::Int=100, callback=nothing, kwargs...)
    return optimiser(gp, y; n_iter, kwargs...)
end 


abstract type GPOptimiser end

struct AnalyticVI <: GPOptimiser

end

(opt::AnalyticVI)(gp::VariationalGP, y; n_iter, kwargs...)
    aug_var = create_augmented_variables(gp.f)
    for i in 1:n_iter
        update_aug_var!(aug_var, gp.f, gp.μ, gp.Σ)
        update_var!(gp, y, aug_var)
    end
    return gp
end

function create_augmented_variables(f::LatentFiniteGP)
    return create_augmented_variables(f.lik, length(f))
end

function create_augmented_variables(l::BernoulliLikelihood{<:LogitLink}, n::Int)
    return (θ=rand(n), c=rand(n))
end

function update_aug_var!(aug_var, f, μ, Σ)
    update_aug_var!(aug_var, f.lik, μ, diag(Σ))
end

function update_aug_var!(aug_var, l::BernoulliLikelihood, μ, diagΣ)
    @. aug_var.c = sqrt(diagΣ + abs2(μ))
    @. aug_var.θ = 0.5 * tanh(0.5 * l.c) / l.c
end

function update_var!(gp, y, aug_var)
    return update_var!(gp.f.lik, y, gp.μ, gp.Σ, cov(gp.f), mean(gp.f), aug_var)
end

function update_var!(l, y, μ, Σ, K, μ₀, aug_var)
    Σ .= inv(Symmetric(inv(K) + 2 * Diagonal(∇E_Σ(l, aug_var, y))))
    μ .= Σ * (∇E_μ(l, y, aug_var) + K \ μ₀)
    return nothing
end

function ∇E_μ(::BernoulliLikelihood{<:LogitLink}, y, aug_var)
    0.5 * y
end

function ∇E_Σ(::BernoulliLikelihood{<:LogitLink}, y, aug_var)
    0.5 * aug_var.θ
end

function elbo(gp::VariationalGP)
    
end

struct SparseVariationalGP{Tgp<:LatentFiniteGP,Tμ<:AbstractVector,TΣ} <: AbstractVariationalGP
    gp::LatentGP
    μ::Tμ
    Σ::TΣ
end

end
