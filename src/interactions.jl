#=
Created on Friday 21 Feb 2020
Last update: Wednesday 15 April 2020

@author: Michiel Stock
michielfmstock@gmail.com

Model the interactions between bacteria and their phages
=#

export AbstractInteractionRules, InteractionRules

using Random, Distributions

abstract type AbstractInteractionRules end

"""
Rules for when only a lytic life cycle is present.
"""
struct InteractionRules{TI,TB,TLG} <: AbstractInteractionRules
    Pinf::TI  # matrix to determine the chance of infection host-phage
    burstsize::TB  # average number of virons
    plysogeny::TLG  # probability of entering the lysogentic cycle
    plysis::Float64  # probability that an infected bacterium lyses
    R::Int  # radius for quorum sensing
    function InteractionRules(Pinf, burstsize::Number,
            plysogeny::Union{BP,AbstractVector{BP}}=false,
            plysis=1.0, R=0)  where {BP <: Union{Bool, AbstractFloat}}
        @assert all(0 .≤ plysogeny .≤ 1) "All values for `plysogeny` should be Booleans or probabilities"
        @assert !(plysogeny isa AbstractVector) || size(Pinf, 2) == length(plysogeny) "`plysogeny` has incorrect size"
        @assert 0 ≤ plysis ≤ 1 "`plysis` should be a valid probability"
        @assert R ≥ 0 "radius for quorum sensing `R` should be nonzero"
        new{typeof(Pinf),typeof(burstsize),typeof(plysogeny)}(Pinf, burstsize, plysogeny,plysis)
    end
end

"""
    lysogenic(bact::AbstractBacterium, phagetype::Int,
                    interactionrules::AbstractInteractionRules)

Determines whether a phage will enter a lysogentic phase with the host as a
prophage or whether it will enter a lytic phase and kill its host immediately.
"""
function lysogenic(bact::AbstractBacterium, phagetype::Int,
                    interactionrules::AbstractInteractionRules)
    plysogeny = interactionrules.plysogeny
    plysogeny isa Bool && return plysogeny
    plysogeny isa Number && return rand() < plysogeny
    plysogeny isa AbstractVector{Bool} && return plysogeny[phagetype]
    plysogeny isa AbstractVector && return rand() < plysogeny[phagetype]
    error("behaviour of `plysogeny` not defined (type is $(typeof(plysogeny)))")
end

"""
    lyses(bact::AbstractBacterium, grid::BactGrid, I::CartesianIndex,
                                interactionrules::AbstractInteractionRules)

Determines whether a bacterium with a latent phage will lyse. Is dependent of
the local density of the bacterial species and the `interactionrules`.
"""
function lyses(bact::AbstractBacterium, grid::BactGrid, I::CartesianIndex,
                            interactionrules::AbstractInteractionRules)
    # compute the local density
    R = interactionrules.R
    sp = species(bact)
    ρ = density(grid, I, sp, R=R)
    # so the probability of lysis is proportional to the density
    # THE BACTERIUM ITSELF INCLUDED!
    return rand() < ρ * interactionrules.plysis
end


"""
    infects(bact::AbstractBacterium, phagetype::Int, nphages,
            interactionrules::AbstractInteractionRules)

Determine whether an infection by `phagetype` will take place.
"""
function infects(bact::AbstractBacterium, phagetype::Int, nphages,
            interactionrules::AbstractInteractionRules)
    haslatent(bact) && return false  # if bactium has a latent phage it is protected
    p = interactionrules.Pinf[species(bact), phagetype]
    p == 0.0 && return false
    # probability that at least one phage succesfully infect the bacterium
    return (1.0 - p)^nphages < rand()
end

"""
Compute the burstsize of an infected bacterium, sampled from a Poisson
distribution.
"""
burstsize(bact::AbstractBacterium, phagetype,
            interactionrules::AbstractInteractionRules) = rand(Poisson(interactionrules.burstsize))
