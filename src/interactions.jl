#=
Created on Friday 21 Feb 2020
Last update: Wednesday 26 March 2020

@author: Michiel Stock
michielfmstock@gmail.com

Model the interactions between bacteria and their phages
=#

export AbstractInteractionRules, InteractionRules
export update!

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
    function InteractionRules(Pinf::AbstractMatrix, burstsize::Number,
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
    plysogeny = interactionrules.lysogeny
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
    ρ = density(bactgrid, I, sp, R=R)
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

# compute the burstsize of an infected bacterium
# sampled by means of a Poisson distribution
burstsize(bact::AbstractBacterium, phagetype,
            interactionrules::AbstractInteractionRules) = rand(Poisson(interactionrules.burstsize))


function update!(grid::BactGrid,
                phagegrids,
                bactrules::AbstractBacteriaRules,
                phagerules::AbstractPhageRules,
                interactionrules::AbstractInteractionRules)
    # number of phages
    nphagetypes = length(phagegrids)
    # list of phage types, used for effecient shuffling
    phagetypes = collect(1:nphagetypes)
    bactcoors = Set(findall(isbacterium, grid))
    nbacts = length(bactcoors)
    C = CartesianIndices(grid)
    Ifirst, Ilast = first(C), last(C)
    IR = oneunit(Ifirst)
    for i in 1:nbacts
        I = rand(bactcoors)
        sI = bact = grid[I]  # state at position I
        delete!(bactcoors, I)
        # first check if the bacterium has a laten phage
        if haslatent(bact)
            # if it has a latent phage, it might lyse now
            if lyses(bact, grid, I, interactionrules::AbstractInteractionRules)
                # bacterium lyses, bye!
                phagetype = prophage(bact)
                phagegrids[phagetype][I] += burstsize(bact, phagetype, interactionrules)
                grid[I] = nothing
                sI = nothing
                continue
            end
        else
            # no latent phage, so check if there are any phages that might attack
            for phagetype in randperm!(phagetypes)
                nphages = phagegrids[phagetype][I]
                # if there are phages, will they infect the bacterium?
                if nphages > 0 && infects(bact, phagetype, nphages, interactionrules)
                    # yes, but lysogentic or lytic?
                    if lysogenic(bact, phagetype, interactionrules)
                        bact = prophage(bact, phagetype)  # add a prophage
                        grid[I] = bact
                        phagegrids[phagetype][I] -= 1
                    else
                        grid[I] = nothing  # kill bacterium
                        phagegrids[phagetype][I] += burstsize(bact, phagetype, interactionrules) - 1
                        break
                    end
                end
            end
        end
        isbacterium(grid[I]) || continue  # if killed, skip the rest
        # determine action
        neighborhood = max(Ifirst, I-IR):min(Ilast, I+IR)
        N = randnonident(I, neighborhood)
        sN = grid[N]
        sI, sN = updatebact(sI, sN, bactrules)
        grid[I], grid[N] = sI, sN
    end
end
