#=
Created on Friday 21 Feb 2020
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

Model the interactions between bacteria and their phages
=#

export AbstractInteractionRules, InteractionRules
export update!

using Random, Distributions

abstract type AbstractInteractionRules end

struct InteractionRules{TA, TB} <: AbstractInteractionRules
    Pinf::TA  # matrix to determine the chance of infection host-phage
    burstsize::TB
end

"""
    infects(bact::AbstractBacterium, phagetype, nphages,
            interactionrules::AbstractInteractionRules)

Determine whether an infection by `phagetype` will take place.
"""
function infects(bact::AbstractBacterium, phagetype, nphages,
            interactionrules::AbstractInteractionRules)
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
        sI = bact = grid[I]
        delete!(bactcoors, I)
        # infections, randomly go over all phage types
        for phagetype in randperm!(phagetypes)
            nphages = phagegrids[phagetype][I]
            if nphages > 0 && infects(bact, phagetype, nphages, interactionrules)
                grid[I] = nothing  # kill bacterium
                phagegrids[phagetype][I] += burstsize(bact, phagetype, interactionrules) - 1
                break
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
