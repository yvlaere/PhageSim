#=
Created on Saturday 28 December 2019
Last update: Friday 14 February 2020

@author: Michiel Stock
michielfmstock@gmail.com

Implementation of the bacteria. For the moment, does not yet contain latent
phages.
=#

export AbstractBacterium, Bacterium, BactGrid
export isbacterium, species, copy, nbacteria, emptybactgrid, phage, haslatent
export BacteriaRules, updatebacteria!


abstract type AbstractBacterium end

struct Bacterium <: AbstractBacterium
    species::Int  # decribes the species of the bacterium
    phage::Int  # either carries a latent phage (i) or not (0)
end

Bacterium(species::Int) = Bacterium(species, 0)

"""Structure of bacterial grid"""
BactGrid = Array{Union{Nothing, Bacterium}}
emptybactgrid(dims...) = Array{Union{Nothing, Bacterium}}(nothing, dims...)

isbacterium(state) = state isa AbstractBacterium
haslatent(bact::Bacterium) = bact.phage != 0
phage(bact::Bacterium) = bact.phage
species(bact::Bacterium) = bact.species
species(bact::Bacterium, i::Int) = bact.species == i
species(::Nothing) = missing
copy(bact::Bacterium) = Bacterium(bact.species, bact.phage)
nbacteria(bactgrid) = count(isbacterium, bactgrid)

abstract type AbstractBacteriaRules end

struct BacteriaRules <: AbstractBacteriaRules
    prepr::Float64
    pmove::Float64
    pdie::Float64
    function BacteriaRules(prepr, pmove, pdie)
        @assert prepr ≥ 0 &&
                pmove ≥ 0 &&
                pdie ≥ 0 &&
                +(prepr, pmove, pdie) ≤ 1 "behaviour of the bacteria should be valid probabilites"
        new(prepr, pmove, pdie)
    end
end


updatebact(s1::AbstractBacterium, s2::AbstractBacterium, bacteriarules::BacteriaRules) = s1, s2

"""
    updatebact(s1::AbstractBacterium, s2::Nothing, bacteriarules::BacteriaRules)

Rules for updating the bacteria without the phage component.
"""
function updatebact(s1::AbstractBacterium, s2::Nothing, bacteriarules::BacteriaRules)
    r = rand()
    if r ≤ bacteriarules.prepr
        # reproduce
        return s1, copy(s1)
    elseif r ≤ bacteriarules.prepr + bacteriarules.pmove
        # move
        return nothing, s1
        # die
    elseif r ≤ bacteriarules.prepr + bacteriarules.pmove + bacteriarules.pdie
        return nothing, nothing
    else
        return s1, s2
    end
end


#=
function updatebacteria!(bactgrid, phagegrid; R=1, pinfect=0.9, burstsize=10)
    pmove = 0.6
    preproduce = 0.2
    bactcoors = Set(findall(isbacterium, bactgrid))
    nbacts = length(bactcoors)
    C = CartesianIndices(bactgrid)
    Ifirst, Ilast = first(C), last(C)
    IR = R * Ifirst
    for i in 1:nbacts
        I = rand(bactcoors)
        delete!(bactcoors, I)
        # check infections
        nphages = phagegrid[I]
        if nphages > 0 && rand() > (1.0-pinfect)^nphages
            bactgrid[I] = nothing
            phagegrid[I] += burstsize - 1  #IDEA: change by Poison?
            continue
        end
        region = max(Ifirst, I-IR):min(Ilast, I+IR)
        # determine action
        r = rand()
        if r < pmove
            Inew = rand(region)
            if bactgrid[Inew] == nothing
                bactgrid[I], bactgrid[Inew] = bactgrid[Inew], bactgrid[I]
            end
        else r < pmove + preproduce
            Inew = rand(region)
            if bactgrid[Inew] == nothing
                bactgrid[Inew] = copy(bactgrid[I])
            end
        end
    end
    return bactgrid
end
=#
