#=
Created on Saturday 28 December 2019
Last update: Friday 31 January 2020

@author: Michiel Stock
michielfmstock@gmail.com

Implementation of the bacteria. For the moment, does not yet contain latent
phages.
=#

export AbstractBacterium, Bacterium, BactGrid
export isbacterium, species, copy, nbacteria, emptybactgrid
export updatebacteria!


abstract type AbstractBacterium end

struct Bacterium <: AbstractBacterium
    species::Int  # decribes the species of the bacterium
    #phage::Int  # either carries a latent phage (i) or not (0)
end

BactGrid = Array{Union{Nothing, Bacterium}}

#Bacterium(species::Int) = Bacterium(species)
isbacterium(state) = state isa AbstractBacterium

#haslatent(bact::Bacterium) = !isnothing(bact.phage)
#phage(bact::Bacterium) = bacterium.phage
species(bact::Bacterium) = bacterium.species

copy(bact::Bacterium) = Bacterium(bact.species)

"""
Count the number of bacteria in a grid.
"""
nbacteria(bactgrid) = count(isbacterium, bactgrid)



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
