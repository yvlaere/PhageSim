#=
Created on Saturday 28 December 2019
Last update: Monday 30 March 2020

@author: Michiel Stock
michielfmstock@gmail.com

Implementation of the bacteria. For the moment, does not yet contain latent
phages.
=#

export AbstractBacterium, Bacterium, BactGrid
export isbacterium, species, copy, nbacteria, emptybactgrid, prophage, haslatent
export density
export AbstractBacteriaRules, BacteriaRules, updatebacteria!


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
prophage(bact::Bacterium) = bact.phage
prophage(bact::Bacterium, i::Int) = Bacterium(bact.species, i)
species(bact::Bacterium) = bact.species
species(bactgrid::BactGrid) = bactgrid .|> species |> skipmissing |> unique
species(bact::Bacterium, i::Int) = bact.species == i
species(::Nothing) = missing
copy(bact::Bacterium) = Bacterium(bact.species, bact.phage)
nbacteria(bactgrid::BactGrid) = count(isbacterium, bactgrid)
nbacteria(bactgrid::BactGrid, sp::Int) = count(b->isbacterium(b) && species(b)==sp, bactgrid)

# computing the global and local density of the bacteria

"""
    density(bactgrid::BactGrid)

Computes the density of the bacteria in the grid, i.e., number of bacteria divided
by the total size of the grid.
"""
density(bactgrid::BactGrid) = nbacteria(bactgrid) / length(bactgrid)

"""
    density(bactgrid::BactGrid, sp::Int)

Computes the density of bacteria of species `sp` in the grid.
"""
density(bactgrid::BactGrid, sp::Int) = nbacteria(bactgrid, sp) / length(bactgrid)

"""
    density(bactgrid::BactGrid, I::CartesianIndex,
                        sp::Union{Nothing,Int}=nothing; R::Int=1)

Computes the local density of bacteria (or bacteria of species `sp` if provided)
in the region with radius `R` around position `I`.
"""
function density(bactgrid::BactGrid, I::CartesianIndex,
                    sp::Union{Nothing,Int}=nothing; R::Int=1)
    C = CartesianIndices(bactgrid)
    Ifirst, Ilast = first(C), last(C)
    IR = R * oneunit(Ifirst)
    region = max(I-IR, Ifirst):min(I+IR, Ilast)
    sp isa Nothing && return density(bactgrid[region])
    density(bactgrid[region], sp)
end

abstract type AbstractBacteriaRules end

struct BacteriaRules <: AbstractBacteriaRules
    prepr::Float64
    pmove::Float64
    pdie::Float64
    function BacteriaRules(prepr, pmove, pdie)
        @assert +(prepr, pmove, pdie) ≤ 1.0 &&
                all((prepr, pmove, pdie) .≥ 0) "behaviour of the bacteria should be valid probabilites"
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
