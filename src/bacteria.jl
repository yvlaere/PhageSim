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
export AbstractBacteriaRules, BacteriaRules, HeteroBacteriaRules, ProphageBacteriaRules
export bacteriaprobs, updatebacteria!


abstract type AbstractBacterium end

struct Bacterium <: AbstractBacterium
    species::Int  # decribes the species of the bacterium
    phage::Int  # either carries a latent phage (i) or not (0)
end

Bacterium(species::Int) = Bacterium(species, 0)

"""Structure of bacterial grid"""
BactGrid = Array{Union{Nothing, Bacterium}}
emptybactgrid(dims...) = Array{Union{Nothing, Bacterium}}(nothing, dims...)

"""Check if a state is a `AbstractBacterium`"""
isbacterium(state) = state isa AbstractBacterium

"""
    haslatent(bact::Bacterium)

Test whether a bacterium has a laten prophage
"""
haslatent(bact::Bacterium) = bact.phage != 0

"""
    prophage(bact::Bacterium)

Return the prophage of a bacterium. Return `missing` if bact has no prophage.
"""
prophage(bact::Bacterium) = bact.phage == 0 ? missing : bact.phage

"""
    prophage(bact::Bacterium, i::Int)

Returns a NEW bacterium with a prophage of type `i`.
"""
prophage(bact::Bacterium, i::Int) = Bacterium(bact.species, i)

"""
    species(bact::Bacterium)

Returns the species of a bacterium.
"""
species(bact::Bacterium) = bact.species

"""
    species(bactgrid::BactGrid)

Returns a list of the species of bacteria in a grid.
"""
species(bactgrid::BactGrid) = bactgrid .|> species |> skipmissing |> unique |> sort!

"""
    species(bact::Bacterium, sp::Int)

Test if `bact` if of species `sp`.
"""
species(bact::Bacterium, sp::Int) = bact.species == sp

species(::Nothing) = missing

"""
    copy(bact::Bacterium)

Make a copy of a bacterium.
"""
copy(bact::Bacterium) = Bacterium(bact.species, bact.phage)

"""
    nbacteria(bactgrid::BactGrid)

Counts the number of bacteria in `bactgrid`.
"""
nbacteria(bactgrid::BactGrid) = count(isbacterium, bactgrid)

"""
    nbacteria(bactgrid::BactGrid, sp::Int)

Counts the number of bacteria of species `sp` in `bactgrid`.
"""
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

"""
Rules for when all bacteria behave the same.
"""
struct BacteriaRules <: AbstractBacteriaRules
    prepr::Float64
    pmove::Float64
    pdie::Float64
    function BacteriaRules(prepr::Float64, pmove::Float64, pdie::Float64)
        @assert +(prepr, pmove, pdie) ≤ 1.0 &&
                all((prepr, pmove, pdie) .≥ 0) "behaviour of the bacteria should be valid probabilites"
        new(prepr, pmove, pdie)
    end
end

"""
    bacteriaprobs(br::BacteriaRules, bact::AbstractBacterium)

Get the behaviour parameters for a specific bacterium.
"""
bacteriaprobs(br::BacteriaRules, bact::AbstractBacterium) = (br.prepr, br.pmove, br.pdie)

"""
Rules for when species of bacteria might show different behaviours.
"""
struct HeteroBacteriaRules <: AbstractBacteriaRules
    prepr::Array{Float64,1}
    pmove::Array{Float64,1}
    pdie::Array{Float64,1}
end

function BacteriaRules(prepr::Array{Float64,1}, pmove::Array{Float64,1}, pdie::Array{Float64,1})
    @assert all(.+(prepr, pmove, pdie) .≤ 1.0) && all(prepr .≥ 0) && all(pmove .≥ 0) &&
             all(pdie .≥ 0) "behaviour of the bacteria should be valid probabilites"
    HeteroBacteriaRules(prepr, pmove, pdie)
end

"""
Rules for when a bacterium with a prophage exhibits a different behaviour from
a bacterium without a prophage.
"""
struct ProphageBacteriaRules <: AbstractBacteriaRules
    probsnoinf::Tuple{Float64,Float64,Float64}
    probsinf::Tuple{Float64,Float64,Float64}
    precover::Float64
end

"""
    bacteriaprobs(br::BacteriaRules, bact::AbstractBacterium)

Get the behaviour parameters for a specific bacterium.
"""
function bacteriaprobs(br::HeteroBacteriaRules, bact::AbstractBacterium)
    sp = species(bact)
    return (br.prepr[sp], br.pmove[sp], br.pdie[sp])
end

function BacteriaRules(probsnoinf::Tuple{Float64,Float64,Float64},
                probsinf::Tuple{Float64,Float64,Float64},
                precover::Float64=0.0)
    @assert sum(probsnoinf) ≤ 1.0 &&
            all(probsnoinf .≥ 0) "behaviour of the bacteria should be valid probabilites"
    @assert sum(probsinf) ≤ 1.0 &&
            all(precover .≥ 0) "behaviour of the bacteria should be valid probabilites"
    @assert 0 ≤ precover ≤ 1 "`precover` should be a valid probability"
    ProphageBacteriaRules(probsnoinf, probsinf, precover)
end

"""
    bacteriaprobs(br::ProphageBacteriaRules, bact::AbstractBacterium)

Get the behaviour parameters depending on whether the bacterium has a prophage
or not.
"""
function bacteriaprobs(br::ProphageBacteriaRules, bact::AbstractBacterium)
    if haslatent(bact)
        return br.probsinf
    else
        return br.probsnoinf
    end
end

function updatebact(s1::AbstractBacterium, s2::AbstractBacterium, bacteriarules::AbstractBacteriaRules)
    prepr, pmove, pdie = bacteriaprobs(bacteriarules, s1)
    if rand() < pdie
        return nothing, s2
    else
        return s1, s2
    end
end


"""
    updatebact(s1::AbstractBacterium, s2::Nothing, bacteriarules::AbstractBacteriaRules)

Rules for updating the bacteria without the phage component.
"""
function updatebact(s1::AbstractBacterium, s2::Nothing, bacteriarules::AbstractBacteriaRules)
    # get the rules
    prepr, pmove, pdie = bacteriaprobs(bacteriarules, s1)
    r = rand()
    if r ≤ prepr
        # reproduce
        return s1, copy(s1)
    elseif r ≤ prepr + pmove
        # move
        return nothing, s1
        # die
    elseif r ≤ prepr + pmove + pdie
        return nothing, nothing
    else
        return s1, s2
    end
end
