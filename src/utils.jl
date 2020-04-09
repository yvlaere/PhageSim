#=
Created on Wednesday 26 March 2020
Last update: Thursday 9 April 2020

@author: Michiel Stock
michielfmstock@gmail.com

General utilies for modelling.
=#

export insimplex, species_counts, species_counts_lat

"""
    insimplex(x::AbstractFloat...)

Test if values of `x` are in a simplex, i.e., non-negative and sum to one.
"""
insimplex(x::AbstractFloat...) = all(x .≥ 0.0) && sum(x) ≈ 1.0


"""
    species_counts(bactgrid::BactGrid, phagegrids; nspecies::Int)

Returns a tuple with a list counting all the bacteria of each species and a
list counting all phages of each species.
"""
function species_counts(bactgrid::BactGrid, phagegrids; nspecies::Int)
    # first count the bacteria of each species
    bactcount = [nbacteria(bactgrid, sp) for sp in 1:nspecies]
    phagecount = sum.(phagegrids)
    return bactcount, phagecount
end

"""
    species_counts_lat(bactgrid::BactGrid, phagegrids; nspecies::Int)

Returns a tuple with a list counting all the bacteria of each species, a
list counting all phages of each species and a matrix containing the number
latent phages per bacterium species.
"""
function species_counts_lat(bactgrid::BactGrid, phagegrids; nspecies::Int)
    # first count the bacteria of each species
    nphages = length(phagegrids)
    bactcount = [nbacteria(bactgrid, sp) for sp in 1:nspecies]
    phagecount = sum.(phagegrids)
    bacteria = filter(isbacterium, bactgrid)
    prophagescount = [count(
            b -> species(b) == i && haslatent(b) && prophage(b)==j, bacteria)
            for i in 1:nspecies, j in 1:nphages]
    return bactcount, phagecount, prophagescount
end
