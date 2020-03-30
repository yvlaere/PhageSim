#=
Created on Wednesday 26 March 2020
Last update: Friday 30 March 2020

@author: Michiel Stock
michielfmstock@gmail.com

General utilies for modelling.
=#

export insimplex, species_counts

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
