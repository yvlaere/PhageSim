#=
Created on Friday 27 December 2019
Last update: Friday 14 February 2019

@author: Michiel Stock
michielfmstock@gmail.com

Functions to model the spread of the phages.
=#

using Distributions

export AbstractPhageRules, PhageRules, phagedecay, updatephages!

abstract type AbstractPhageRules end

struct PhageRules <: AbstractPhageRules
    pdecay::Float64
    R::Int
    function PhageRules(pdecay, R=1)
        @assert 0.0 ≤ pdecay ≤ 1.0 "`pdecay` should be in [0, 1], got $pdecay"
        @assert R > 0 "R should be a positive integer"
        return new(pdecay, R)
    end
end

"""
    phagedecay(nphages, phagerules::PhageRules)

Samples the number of phages after applying decay
"""
function phagedecay(nphages, phagerules::PhageRules)
    nphages == 0 && return nphages
    psurv = 1.0 - phagerules.pdecay
    psurv == 0.0 && return 0
    psurv == 1.0 && return nphages
    return rand(Binomial(nphages, psurv))
end


"""
updatephages!(grid::AbstractArray{T} where {T<:Integer},
                        phagerules::AbstractPhageRules;
                        poissonapprox=false,
                        nsteps::Union{Nothing,Int}=nothing
                        )

Updating the spread of the phages by means of random diffusion and decay as specified by `phagerules`.
By setting `poissonapprox` to true, a Poisson approximation is used instead of the exact,
but slower Multinomial distribution.
"""
function updatephages!(grid::AbstractArray{T} where {T<:Integer},
                        phagerules::AbstractPhageRules;
                        poissonapprox=false,
                        nsteps::Union{Nothing,Int}=nothing
                        )
    ncells = length(grid)
    R, pdecay = phagerules.R, phagerules.pdecay
    C = CartesianIndices(grid)
    Ifirst, Ilast = first(C), last(C)
    IR = R * oneunit(Ifirst)
    neigsize = (2R + 1)^2
    if nsteps isa Nothing
        nsteps = ncells ÷ neigsize
    end
    for i in 1:nsteps
        # pick a cell
        I = rand(C)
        region = max(Ifirst, I - IR):min(Ilast, I + IR)
        regsize = length(region)
        nphages = sum(grid[region])
        nphages > 0 || continue  # only do the next if there are phages
        if !poissonapprox
            if pdecay > 0.0
                nphages = rand(Binomial(nphages, 1.0 - pdecay))
            end
            @inbounds grid[region] = rand(Multinomial(nphages, regsize))
        else
            avphages = (1.0 - pdecay) * nphages / regsize
            @inbounds grid[region] = rand(Poisson(avphages), regsize)
        end
    end
    return grid
end
