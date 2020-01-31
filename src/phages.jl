#=
Created on Friday 27 December 2019
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

Functions to model the spread of the phages.
=#

using Distributions

"""
updatephages!(grid::AbstractArray{T} where {T<:Integer};
                        R::Integer=1,
                        pdecay::Real=0.0,
                        poissonapprox=false)

Updating the spread of the phages by means of random diffusion and decay (`pdecay`).
By setting `poissonapprox` to true
a Poisson approximation is used instead of the exact, but slower Multinomial distribution.
"""
function updatephages!(grid::AbstractArray{T} where {T<:Integer};
                        R::Integer=1,
                        pdecay::Real=0.0,
                        poissonapprox=false,
                        nsteps::Union{Nothing,Int}=nothing
                        )
    @assert 0.0 ≤ pdecay ≤ 1.0 "`pdecay` should be in [0, 1], got $pdecay"
    @assert R > 0 && R < minimum(size(grid))
    pdecay == 1.0 && return (grid .= 0)
    ncells = length(grid)
    C = CartesianIndices(grid)
    Ifirst, Ilast = first(C), last(C)
    IR = R * oneunit(Ifirst)
    neigsize = (2R + 1)^2
    if nsteps isa Nothing
        nsteps = ncells
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
                # QUESTION: scale decay with neighsize
                nphages = rand(Binomial(nphages, 1.0 - pdecay / regsize))
            end
            grid[region] = rand(Multinomial(nphages, regsize))
        else
            avphages = (1.0 - pdecay / regsize) * nphages / regsize
            grid[region] = rand(Poisson(avphages), regsize)
        end
    end
    return grid
end

#=
updatephages!(grid::AbstractArray{T} where {T<:Integer};
                        R::Integer=1,
                        pdecay::Real=0.0,
                        poissonapprox=false) =
     updatephages!(grid, R=R, pdecay=pdecay, poissonapprox=poissonapprox,
                                        nsteps=length(grid))
=#
