#=
Created on Friday 27 December 2019
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

Functions to model the spread of the phages.
=#

using Distributions

"""
updatephages!(grid::AbstractArray{T} where {T<:Integer},
                        R;
                        pdecay::Real=0.0,
                        poissonapprox=false)

Updating the spread of the phages by means of random diffusion and decay (`pdecay`).
Needs the array `R` containing the regions. By setting `poissonapprox` to true
a Poisson approximation is used instead of the exact, but slower Multinomial distribution.
"""
function updatephages!(grid::AbstractArray{T} where {T<:Integer},
                        R;
                        pdecay::Real=0.0,
                        poissonapprox=false)
    @assert 0.0 ≤ pdecay ≤ 1.0 "`pdecay` should be in [0, 1], got $pdecay"
    pdecay == 1.0 && return (grid .= 0)
    ncells = length(grid)
    C = CartesianIndices(grid)
    neigsize = length(first(R))
    # every cell will appear several times (on average)
    pdecay /= neigsize
    for i in 1:ncells
        # pick a cell
        I = rand(C)
        nphages = sum(grid[R[I]])
        nphages > 0 || continue  # only do the next if there are phages
        if !poissonapprox
            if pdecay > 0.0
                nphages = rand(Binomial(nphages, 1.0-pdecay))
            end
            grid[R[I]] .= rand(Multinomial(nphages, neigsize))
        else
            avphages = (1.0 - pdecay) * nphages / neigsize
            grid[R[I]] .= rand(Poisson(avphages), neigsize)
        end
    end
    return grid
end
