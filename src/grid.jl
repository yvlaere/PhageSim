#=
Created on Friday 27 December 2019
Last update: Friday 21 Feb 2020

@author: Michiel Stock
michielfmstock@gmail.com

Function for working with grids
=#

export randnonident
export periodic, neighbors, regions

CI = CartesianIndex

"""
    periodic(I::CI, Ifirst::CI, Ilast::CartesianIndex)

Wraps I around a periodic boundary.
"""
function periodic(I::CI, Ifirst::CI, Ilast::CartesianIndex)
    CI([i < ifir ? i + ila : (i > ila ? i - ila : i) for (i, ifir, ila) in zip(Tuple.([I, Ifirst, Ilast])...)]...)
end

"""
    neighbors(A::AbstractArray)

Computes all the Moore neighbors of A, given periodic boundary conditions.
"""
function neighbors(A::AbstractArray)
    C = CartesianIndices(A)
    Ifirst, Ilast = first(C), last(C)
    I1, I0 = oneunit(Ifirst), zero(Ifirst)
    D = [d for d in reshape(-I1:I1, :) if d != I0]
    per = I -> periodic(I, Ifirst, Ilast)
    N = [per.([I] .+ D) for I in C]
    return N
end

"""
    region(A::AbstractArray)

Computes the Moore region for every position of A, given periodic boundary
conditions. E.g. this returns the Cartesian indices of the neighbors AND the
index of the cell itself.
"""
function regions(A::AbstractArray)
    C = CartesianIndices(A)
    Ifirst, Ilast = first(C), last(C)
    I1, I0 = oneunit(Ifirst), zero(Ifirst)
    D = [d for d in reshape(-I1:I1, :)]
    per = I -> periodic(I, Ifirst, Ilast)
    R = [per.([I] .+ D) for I in C]
    return R
end


"""
    randnonident(I, set)

Returns a random element `N` from `set` such that `N != I`.
Used to sample from a region.
"""
function randnonident(I, set)
    while true
        N = rand(set)
        N != I && return N
    end
end
