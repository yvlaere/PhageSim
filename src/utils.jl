#=
Created on Wednesday 26 March 2020
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

General utilies for modelling.
=#

export insimplex

"""
    insimplex(x::AbstractFloat...)

Test if values of `x` are in a simplex, i.e. non-negative and normalized.
"""
insimplex(x::AbstractFloat...) = all(x .≥ 0.0) && sum(x) ≈ 1.0
