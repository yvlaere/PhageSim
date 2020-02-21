#=
Created on Friday 27 December 2019
Last update: Friday 21 Feb 2020

@author: Michiel Stock
michielfmstock@gmail.com

Illustration of the phage movement.
=#

using DrWatson
quickactivate(@__DIR__, "PhageSim")

using PhageSim, Plots

phagerules = PhageRules(0.01, 2)

grid = zeros(Int, 50, 50)
grid[rand(CartesianIndices(grid), 25)] .= 10

nsteps = 100

anim = @animate for t=1:nsteps+1
    heatmap(grid, title="step $t\n $(sum(grid)) particles",
                clims=(0, 10), color=:viridis, framestyle=:none)
    updatephages!(grid, phagerules)
end

gif(anim, plotsdir("phagespread.gif"), fps = 10)
