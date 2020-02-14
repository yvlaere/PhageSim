#=
Created on Friday 27 December 2019
Last update: Friday 31 January 2020

@author: Michiel Stock
michielfmstock@gmail.com

Illustration of the phage movement.
=#

using DrWatson
quickactivate(@__DIR__, "PhageSim")

using PhageSim, Plots

pdecay = 0.01
nsteps = 500

grid = zeros(Int, 50, 50)
grid[rand(CartesianIndices(grid), 25)] .= 10

anim = @animate for t=1:nsteps+1
    heatmap(grid, title="step $t\n $(sum(grid)) particles",
                clims=(0, 10), color=:viridis)
    updatephages!(grid, R=2, pdecay=pdecay)
end

gif(anim, plotsdir("phagespread.gif"), fps = 10)
