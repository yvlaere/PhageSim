#=
Created on Friday 27 December 2019
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

Illustration of the phage movement.
=#

using PhageSim, Plots

pdecay = 0.01
nsteps = 100

grid = zeros(Int, 50, 50)
grid[rand(CartesianIndices(grid), 25)] .= 1000
R = regions(grid)

anim = @animate for t=1:nsteps+1
    heatmap(grid, title="step $t\n $(sum(grid)) particles")
    updatephages!(grid, R, pdecay=pdecay)
end

gif(anim, "phagespread.gif", fps = 10)
