#=
Created on Friday 31 January 2020
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

Simple simulation of bactium-phage interactions.
=#

using DrWatson
quickactivate(@__DIR__, "PhageSim")
using PhageSim, Plots

pdecay = 0.1
tsteps = 500
N = 50
ninitbact = 25
ninitphages = 10

# init the bacteria
bactgrid = BactGrid(nothing, N, N)
bactgrid[rand(CartesianIndices(bactgrid), ninitbact)] .= [Bacterium(1) for i in 1:ninitbact]

# init the phages
phagegrid = zeros(Int, N, N)
phagegrid[rand(CartesianIndices(phagegrid), ninitphages)] .= 10

anim = @animate for t=1:tsteps+1
    heatmap(1isbacterium.(bactgrid), color=:Blues, colorbar=:none, framestyle=:none)
    locphage = Tuple.(findall(x-> x> 0, phagegrid))
    scatter!([x for (x, y) in locphage], [y for (x, y) in locphage],
        color=:orange, alpha=0.5, label="", markersize=0.2)
    nbacts = nbacteria(bactgrid)
    nphages = sum(phagegrid)
    status = "step $t\n $(nbacts) bacteria, $(nphages) particles"
    println(status)
    title!(status)
    updatephages!(phagegrid, R=1, pdecay=pdecay)
    updatebacteria!(bactgrid, phagegrid, burstsize=3, pinfect=0.8)
end

gif(anim, plotsdir("bactphage.gif"), fps = 10)
