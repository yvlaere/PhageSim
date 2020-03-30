#=
Created on Friday 31 January 2020
Last update: Friday 21 Feb 2020

@author: Michiel Stock
michielfmstock@gmail.com

Simple simulation of bactium-phage interactions.
=#

using DrWatson
quickactivate(@__DIR__, "PhageSim")
using PhageSim, Plots

pdecay = 0.1
tsteps = 1000
N = 100
ninitbact = 25
ninitphages = 10

bactrules = BacteriaRules(0.4, 0.2, 0.1)
phagerules = PhageRules(0.2, 2)
interactrules = InteractionRules(reshape([0.8], 1, 1), 5)

# init the bacteria
bactgrid = BactGrid(nothing, N, N)
bactgrid[rand(CartesianIndices(bactgrid), ninitbact)] .= [Bacterium(1) for i in 1:ninitbact]

# init the phages
phagegrid = zeros(Int, N, N)
phagegrid[rand(CartesianIndices(phagegrid), ninitphages)] .= 10

anim = @animate for t=1:tsteps+1
    heatmap(isbacterium.(bactgrid), color=:Blues, colorbar=:none, framestyle=:none)
    locphage = Tuple.(findall(x-> x> 0, phagegrid))
    scatter!([x for (x, y) in locphage], [y for (x, y) in locphage],
        color=:orange, alpha=0.5, label="", markersize=0.2)
    nbacts = nbacteria(bactgrid)
    nphages = sum(phagegrid)
    status = "step $t\n $(nbacts) bacteria, $(nphages) particles"
    println(status)
    title!(status)
    step_phages!(phagegrid, phagerules)
    step_bacteria!(bactgrid, [phagegrid], bactrules, phagerules, interactrules)
end

gif(anim, plotsdir("bactphage.gif"), fps = 10)
