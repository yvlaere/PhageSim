#=
Created on Friday 31 January 2020
Last update: Friday 21 Feb 2020

@author: Michiel Stock
michielfmstock@gmail.com

Simulation of bactium-phage interactions with three species of bacteria and
three species of phages.
=#

using DrWatson
quickactivate(@__DIR__, "PhageSim")
using PhageSim, Plots

pdecay = 0.1
pinfect = 0.1
tsteps = 500
N = 100
ninitbact = 100
ninitphages = 10
bustsize = 20

bactrules = BacteriaRules(0.4, 0.4, 0.1)
phagerules = PhageRules(pdecay, bustsize)
interactrules = InteractionRules([pinfect 0 0; 0 pinfect 0; 0 0 pinfect], bustsize)

# init the bacteria
bactgrid = BactGrid(nothing, N, N)
bactgrid[rand(CartesianIndices(bactgrid), ninitbact)] .= [Bacterium(rand(1:3)) for i in 1:ninitbact]

# init the phages
phagegrids = Array{Int,2}[]
for i in 1:3
    phagegrid = zeros(Int, N, N)
    phagegrid[rand(CartesianIndices(phagegrid), ninitphages)] .= 10
    push!(phagegrids, phagegrid)
end

anim = @animate for t=0:tsteps
    pbact = heatmap(species.(bactgrid), color=:darktest, colorbar=:none, framestyle=:none,
                title="bacteria (step $t)");
    pphagues1 = heatmap(phagegrids[1], color=:Reds, framestyle=:none,
                title="Phage sp. 1\n $(sum(phagegrids[1])) particles");
    pphagues2 = heatmap(phagegrids[2], color=:Greens, framestyle=:none,
                title="Phage sp. 2\n $(sum(phagegrids[2])) particles");
    pphagues3 = heatmap(phagegrids[3], color=:Blues, framestyle=:none,
                title="Phage sp. 3\n $(sum(phagegrids[3])) particles");
    plot(pbact, pphagues1, pphagues2, pphagues3)
    nbacts = nbacteria(bactgrid)
    #status = "step $t\n $(nbacts) bacteria, $(nphages) particles"
    #println(status)
    #title!(status)
    for phagegrid in phagegrids
        updatephages!(phagegrid, phagerules)
    end
    update!(bactgrid, phagegrids, bactrules, phagerules, interactrules)
end

gif(anim, plotsdir("bactphage3sp.gif"), fps = 10)
