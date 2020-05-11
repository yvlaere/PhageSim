#alle figuren voor onderzoeksvraag1
#worden in 1 keer aangemaakt

using DrWatson
quickactivate(@__DIR__, "PhageSim")
using PhageSim
using Plots
using Statistics
using DataFrames
import BSON: @save, @load

#scenario
#ninitbact, initphages, pdie, prepr
ldens = [50, 500, 0.15, 0.30]
mdens = [200, 2000, 0.10, 0.35]
hdens = [500, 5000, 0.05, 0.40]
dens = [ldens, mdens, hdens]

#situatie
#plysogeny, plysis
#plysogeny = bij contact kans op lysogenie, rest voor lyse
#plyse = bij lysogenie, kans op lyse
sit1 = [0.0, 1.0]
sit2 = [0.3, 0.7]
sit3 = [0.5, 0.5]
sit4 = [0.6, 0.4]
sit5 = [0.7, 0.3]
sit6 = [0.8, 0.2]
sit7 = [0.85, 0.1]
sit = [sit1, sit2, sit3, sit4, sit5, sit6, sit7]

#aantal simulaties van 1 situatie
nsims = 10
overzicht = zeros(Float32, length(dens)*length(sit), 2)
ratio = zeros(Float32, length(sit))
plotname = "degelijke test2"

for i = 1:length(dens)
    for j = 1:length(sit)

        nbact = zeros(Int64, nsims, 3)
        nfaag = zeros(Int64, nsims, 3)

        name = "scenario" * string(i) * "situatie" * string(j)  # give a name to your simulation
        println(name)
        # PARAMETERS OF THE SIMULATION
        # ----------------------------

        # bacteria
        # --------

        prepr = dens[i][4]
        pmove = 0.4
        pdie = dens[i][3]

        # phages
        # -----

        pdecay = 0.02  # probability that a viron decays
        Rphages = 5  # mixing radius for the dispersion of the phages

        # interactions
        # ------------

        # matrix to determine the chance of infection host-phage, nbact x nphages
        Pinf =  [0.2 0.01 0.01; 0.01 0.2 0.01; 0.01 0.01 0.2]
        burstsize = 10  # average number of virons
        plysogeny = sit[j][1]  # probability of entering the lysogentic cycle, either probability
        # for all phages or a list of probabilities, one for each phage
        plysis = sit[j][2]  # probability that an infected bacterium lyses in every step
        Rqs = 0  # radius for quorum sensing, lysis probability is proportional to the local density
        # 0 means that it is independent

        # simulation stuff
        # ----------------

        nsteps = 200
        D = 100  # size of the grid
        ninitbact = floor(Int64, dens[i][1])
        ninitphages = floor(Int64, dens[i][2])

        # PROCESSING THE PARAMETERS
        # -------------------------

        nbactsp, nphagesp = size(Pinf)

        simpars = @dict prepr pdie  plysogeny plysis ninitbact

        bactrules = BacteriaRules(prepr, pmove, pdie)
        phagerules = PhageRules(pdecay, Rphages)
        interactrules = InteractionRules(Pinf, burstsize, plysogeny, plysis, Rqs)

        # SETTING UP THE INITIAL CONDITIONS
        # ---------------------------------

        # init the bacteria
        bactgrid = BactGrid(nothing, D, D)  # empty grid
        bactgrid[rand(CartesianIndices(bactgrid), ninitbact)] .= [Bacterium(rand(1:nbactsp)) for i in 1:ninitbact]

        # init the phages
        phagegrids = Array{Int,2}[]
        for i in 1:nphagesp
            phagegrid = zeros(Int, D, D)
            phagegrid[rand(CartesianIndices(phagegrid), ninitphages)] .= 10  # place 10 phages
            push!(phagegrids, phagegrid)
        end

        # RUNNING THE SIMULATION
        # ----------------------

        # perform the simulations, the function `species_counts`, with keyword arguments
        # `nspecies` for the number of bacterium species, count the number of bacteria
        # and the number of phage sp in every time step
        for k = 1:nsims
            results = simulate!(bactgrid, phagegrids, bactrules, phagerules, interactrules, nsteps; resultfun=species_counts, nspecies=nbactsp)
            nbact[k, 1:3] = last(results)[1]
            nfaag[k, 1:3] = last(results)[2]
        end

        overzicht[(i - 1)*length(sit) + j, 1] = mean(mean(nbact, dims = 1))
        overzicht[(i - 1)*length(sit) + j, 2] = mean(mean(nfaag, dims = 1))
        # DATA STORAGE AND VISUALIZATION
        # ------------------------------

        # generate a unique name
        #fname = savename(name, simpars)

        # save the simulation in a BSON file
        #@save datadir(fname) * ".bson" results

        # can be loaded using
        # @load datadir(fname) * ".bson"

        # plotting

        #bactres = [res[1][sp] for res in results, sp in 1:nbactsp]
        #phageres = [res[2][sp] for res in results, sp in 1:nphagesp]

        #p = plot(plot(bactres, labels=["sp. 1" "sp. 2" "sp. 3"], xlabel="step",
        #        ylabel="number of bacteria", title="Bacteria composition"),
        #    plot(phageres, labels=["sp. 1" "sp. 2" "sp. 3"], ls=:dash, xlabel="step",
        #            ylabel="number of phages", title="Phage composition"))

        #savefig(p, plotsdir(fname) * ".png")
    end
end

print(convert(DataFrame, overzicht))

#plysogeny/plyse
for k = 1:length(sit)
    ratio[k] = sit[k][1]/sit[k][2]
end

bactplotdata = [overzicht[1:length(sit), 1]/mean(overzicht[1:length(sit), 1]),
 overzicht[length(sit) + 1:2*length(sit), 1]/mean(overzicht[length(sit) + 1:2*length(sit), 1]),
 overzicht[2*length(sit) + 1:3*length(sit), 1]/mean(overzicht[2*length(sit) + 1:3*length(sit), 1])]
faagplotdata = [overzicht[1:length(sit), 2]/mean(overzicht[1:length(sit), 2]),
 overzicht[length(sit) + 1:2*length(sit), 2]/mean(overzicht[length(sit) + 1:2*length(sit), 2]),
 overzicht[2*length(sit) + 1:3*length(sit), 2]/mean(overzicht[2*length(sit) + 1:3*length(sit), 2])]
p = plot(plot(ratio, bactplotdata, labels=["milieu1" "milieu2" "milieu3"], xlabel = "plysogeny/plyse", ylabel="nr. bact/gem(nr. bact)", title="Bacteria composition"),
 plot(ratio, faagplotdata, labels=["milieu1" "milieu2" "milieu3"], xlabel = "plysogeny/plyse", ylabel="nr. faag/gem(nr. faag)", title="Faag composition"))

simpars = @dict nsims nsteps
fname = savename(plotname, simpars)
savefig(p, plotsdir(fname) * ".png")
