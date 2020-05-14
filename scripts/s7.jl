
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

plotname = "s7plotsfinal2"
nsteps = 200
overzicht = zeros(Float64, nsteps + 1, 2*length(dens))

for i = 1:length(dens)
        j = 7
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

        results = simulate!(bactgrid, phagegrids, bactrules, phagerules, interactrules, nsteps; resultfun=species_counts, nspecies=nbactsp)
        bactres = [res[1][sp] for res in results, sp in 1:nbactsp]./ninitbact
        phageres = [res[2][sp] for res in results, sp in 1:nphagesp]./ninitphages
        print(bactres)
        print(mean(bactres, dims = 2))
        overzicht[1:nsteps + 1, i] = mean(bactres, dims = 2)
        overzicht[1:nsteps + 1, i + 3] = mean(phageres, dims = 2)

        print(convert(DataFrame, overzicht))
        # DATA STORAGE AND VISUALIZATION
        # ------------------------------

        # generate a unique name
        #fname = savename(name, simpars)

        # save the simulation in a BSON file
        #@save datadir(fname) * ".bson" results

        # can be loaded using
        # @load datadir(fname) * ".bson"

        # plotting

        #bactres = [res[1][sp] for res in results, sp in 1:nbactsp]./ninitbact
        #phageres = [res[2][sp] for res in results, sp in 1:nphagesp]./ninitphages

        #p = plot(plot(bactres, labels=["sp. 1" "sp. 2" "sp. 3"], xlabel="step",
        #        ylabel="number of bacteria", title="Bacteria composition"),
        #    plot(phageres, labels=["sp. 1" "sp. 2" "sp. 3"], ls=:dash, xlabel="step",
        #            ylabel="number of phages", title="Phage composition"))

        #savefig(p, plotsdir(fname) * ".png")
end

p = plot(plot(overzicht[:, 1:3], labels=["milieu1" "milieu2" "milieu3"], xlabel="step",
        ylabel="number of bacteria", title="Bacteria composition"),
    plot(overzicht[:, 4:6], labels=["milieu1" "milieu2" "milieu3"], ls=:dash, xlabel="step",
            ylabel="number of phages", title="Phage composition"))

#p = plot(plot(overzicht[1], ),
# plot(overzicht[2]))

simpars = @dict nsteps
fname = savename(plotname, simpars)
savefig(p, plotsdir(fname) * ".png")
