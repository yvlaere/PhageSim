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
ldens = [100, 100, 0.15, 0.30]
mdens = [100, 100, 0.10, 0.35]
hdens = [100, 100, 0.05, 0.40]
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
nsims = 1
overzicht = zeros(Float32, length(dens)*length(sit), 2)
ratio = zeros(Float32, length(sit))
plotname = "plotloglog7"
nsteps = 10

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
            println("sim done")
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

#lysogeny/(1 - plysogeny
for k = 1:length(sit)
    ratio[k] = round(sit[k][1]/(1 - sit[k][1]), digits = 5)#sit[k][2]
end

#logoverzicht = log.(10, overzicht)
#for l = 1:2
#    for m = 1:length(dens)*length(sit)
#        if logoverzicht[m, l] == -Inf
#            logoverzicht[m, l] = 0
#        end
#    end
#end
#print(convert(DataFrame, logoverzicht))

#logbactplotdata = [logoverzicht[1:length(sit), 1]/mean(logoverzicht[1:length(sit), 1]),
# logoverzicht[length(sit) + 1:2*length(sit), 1]/mean(logoverzicht[length(sit) + 1:2*length(sit), 1]),
# logoverzicht[2*length(sit) + 1:3*length(sit), 1]/mean(logoverzicht[2*length(sit) + 1:3*length(sit), 1])]
#logfaagplotdata = [logoverzicht[1:length(sit), 2]/mean(logoverzicht[1:length(sit), 2]),
# logoverzicht[length(sit) + 1:2*length(sit), 2]/mean(logoverzicht[length(sit) + 1:2*length(sit), 2]),
# logoverzicht[2*length(sit) + 1:3*length(sit), 2]/mean(logoverzicht[2*length(sit) + 1:3*length(sit), 2])]
#p1 = plot(plot(ratio, logbactplotdata, labels=["arm milieu" "matig milieu" "rijk milieu"], xlabel = "lysogeny/lyse", ylabel="nr. bact/gem(nr. bact)", title="Voorkomen bacterien"),
# plot(ratio, logfaagplotdata, labels=["arm milieu" "matig milieu" "rijk milieu"], xlabel = "lysogeny/lyse", ylabel="nr. faag/gem(nr. faag)", title="Voorkomen fagen"))

 bactplotdata = [overzicht[1:length(sit), 1]/mean(overzicht[1:length(sit), 1]),
  overzicht[length(sit) + 1:2*length(sit), 1]/mean(overzicht[length(sit) + 1:2*length(sit), 1]),
  overzicht[2*length(sit) + 1:3*length(sit), 1]/mean(overzicht[2*length(sit) + 1:3*length(sit), 1])]
 faagplotdata = [overzicht[1:length(sit), 2]/mean(overzicht[1:length(sit), 2]),
  overzicht[length(sit) + 1:2*length(sit), 2]/mean(overzicht[length(sit) + 1:2*length(sit), 2]),
  overzicht[2*length(sit) + 1:3*length(sit), 2]/mean(overzicht[2*length(sit) + 1:3*length(sit), 2])]

for l = 1:3
    for m = 1:7
        bactplotdata[l][m] = round(bactplotdata[l][m], digits = 5)
        faagplotdata[l][m] = round(faagplotdata[l][m], digits = 5)
    end
end

p = plot(plot(ratio, bactplotdata, labels=["arm milieu" "matig milieu" "rijk milieu"], xlabel = "lysogeny/lyse", ylabel="nr. bact/gem(nr. bact)", title="Voorkomen bacterien", legend=:topleft),
 plot(ratio, faagplotdata, labels=["arm milieu" "matig milieu" "rijk milieu"], xlabel = "lysogeny/lyse", ylabel="nr. faag/gem(nr. faag)", title="Voorkomen fagen", legend=:topleft))

simpars = @dict nsims nsteps
fname = savename(plotname, simpars)
savefig(p, plotsdir(fname) * ".png")
