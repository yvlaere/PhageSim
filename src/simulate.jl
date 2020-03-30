#=
Created on Wednesday 26 March 2020
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

Functionality to simulate bacterium-phage interactions.
=#

export step_bacteria!, simulate!


"""
    step_bacteria!(grid::BactGrid,
                    phagegrids,
                    bactrules::AbstractBacteriaRules,
                    phagerules::AbstractPhageRules,
                    interactionrules::AbstractInteractionRules)

Perform a single step in the bacteria-phage model.
"""
function step_bacteria!(grid::BactGrid,
                phagegrids,
                bactrules::AbstractBacteriaRules,
                phagerules::AbstractPhageRules,
                interactionrules::AbstractInteractionRules)
    # number of phages
    nphagetypes = length(phagegrids)
    # list of phage types, used for effecient shuffling
    phagetypes = collect(1:nphagetypes)
    bactcoors = Set(findall(isbacterium, grid))
    nbacts = length(bactcoors)
    C = CartesianIndices(grid)
    Ifirst, Ilast = first(C), last(C)
    IR = oneunit(Ifirst)
    for i in 1:nbacts
        I = rand(bactcoors)
        sI = bact = grid[I]  # state at position I
        delete!(bactcoors, I)
        # first check if the bacterium has a laten phage
        if haslatent(bact)
            # if it has a latent phage, it might lyse now
            if lyses(bact, grid, I, interactionrules::AbstractInteractionRules)
                # bacterium lyses, bye!
                phagetype = prophage(bact)
                phagegrids[phagetype][I] += burstsize(bact, phagetype, interactionrules)
                grid[I] = nothing
                sI = nothing
                continue
            end
        else
            # no latent phage, so check if there are any phages that might attack
            for phagetype in randperm!(phagetypes)
                nphages = phagegrids[phagetype][I]
                # if there are phages, will they infect the bacterium?
                if nphages > 0 && infects(bact, phagetype, nphages, interactionrules)
                    # yes, but lysogentic or lytic?
                    if lysogenic(bact, phagetype, interactionrules)
                        bact = prophage(bact, phagetype)  # add a prophage
                        grid[I] = bact
                        phagegrids[phagetype][I] -= 1
                    else
                        grid[I] = nothing  # kill bacterium
                        phagegrids[phagetype][I] += burstsize(bact, phagetype, interactionrules) - 1
                        break
                    end
                end
            end
        end
        isbacterium(grid[I]) || continue  # if killed, skip the rest
        # determine action
        neighborhood = max(Ifirst, I-IR):min(Ilast, I+IR)
        N = randnonident(I, neighborhood)
        sN = grid[N]
        sI, sN = updatebact(sI, sN, bactrules)
        grid[I], grid[N] = sI, sN
    end
end

"""
    simulate!(grid::BactGrid,
                    phagegrids,
                    bactrules::AbstractBacteriaRules,
                    phagerules::AbstractPhageRules,
                    interactionrules::AbstractInteractionRules,
                    nsteps::Int;
                    resultfun::union{Function,Nothing}=nothing,
                    store_every::Int=1
                    )

Perform a simulation of `nsteps`. Optionally provide a function `resultfun`
with as inputs the bacteria grid and the phagegrids to store intermediate
results every `store_every` steps.
"""
function simulate!(grid::BactGrid,
                phagegrids,
                bactrules::AbstractBacteriaRules,
                phagerules::AbstractPhageRules,
                interactionrules::AbstractInteractionRules,
                nsteps::Int;
                resultfun::Union{Function,Nothing}=nothing,
                store_every::Int=1
                )
    # should we store any results?
    if !isnothing(resultfun)
        # compute the first result
        res0 = resultfun(grid, phagegrids)
        results = Array{typeof(res0)}(undef, nsteps ÷ store_every + 1)
        results[1] = res0
    end
    for t in 1:nsteps
        # update the bacteria
        step_bacteria!(grid, phagegrids, bactrules, phagerules, interactionrules)
        # update the phages
        for phagegrid in phagegrids
            step_phages!(phagegrid, phagerules)
        end
        !isnothing(resultfun) && (t % store_every == 0) && (
                results[t ÷ store_every + 1] = resultfun(grid, phagegrids)
        )
    end
    !isnothing(resultfun) && return results
end
