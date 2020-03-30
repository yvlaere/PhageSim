@testset "Bacteria" begin
    bactgrid = emptybactgrid(10, 10)
    @test nbacteria(bactgrid) == 0

    bactgrid[3, 3] = Bacterium(1)
    bactgrid[4, 4] = Bacterium(1, 1)
    bactgrid[4, 6] = Bacterium(2)

    @test isbacterium(bactgrid[4, 4])
    @test !isbacterium(bactgrid[4, 5])
    @test nbacteria(bactgrid) == 3
    @test nbacteria(bactgrid, 1) == 2
    @test haslatent(bactgrid[4, 4])
    @test !haslatent(bactgrid[4, 6])
    @test prophage(bactgrid[4, 4]) == 1
    @test prophage(bactgrid[3,3]) === missing
    # creates new bacteria with a phage
    @test prophage(bactgrid[3, 3], 2) |> prophage == 2
    @test species(bactgrid[4, 4]) == 1
    @test Set(species(bactgrid)) == Set([1, 2])

    @test density(bactgrid) ≈ 3 / 100
    @test density(bactgrid, 1) ≈ 2 / 100
    @test density(bactgrid, CartesianIndex(8, 2)) == 0
    @test density(bactgrid, CartesianIndex(4, 4), R=1) ≈ 2 / 9
    @test density(bactgrid, CartesianIndex(4, 4), R=0) ≈ 1
    @test density(bactgrid, CartesianIndex(4, 4), R=2) ≈ 3 / 25
    @test density(bactgrid, CartesianIndex(4, 4), 1, R=2) ≈ 2 / 25

    bactrules = BacteriaRules(0.1, 0.3, 0.0)

end
