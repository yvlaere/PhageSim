@testset "Bacteria" begin
    bactgrid = emptybactgrid(10, 10)
    @test nbacteria(bactgrid) == 0

    bactgrid[3, 3] = Bacterium(1)
    bactgrid[4, 4] = Bacterium(1, 1)
    bactgrid[4, 6] = Bacterium(2)

    @test isbacterium(bactgrid[4, 4])
    @test !isbacterium(bactgrid[4, 5])
    @test nbacteria(bactgrid) == 3
    @test haslatent(bactgrid[4, 4])
    @test !haslatent(bactgrid[4, 6])
    @test phage(bactgrid[4, 4]) == 1
    @test species(bactgrid[4, 4]) == 1

    bactrules = BacteriaRules(0.1, 0.3, 0.0)

end
