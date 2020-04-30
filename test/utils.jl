@testset "utils" begin

    @testset "simplex" begin
        @test insimplex(0.1, 0.7, 0.2)
        @test !insimplex(0.1, 0.8, 0.2)
        @test !insimplex(-0.1, 0.8, 0.3)
    end

    @testset "summaries" begin
        bactgrid = emptybactgrid(10, 10)
        bactgrid[3, 3] = Bacterium(1)
        bactgrid[4, 4] = Bacterium(1, 1)
        bactgrid[4, 6] = Bacterium(2)
        bactgrid[1, 2] = Bacterium(2, 2)

        phagegrids = [ones(10, 10), zeros(10, 10)]

        bactcount, phagecount = species_counts(bactgrid, phagegrids; nspecies=2)
        @test bactcount == [2, 2]
        @test phagecount == [100, 0]

        bactcount, phagecount, prophagescount = species_counts_lat(bactgrid, phagegrids; nspecies=2)

        @test prophagescount == [1 0; 0 1]

    end
end
