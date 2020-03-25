@testset "utils" begin

    @testset "simplex" begin
        @test insimplex(0.1, 0.7, 0.2)
        @test !insimplex(0.1, 0.8, 0.2)
        @test !insimplex(-0.1, 0.8, 0.3)
    end
end
