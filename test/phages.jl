
@testset "updating phages" begin

    grid = rand(1:100, 50, 50)

    nphages = sum(grid)

    # check mass balances
    @test sum(updatephages!(grid, PhageRules(0.0))) == nphages


    @test_throws AssertionError updatephages!(grid, PhageRules(-0.01))
    @test sum(updatephages!(grid, PhageRules(0.6))) < nphages
    @test sum(updatephages!(grid, PhageRules(0.1), poissonapprox=true)) < nphages
end
