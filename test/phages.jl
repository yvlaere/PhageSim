
@testset "updating phages" begin

    grid = rand(1:100, 50, 50)

    nphages = sum(grid)

    # check mass balances
    @test sum(updatephages!(grid)) == nphages


    @test_throws AssertionError updatephages!(grid, pdecay=1.2)
    @test sum(updatephages!(grid, pdecay=0.6)) < nphages
    @test sum(updatephages!(grid, poissonapprox=true, pdecay=0.1)) < nphages
    @test sum(updatephages!(grid, pdecay=1.0)) == 0
end
