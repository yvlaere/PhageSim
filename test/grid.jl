

@testset "grid" begin

    A = ones(4, 5)
    C = CartesianIndices(A)
    Ifirst, Ilast = first(C), last(C)

    @testset "periodic" begin
        # inbounds
        I = CartesianIndex(2, 3)
        @test periodic(I, Ifirst, Ilast) == I

        #out
        I = CartesianIndex(1, 6)
        @test periodic(I, Ifirst, Ilast) == CartesianIndex(1, 1)

        I = CartesianIndex(0, 6)
        @test periodic(I, Ifirst, Ilast) == CartesianIndex(4, 1)

        I = CartesianIndex(0, 2)
        @test periodic(I, Ifirst, Ilast) == CartesianIndex(4, 2)
    end

    @testset "neigbors" begin
        N = neighbors(A)

        @test size(N) == size(A)  # eight neigbors for each cell
        @test size(N[1]) == (8,)

        I = CartesianIndex(2, 3)

        @test I ∉ N[I]
    end

    @testset "regions" begin
        R = regions(A)

        @test size(R) == size(A)  # eight neigbors for each cell
        @test size(R[1]) == (9,)

        I = CartesianIndex(2, 3)

        @test I ∈ R[I]
    end



end
