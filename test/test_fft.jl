@testset "allocating - some tests                " begin
    # make grid
    n, m = 4, 4
    x, y = make_grid(m)

    # constant field is reset by _apply_mask
    U = FTField(n, m)
    U[WaveNumber(0, 0)] = 1
    u = IFFT(U)
    @test all(u.data .== 0.0)

    # sin(y)
    U = FTField(n, m)
    U[WaveNumber( 1, 0)] =      -0.5im
    U[WaveNumber(-1, 0)] = conj(-0.5im)
    u = IFFT(U)
    @test maximum(abs, u.data - sin.(0.0*x .+ 1.0*y)) < 1e-15

    # cos(2y)
    U = FTField(n, m)
    U[WaveNumber( 2, 0)] =  0.5
    U[WaveNumber(-2, 0)] =  0.5
    u = IFFT(U)
    @test maximum(abs, u.data - cos.(0.0*x .+ 2.0*y)) < 1e-15

    # cos(2x)
    U = FTField(n, m)
    U[WaveNumber(0, 2)] = 0.5
    u = IFFT(U)
    @test maximum(abs, u.data - cos.(2.0*x .+ 0.0*y)) < 1e-15

    # cos(2x + 2y)
    U = FTField(n, m)
    U[WaveNumber(2, 2)] = 0.5
    u = IFFT(U)
    @test maximum(abs, u.data - cos.(2.0*x .+ 2.0*y)) < 1e-15
end

@testset "utils                                  " begin
    @testset "up                                  " begin
        # 0 1 x
        # 0 1 x
        @test up_dealias_size(1) == 1

        # 0 1 2 x
        # 0 1 2 3 x
        @test up_dealias_size(2) == 3

        # 0 1 2 3 x
        # 0 1 2 3 4 x
        @test up_dealias_size(3) == 4

        # 0 1 2 3 4 x
        # 0 1 2 3 4 5 6 x
        @test up_dealias_size(4) == 6

        # 0 1 2 3 4 5 x
        # 0 1 2 3 4 5 6 7 x
        @test up_dealias_size(5) == 7

        # 0 1 2 3 4 5 6 x
        # 0 1 2 3 4 5 6 7 8 9 x
        @test up_dealias_size(6) == 9
    end
    @testset "down                                  " begin
        # 0 1 x
        # 0 1 2 x
        @test down_dealias_size(2) == 1

        # 0 1 2 x
        # 0 1 2 3 x
        @test down_dealias_size(3) == 2

        # 0 1 2 3 x
        # 0 1 2 3 4 x
        @test down_dealias_size(4) == 3

        # 0 1 2 3 x
        # 0 1 2 3 4 5 x
        @test down_dealias_size(5) == 3

        # 0 1 2 3 4 x
        # 0 1 2 3 4 5 6 x
        @test down_dealias_size(6) == 4

        # 0 1 2 3 4 5 x
        # 0 1 2 3 4 5 6 7 x
        @test down_dealias_size(7) == 5

        # 0 1 2 3 4 5 x
        # 0 1 2 3 4 5 6 7 8 x
        @test down_dealias_size(8) == 5
    end
end