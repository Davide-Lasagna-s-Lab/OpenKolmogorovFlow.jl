@testset "constructors                           " begin
    @test_throws ArgumentError Field(randn(4, 5))
    @test_throws ArgumentError Field(randn(4, 6))
    @test_throws ArgumentError Field(randn(7, 7))

    u = Field(randn(6, 6))
    @test typeof(u) == Field{2, Float64, Matrix{Float64}}

    u = Field(4, Float32)
    @test typeof(u) == Field{4, Float32, Matrix{Float32}}

    u = Field(4, (x, y)->1)
    @test all(u.data .== 1)

    u = Field(4, (x, y)->sin(x))
    @test abs(u[1, 1]) < 2e-16
    @test abs(u[1, 6]) < 2e-16
end

@testset "indexing                               " begin
    @testset "getindex                               " begin
        u = Field([0  4  8 12
                   1  5  9 13
                   2  6 10 14
                   3  7 11 15])
        # lower left corner of domain, at
        # (x, y) = (1, 1)
        for i = 1:4, j = 1:4
            @test u[i, j] ==  parent(u)[i, j]
        end
        # linear indexing
        for i = 1:16
            @test u[i] == i-1
        end
    end
    @testset "setindex                               " begin
        u = Field([0  4  8 12
                   1  5  9 13
                   2  6 10 14
                   3  7 11 15])
        # lower left corner of domain, at
        # (x, y) = (1, 1)
        for i = 1:4, j = 1:4
            u[i, j] = i*j
            @test u[i, j] == i*j
        end
        # linear indexing
        for i = 1:16
            u[i] = i
            @test u[i] == i
        end
    end
end

@testset "similar/copy                           " begin
    u = Field(rand(6, 6)); u[1, 1] = 1
    v = similar(u)
    w = copy(u)
    @test typeof(v) == Field{2, Float64, Matrix{Float64}}
    @test w[1, 1] == 1
end