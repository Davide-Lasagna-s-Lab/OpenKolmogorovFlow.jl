@testset "dot product                            " begin
    n, m = 5, 10
    x, y = make_grid(m)

    # tolerance on integrals
    TOL = 1e-14

    # init rng
    Random.seed!(0)

    # all cosine waves fit properly on the grid
    for j = -n:n, k =-n:n
        a, b = rand(), rand()
        u = a.*cos.(j.*x.+k.*y); U = FFT(Field(u), n)
        v = b.*cos.(j.*x.+k.*y); V = FFT(Field(v), n)
        j == k == 0 || @test abs(dot(U, V) - a*b/2) < TOL
    end

    # orthogonal fields have zero dot product
    u = sin.(1.0*x.+1.0*y); U = FFT(Field(u), n)
    v = sin.(1.0*x.+2.0*y); V = FFT(Field(v), n)
    @test abs(dot(U, V) - 0) < TOL

    u = cos.(5.0*x.+1.0*y) .+ sin.(4.0*x.+1.0*y); U = FFT(Field(u), n)
    v = cos.(1.0*x.+2.0*y) .+ sin.(1.0*x.+2.0*y); V = FFT(Field(v), n)
    @test abs(dot(U, V) - 0) < TOL

    # count only what is not orthogonal
    u = cos.(1.0*x.+0.0*y) .+ sin.(1.0*x.+1.0*y); U = FFT(Field(u), n)
    v = cos.(1.0*x.+0.0*y) .+ sin.(1.0*x.+2.0*y); V = FFT(Field(v), n)
    @test abs(dot(U, V) - 0.5) < TOL

    u = cos.(1.0*x.+0.0*y) .+ sin.(1.0*x.+1.0*y); U = FFT(Field(u), n)
    v = cos.(1.0*x.+0.0*y) .+ sin.(1.0*x.+1.0*y); V = FFT(Field(v), n)
    @test abs(dot(U, V) - 2*0.5) < TOL

    u = cos.(1.0*x.+5.0*y) .+ sin.(1.0*x.+5.0*y); U = FFT(Field(u), n)
    v = cos.(1.0*x.+5.0*y) .+ sin.(1.0*x.+5.0*y); V = FFT(Field(v), n)
    @test abs(dot(U, V) - 2*0.5) < TOL

    u = 0.4*cos.(5.0*x.+1.0*y) + 0.4*sin.(5.0*x.+2.0*y); U = FFT(Field(u), n)
    v = 0.3*cos.(5.0*x.+1.0*y) + 0.3*sin.(5.0*x.+2.0*y); V = FFT(Field(v), n)
    @test abs(dot(U, V) - 0.24*0.5) < TOL
end

@testset "dot product performance                " begin
    m, n = 49, 49
    U = FFT(Field(m, (x, y)->rand()), n)
    V = FFT(Field(m, (x, y)->rand()), n)
    @test minimum([@elapsed dot(U, V) for i = 1:100000]) < 8*10.0^(-6)
end

@testset "norm                                   " begin
    n, m = 10, 10
    x, y = make_grid(m)

    # tolerance on integrals
    TOL = 1e-14

    for (val, u) in [(sqrt(2)*π, sin.(1.0*x.+1.0*y)),
                     (2*π,       cos.(1.0*x.+2.0*y) .+ sin.(1.0*x.+2.0*y)),
                     (2*π,       cos.(0.0*x.+5.0*y) .+ sin.(1.0*x.+5.0*y)),
                     (2*π,       cos.(0.0*x.+5.0*y) .+ sin.(0.0*x.+5.0*y))]
        U = FFT(Field(u), n)
        @test abs(norm(U) - val/2π) < TOL
        @test_throws ArgumentError norm(U, 1)
    end
end

@testset "diff                                   " begin
    n, m = 5, 10
    x, y = make_grid(m)

    # tolerance on integrals
    TOL = 1e-14

    u = sin.(1.0*x.+1.0*y); U = FFT(Field(u), n)
    v = sin.(1.0*x.+1.0*y); V = FFT(Field(v), n)
    @test normdiff(U, V) == 0

    u = 2*sin.(1.0*x.+1.0*y); U = FFT(Field(u), n)
    v =   sin.(1.0*x.+1.0*y); V = FFT(Field(v), n)
    @test abs(normdiff(U, V) - 0.5) < TOL
end

@testset "minnormdiff                            " begin
    n, m = 41, up_dealias_size(41)
    U = FFT(Field(randn(2m+2, 2m+2)), n)
    V = copy(U)

    # number of shifts along x to test
    TMP = copy(U)
    N = 20

    # apply random shift to V, then min distance should be zero
    s = 3*2*π/20
    m = 2
    shift!(U, s, m)

    dmin, (smin, mmin) = minnormdiff(V, U, TMP, N)
    @test mmin == m
    @test smin == s
    @test normdiff(shift!(V, smin, mmin), U) < 1e-16
end