@testset "norm is shift invariant                " begin
    n, m = 41, up_dealias_size(41)
    s, _m = 0.123, 3
    U = FFT(Field(m, (x, y) -> randn()), n)
    V = shift!(copy(U), s, _m)
    @test norm(U) ≈ norm(V)
end

@testset "shift cos right by π/2 get a sine      " begin
    m, n = 40, 40
    Ucos = FFT(Field(m, (x, y)->cos(x)), n)
    Usin = FFT(Field(m, (x, y)->sin(x)), n)

    @test normdiff(Usin, xshift!(Ucos, -π/2)) < 1e-20
end

@testset "do it up, by π/2                       " begin
    n, m = 41, up_dealias_size(41)
    Ucos = FFT(Field(m, (x, y)->cos(y)), n)
    Usin = FFT(Field(m, (x, y)->sin(y)), n)

    @test normdiff(Usin, yshift!(Ucos, -1)) < 1e-20
end