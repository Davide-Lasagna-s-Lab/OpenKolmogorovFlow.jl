@testset "dissipation                            " begin
    # test on laminar flow
    Re = 1.234567
    kf = 4
    @test_throws ArgumentError laminarflow(2, 4, Re, 6)
    
    Ω = laminarflow(10, 20, Re, kf)
    @test dissrate(Ω, Re) ≈ Re/2/kf^2

    # phase shift does not change dissipation
    Ωa = FFT(Field(10, (x, y)->cos(x+y)), 5)
    Ωb = FFT(Field(10, (x, y)->sin(x+y)), 5)
    @test abs(dissrate(Ωa, 1.0) - dissrate(Ωb, 1.0)) < 4e-16

    # for any wave dissipation is 1/2/Re
    d = 5
    for j = -d:d, k=-d:d
        if !(j == 0 && k == 0)
            Re = randn()
            Ω = FFT(Field(10, (x, y)->cos(j*x+k*y)), d)
            @test dissrate(Ω, Re) ≈ 1/2/Re
        end
    end

    # test no allocations
    @test (@allocated dissrate(Ω, 1.0)) == 16

    # warm up with int
    dissrate(Ω, 1)
    @test (@allocated dissrate(Ω, 1)) == 16
end