  
using OpenKolmogorovFlow


@testset "radial mean                            " begin
    # create some data (with appropriate symmetries)
    Ω = FTField([ 0+0im -2-6im -1+0im;
                 -5+7im -0-2im  2+0im;
                  2+0im -1+2im  3+0im;
                 -5-7im  2+7im  2-0im])

    # mean of abs
    out = radial_mean(abs, Ω)

    @test_broken out[1] == mean(abs, [Ω[1, 1], Ω[1, -1], Ω[-1, 1], Ω[-1, -1],
                                     Ω[1, 0], Ω[0,  1], Ω[0, -1], Ω[-1,  0]])

    @test_broken out[2] == mean(abs, [Ω[-1,  2], Ω[ 0,  2], Ω[ 1,  2],  # top
                                     Ω[-1, -2], Ω[ 0, -2], Ω[ 1, -2],  # bottom
                                     Ω[-2,  1], Ω[-2,  0], Ω[-2,  1],  # left
                                     Ω[ 2,  1], Ω[ 2,  0], Ω[-2, -1]]) # right

    @test_broken out[3] == mean(abs, [Ω[2, 2], Ω[2, -2], Ω[-2, 2], Ω[-2, 2]])
end      