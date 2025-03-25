@testset "test WaveNumberForcing         " begin

    # dimension
    n, m = 2, 4

    # fields
    U = FTField(n, m)
    V    = FTField(n, m)
    dVdt = FTField(n, m)

    # define grid
    x = range(0, stop=2π, length=2m+3)[1:2m+2]

    # test these ones
    tests = ( ((1, 0), cos.(1.0.*x .+ 0.0.*x'),  1),
              ((0, 1), cos.(0.0.*x .+ 1.0.*x'),  1),
              ((1, 1), cos.(1.0.*x .+ 1.0.*x'),  1),
              ((0, 1), sin.(0.0.*x .+ 1.0.*x'), -im),
              ((1, 0), sin.(1.0.*x .+ 0.0.*x'), -im),
              ((1, 1), sin.(1.0.*x .+ 1.0.*x'), -im))

    for ((k, j), expected, val) in tests
        # define
        f = WaveNumberForcing(n, k, j, val)

        # apply
        dVdt .= 0; f(0, U, V, dVdt)

        # check this is the same
        @test norm(vec(parent(IFFT(dVdt)) .- expected)) < 1e-12
    end
end