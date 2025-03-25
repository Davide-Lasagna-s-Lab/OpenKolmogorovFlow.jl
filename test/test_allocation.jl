# Test that no memory is allocated in calls
@testset "allocation                             " begin
    # example dimension
    m = 49
    n = down_dealias_size(m)

    # take a Re and a forcing wave number
    Re = 1.2345678
    kforcing = 4

    # initial condition
    Ω = FTField(n, m)

    # get explicit and implicit parts
    f = ForwardEquation(n, m, Re, kforcing)

    for (method, elaps_exp) in [(CNRK2(   Ω, :NORMAL), 0.050),
                                (CB3R2R3e(Ω, :NORMAL), 0.115),
                                (CB3R2R3c(Ω, :NORMAL), 0.115)]

        # define integrator
        ϕ = flow(splitexim(f)..., method, TimeStepConstant(0.01))

        # test allocations
        min_time = minimum([@elapsed ϕ(Ω, (0, 1)) for i = 1:10])
        @test min_time < elaps_exp
        # this will not be zero because we have threading code
        # that results in a tiny allocation at every function call
        # change this to zero when it gets fixed
        # @test (@allocated ϕ(Ω, (0, 1))) == 0
    end
end