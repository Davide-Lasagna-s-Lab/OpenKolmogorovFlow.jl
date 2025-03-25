# Test that solution converges to the
# laminar flow for small Reynolds numbers
@testset "laminar flow                           " begin
    # example dimension
    n = 10
    m = n

    # take a Re and a forcing wave number
    Re = 1.2345678
    kforcing = 4

    # initial condition
    Ω = FTField(n, m)

    # get explicit and implicit parts
    f = ForwardEquation(n, m, Re, kforcing)

    # define integrator
    ϕ = flow(splitexim(f)..., CB3R2R3e(Ω, :NORMAL), TimeStepConstant(0.01))

    # start from some non zero initial condition
    Ω .= 0.01; 
    Ω[WaveNumber(0, 0)] = 0

    # monitor the state
    mon = Monitor(Ω, copy)

    # map forward
    ϕ(Ω, (0, 50),  mon)

    # test final value is that predicted by explicit equation
    @test normdiff(samples(mon)[end], laminarflow(n, m, Re, kforcing)) < 1e-15
end