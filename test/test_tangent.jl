@testset "nonlinear vs tangent equations         " begin
    # seed random number generator
    Random.seed!(0)

    # take a simple setup
    n = 48
    m = up_dealias_size(n)
    Re = 20
    kforcing = 4

    # initial condition
    Ω = FTField(n, m)
    Ω[WaveNumber(1, 1)]  = 0.0 + 0.5*im
    Ω[WaveNumber(4, 3)]  = 0.1 + 0.1*im

    # get explicit and implicit parts
    EX_nl, IM_nl = splitexim(ForwardEquation(n, m, Re, kforcing))

    # define integrator
    ϕ = flow(EX_nl, IM_nl, CB3R2R3e(Ω, :NORMAL), TimeStepConstant(0.01))

    # Run for a while to land on attractor.
    ϕ(Ω, (0, 10))

    # now perturb Ω along a particular direction
    Ω[WaveNumber(1, 1)] += 1e-6 * (1 + im) 
    Ω_a = ϕ(copy(Ω), (0, 1))
    Ω[WaveNumber(1, 1)] -= 2e-6 * (1 + im)
    Ω_b = ϕ(copy(Ω), (0, 1))
    Ω[WaveNumber(1, 1)] += 1e-6 * (1 + im)

    # now take the difference between the two fields from the
    # nonlinear simulations and use that as reference. This
    # difference represents the evolution of the perturbation over
    # the reference trajectory.
    Ω_b .= (Ω_a.-Ω_b)./2e-6

    # Now define integration setup for variational equations
    EX_lin, IM_lin = splitexim(LinearisedEquation(n, m, Re, TangentMode()))

    # define initial condition, using the previous one
    Λ = FTField(n, m)

    # define integrator
    ψ = flow(couple(EX_nl, EX_lin),
             couple(IM_nl, IM_lin),
             CB3R2R3e(couple(Ω, Λ), :NORMAL),
             TimeStepConstant(0.01))
    
    # perturb mode (1, 1)
    Λ[WaveNumber(1, 1)] = 1 + im

    # integrate in time
    ψ(couple(Ω, Λ), (0, 1))

    # check that amplitude of mode (1, 1) of perturbation is similar
    # to the difference of two perturbed nonlinear simulations. We cannot
    # set too low a tolerance here, because nonlinear effects will kick in
    # at some point reducing the validity of the linearised equations.
    # this checks that we are basically not doing rubbish.
    a = Ω_b[WaveNumber(1, 1)]
    b = Λ[WaveNumber(1, 1)]

    @test abs(a-b)/abs(a) < 1e-8
end