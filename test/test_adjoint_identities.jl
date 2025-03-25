@testset "adjoint                           " begin
    
    # seed rng
    Random.seed!(0)

    # setup
    n, m, Re, Δt = 31, 63, 40, 0.01

    # states
    Ω = FFT(Field(0.1*randn(2m+2, 2m+2)), n)
    #A = FFT(Field(0.1*randn(2m+2, 2m+2)), n)
    #B = FFT(Field(0.1*randn(2m+2, 2m+2)), n)
    A = FTField(n, m); A[WaveNumber(1, 2)] = 3 + 2*im
    B = FTField(n, m); B[WaveNumber(3, 4)] = 1 + 4*im

    # equations
    LD = LinearisedEquation(n, m, Re, TangentMode(), FFTW.ESTIMATE)
    LA = LinearisedEquation(n, m, Re, AdjointMode(), FFTW.ESTIMATE)

    v1 = dot(A, LD(0.0, Ω, copy(B), similar(Ω)))
    v2 = dot(B, LA(0.0, Ω, copy(A), similar(Ω)))
    @test abs(v1 - v2)/max(abs(v1), abs(v2)) < 4e-14


    # # We use the adjoint code to verify the identity 
    # #
    # # (Λ_T, L Ω′_0) = (L⁺Λ_T, Ω′_0)
    # # 
    # # where the subscripts _T and _0 denote the final and initial 
    # # states of the adjoint variables Λ and perturbation Ω′. The 
    # # operator L is the linearisation of the forward map, while L⁺ 
    # # is its adjoint, from integration by parts in space and time.
    # system right hand side
    F = ForwardEquation(n, m, Re, 4, FFTW.ESTIMATE)

    # flow
    ϕ = flow(splitexim(F)..., CNRK2(FTField(n, m), :NORMAL), TimeStepConstant(Δt))

    # define the stage cache
    cache = RAMStageCache(2, FTField(n, m))

    # proceed forward, then store forward solution for a small bit
    ϕ(Ω, (0, 10)); ϕ(Ω, (0, 10), reset!(cache))

    # construct linearised propagators
    ψ_D = flow(splitexim(LD)..., CNRK2(FTField(n, m), :TAN), TimeStepFromCache())
    ψ_A = flow(splitexim(LA)..., CNRK2(FTField(n, m), :ADJ), TimeStepFromCache())
    
    # verify identity
    v1 = dot(A, ψ_A(copy(B), copy(cache)))
    v2 = dot(B, ψ_D(copy(A), copy(cache)))

    @test abs(v1 - v2)/max(abs(v1), abs(v2)) < 4e-14
end