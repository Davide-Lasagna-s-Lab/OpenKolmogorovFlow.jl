@testset "perinterp                              " begin
    for (X, fun, gr, hs) in [([0.0+0.0*im, 0.5+0.0*im, 0.0+0.1*im, 0.0+0.0*im], t->cos(t) - 0.2*sin(2t),   t->-sin(t) - 0.4*cos(2t),  t->-cos(t) + 0.8*sin(2t)),
                             ([1.2+0.0*im, 0.0-0.2*im, 0.0+0.0*im, 0.0+0.0*im], t->1.2+0.4*sin(t),         t->0.4*cos(t),             t->-0.4*sin(t)),
                             ([2.0+0.0*im, 0.0-0.2*im, 0.0+0.5*im, 0.0+0.0*im], t->2.0+0.4*sin(t)-sin(2t), t->0.4*cos(t) - 2*cos(2t), t->-0.4*sin(t) + 4*sin(2t)),
                             ([0.0+0.0*im, 0.0+0.0*im, 0.0+0.0*im, 1.0+0.0*im], t->cos(3t),                t->-3*sin(3t),             t->-9*cos(3t))]
        for t = 0:0.1:2π
            val, grad, hess = OpenKolmogorovFlow._perinterp(X, t)
            @test val  ≈ fun(t)
            @test grad == gr(t)
            @test hess == hs(t)
        end
    end
end

@testset "peroptim                               " begin
    for (X, x_0, x_opt, f_opt) in [([0.0+0.0*im, 0.5+0.0*im, 0.0+0.0*im, 0.0+0.0*im], 0.1, 0.0, 1.0),
                                   ([0.0+0.0*im, 0.0-0.5*im, 0.0+0.0*im, 0.0+0.0*im], 1.5, π/2, 1.0)]
        x_opt_alg, f_opt_alg = OpenKolmogorovFlow._peroptim(X, x_0, 1e-8) 
        @test abs(x_opt - x_opt_alg) < 1e-8
        @test abs(f_opt - f_opt_alg) < 1e-16
    end
end

@testset "distance                               " begin
    @testset "same field                         " begin
        m, n = 31, down_dealias_size(31)
        u = Field(m, (x, y)->0.1*randn())
        U, V = FFT(u, n), FFT(u, n)

        # same cache size
        cache = DistanceCache(n)
        d, (s_opt, m_opt) = distance!(U, V, cache)
        @test abs(d)     < 5e-12
        @test abs(s_opt) < 5e-15
        @test m_opt == 0

        # reduced cache size
        cache = DistanceCache(11)
        d, (s_opt, m_opt) = distance!(U, V, cache)
        @test abs(d)     < 5e-12
        @test abs(s_opt) < 5e-15
        @test m_opt == 0
    end
#     @testset "known shift 1                      " begin
#         n = 64
#         x, y = make_grid(n)
#         u = Field(cos.(x .+ 0.*y)) # peak must be at grid point,
#         v = Field(sin.(x .+ 0.*y)) # else we get sampling artefacts
#         U, V = FFT(u), FFT(v)
#         # same cache size
#         cache = DistanceCache(n)
#         d, (s_opt, m_opt) = distance!(U, V, cache)
#         @test abs(d) < 5e-12
#         @test s_opt ≈ π/2
#         @test m_opt == 0
#         # reduced cache size
#         cache = DistanceCache(48)
#         d, (s_opt, m_opt) = distance!(U, V, cache)
#         @test abs(d) < 5e-12
#         @test s_opt ≈ π/2
#         @test m_opt == 0
#     end
#     @testset "known shift 2                      " begin
#         n, d = 64, 10
#         x, y = make_grid(n)
#         f = zeros(n, n)
#         for k = 0:d, j=0:d
#             if !( (k==d && j==0) || (k==d && j==d) || (k==0 && j==d) || (k==0 && j==0))
#                 f .+= cos.(j.*x .+ k.*y) .+ sin.(j.*x .+ k.*y)
#             end
#         end
#         U = FFT(Field(f))

#         # same cache size
#         cache = DistanceCache(n)

#         # test we recover all possible y shifts
#         for (s, m) in [(1.0, 0), (1.1, 2), (3, 4), (4.1, 6)]
#             V = shift!(deepcopy(U), (s, m))
#             d, (s_opt, m_opt) = distance!(U, V, cache)
#             @test abs(d) < 5e-12
#             @test s_opt ≈ s
#             @test m_opt == m
#         end
#     end
# end

# @testset "distance on fields                     " begin
#     # setup
#     Re, n, Δt = 40, 64, 0.015

#     # initial condition
#     srand(1)
#     Ω = laminarflow(n, Re)
#     for j = 1:5, k=1:5
#         Ω[k, j] = 0.1*(randn() + im*randn())
#     end

#     # integration schemes
#     RK = IMEXRKScheme(IMEXRK3R2R(IMEXRKCB3e, false), Ω)

#     # Get system
#     L, N = imex(ForwardEquation(n, Re; dealias=true))

#     # forward map
#     f  = integrator(N, L, RK,  Δt)

#     # run forward to get to steady state
#     f(Ω, (0, 100))

#     # distance cache
#     cache = DistanceCache(64)

#     for j = 0:70, m = 0:4
#         # shift exactly by a multiple of the grid size
#         Δ = (j*2π/64, 2*m)
#         Ωs = shift!(deepcopy(Ω), Δ)

#         # calculate distance
#         d, (s_opt, m_opt) = distance!(Ω, Ωs, cache)
#         @test abs(d) < 4e-10                  # the distance suffers cancellation
#         @test abs(s_opt - Δ[1] % 2π) < 1e-15  # this reaches machine accuracy
#         @test m_opt == 2*m % 8                # of course exact
#     end

#     # reduced distance cache size has some more error due to re-sampling.
#     cache = DistanceCache(56)

#     for j = 0:70, m = 0:4
#         # shift exactly by a multiple of the grid size
#         Δ = (j*2π/64, 2*m)
#         Ωs = shift!(deepcopy(Ω), Δ)

#         # calculate distance
#         d, (s_opt, m_opt) = distance!(Ω, Ωs, cache)
#         @test abs(d) < 1e-7                   # the distance suffers cancellation
#         @test abs(s_opt - Δ[1] % 2π) < 4e-10  # this reaches machine accuracy
#         @test m_opt == 2*m % 8                # of course exact
#     end

#     # proceed in time ahead
#     Ω2 = f(deepcopy(Ω), (0, 5))

#     # find shift and see if distance matches with the calculation using dotdiff
#     d, (s, m) = distance!(Ω, Ω2, cache)
#     @test abs(d - dotdiff(shifted(Ω, (s, m)), Ω2)) < 4e-9
end