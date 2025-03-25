using OpenKolmogorovFlow
using IMEXRKCB
using BenchmarkTools

# parameters
const Re = 40
const kforcing = 4

# will contain the minimum times
tmin_d = Float64[]
tmin_a = Float64[]

# try these
ns = [32:2:256]

# loop
for n = ns
    # get explicit and implicit parts, for aliased and dealiased calculations
    Ld, Nd = imex(VorticityEquation(n, Re, kforcing; dealias=true))
    La, Na = imex(VorticityEquation(n, Re, kforcing; dealias=false))

    # initial condition
    Ω₀ = FTField(n)

    # define scheme
    scheme = IMEXRKScheme(IMEXRK4R3R(IMEXRKCB4, false), Ω₀)

    # measure time it takes to complete a step
    push!(tmin_d, minimum([@elapsed step!(scheme, Nd, Ld, 0, 0.1, Ω₀) for i = 1:50]))
    push!(tmin_a, minimum([@elapsed step!(scheme, Na, La, 0, 0.1, Ω₀) for i = 1:50]))

    # debug
    @printf "%3d -> %3d - %9.5f ms : %9.5f ms\n" n even_dealias_size(n) 1000*last(tmin_a) 1000*last(tmin_d)
end

using PyPlot; pygui(true)
loglog(ns, 10.0^(-7)*ns.^2.0.*log.(ns), dashes=Any[10, 2], label=L"n^2 log(n)")
loglog(ns, 10.0^(-6)*ns.^2,           dashes=Any[2,  2], label=L"n^2")
loglog(ns, tmin_a, ".-", label=L"aliased")
loglog(ns, tmin_d, "o-", label=L"dealiased")
legend(loc=2)
xlabel(L"n")
ylabel("time per step with IMEXRKCB4 [ms]")
show()