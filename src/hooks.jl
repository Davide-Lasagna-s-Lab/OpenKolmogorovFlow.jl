import Flows: AbstractTimeStepFromHook, Coupled

export CFLHook

"""
    CFLHook([CFL_max=1, Δt_min=0, Δt_max=1])

Time-step selection hook based on the CFL estimate stored by the explicit term.

The hook returns a clamped time step in the interval `(Δt_min, Δt_max)`. For a
single [`ForwardExplicitTerm`](@ref), it reads the most recent advective time
scale from `β[1]`; for a `Flows.Coupled` object, it uses the first equation.
"""
struct CFLHook <: AbstractTimeStepFromHook
       CFL_max::Float64
       Δt_bounds::Tuple{Float64, Float64}
    function CFLHook(CFL_max::Real=1, Δt_min::Real=0, Δt_max::Real=1)
        CFL_max > 0 || throw(ArgumentError("time step must be positive"))
        new(Float64(CFL_max), (Float64(Δt_min), Float64(Δt_max)))
    end
end

# Direct forward integration: the explicit term updates `β[1]` whenever the
# nonlinear right-hand side is evaluated, so the hook can cheaply reuse it here.
(hook::CFLHook)(g::ForwardExplicitTerm, A, z) = 
    clamp(hook.CFL_max * g.β[1], hook.Δt_bounds...)

# integrate coupled equations based on CFL of nonlinear state
(hook::CFLHook)(g::Coupled, A, z) = 
    clamp(hook.CFL_max * g[1].β[1], hook.Δt_bounds...)
