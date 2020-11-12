"""
    Problem

"""
mutable struct Problem <: AbstractProblem
    initial_primal_bound::Union{Nothing, Float64}
    initial_dual_bound::Union{Nothing, Float64}
    original_formulation::Formulation
    re_formulation::Union{Nothing, Reformulation}
    form_counter::Counter # 0 is for original form
    default_optimizer_builder::Function
end

"""
    Problem()

Constructs an empty `Problem`.
"""
function Problem()
    counter = Counter(-1)
    original_formulation = Formulation{Original}(counter)
    return Problem(
        nothing, nothing, original_formulation, nothing, counter,
        no_optimizer_builder
    )
end

set_original_formulation!(m::Problem, of::Formulation) = m.original_formulation = of
set_reformulation!(m::Problem, r::Reformulation) = m.re_formulation = r

get_original_formulation(m::Problem) = m.original_formulation
get_reformulation(m::Problem) = m.re_formulation

set_default_optimizer_builder!(p::Problem, default_opt_builder) = p.default_optimizer_builder = default_opt_builder

set_initial_primal_bound!(p::Problem, value::Real) = p.initial_primal_bound = value
set_initial_dual_bound!(p::Problem, value::Real) = p.initial_dual_bound = value

function get_initial_primal_bound(p::Problem)
    if p.original_formulation === nothing
        error("Cannot retrieve initial primal bound because the problem does not have original formulation.")
    end
    S = getobjsense(get_original_formulation(p))
    if p.initial_primal_bound !== nothing
        return PrimalBound{S}(p.initial_primal_bound)
    end
    return PrimalBound{S}()
end

function get_initial_dual_bound(p::Problem)
    if p.original_formulation === nothing
        error("Cannot retrieve initial dual bound because the problem does not have original formulation.")
    end
    S = getobjsense(get_original_formulation(p))
    if p.initial_dual_bound !== nothing
        return DualBound{S}(p.initial_dual_bound)
    end
    return DualBound{S}()
end
