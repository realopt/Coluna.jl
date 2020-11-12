getvalue(sol::Solution) = ColunaBase.getvalue(sol)
getvalue(bnd::Bound) = ColunaBase.getvalue(bnd)

mutable struct OptimizationState{F<:AbstractFormulation,S<:Coluna.AbstractSense}
    termination_status::TerminationStatus
    incumbents::ObjValues{S}
    max_length_ip_primal_sols::Int
    max_length_lp_primal_sols::Int
    max_length_lp_dual_sols::Int
    insert_function_ip_primal_sols::Function
    insert_function_lp_primal_sols::Function
    insert_function_lp_dual_sols::Function
    ip_primal_sols::Union{Nothing, Vector{PrimalSolution{F}}}
    lp_primal_sols::Union{Nothing, Vector{PrimalSolution{F}}}
    lp_dual_sols::Union{Nothing, Vector{DualSolution{F}}}
end

_sort!(sols::Vector{PrimalSolution{F}}, f::Function) where {F} = sort!(sols, by = x -> f(PrimalBound(x.model, x.value)))
_sort!(sols::Vector{DualSolution{F}}, f::Function) where {F} = sort!(sols, by = x -> f(DualBound(x.model, x.value)))

function bestbound(solutions::Vector{Sol}, max_len::Int, new_sol::Sol) where {Sol<:Solution}
    if length(solutions) < max_len
        push!(solutions, new_sol)
    else
        solutions[end] = new_sol
    end
    _sort!(solutions, MathProg.valueinminsense)
    return
end

"""
    OptimizationState(
        form; 
        termination_status = UNKNOWN_TERMINATION_STATUS,
        ip_primal_bound = nothing, 
        ip_dual_bound = nothing, 
        lp_primal_bound = nothing, 
        lp_dual_bound = nothing,
        max_length_ip_primal_sols = 1, 
        max_length_lp_dual_sols = 1, 
        max_length_lp_dual_sols = 1,
        insert_function_ip_primal_sols = bestbound, 
        insert_function_lp_primal_sols = bestbound, 
        insert_function_lp_dual_sols = bestbound
        )

A convenient structure to maintain and return solutions and bounds of a formulation `form` during an
optimization process. The termination statuses is considered as
unknown by default. You can define the initial incumbent bounds using `ip_primal_bound`,
`ip_dual_bound`, `lp_primal_bound`, and `lp_primal_bound` keyword arguments. Incumbent
bounds are set to infinite (according to formulation objective sense) by default.
You can store three types of solutions `ip_primal_sols`, `lp_primal_sols`, and `lp_dual_sols`.
These solutions are stored in three lists. Keywords `max_length_ip_primal_sols`,
`max_length_lp_primal_sols`, and `max_length_lp_dual_sols` let you define the maximum
size of the lists. Keywords `insert_function_ip_primal_sols`, `insert_function_lp_primal_sols`,
and `insert_function_lp_dual_sols` let you provide a function to define the way
you want to insert a new solution in each list. By default, lists are sorted by
best bound.
"""
function OptimizationState(
    form::F;
    termination_status::TerminationStatus = UNKNOWN_TERMINATION_STATUS,
    ip_primal_bound = nothing,
    ip_dual_bound = nothing,
    lp_primal_bound = nothing,
    lp_dual_bound = nothing,
    max_length_ip_primal_sols = 1,
    max_length_lp_primal_sols = 1,
    max_length_lp_dual_sols = 1,
    insert_function_ip_primal_sols = bestbound,
    insert_function_lp_primal_sols = bestbound,
    insert_function_lp_dual_sols = bestbound
) where {F <: AbstractFormulation}
    incumbents = ObjValues(
        form;
        ip_primal_bound = ip_primal_bound,
        ip_dual_bound = ip_dual_bound,
        lp_primal_bound = lp_primal_bound,
        lp_dual_bound = lp_dual_bound
    )
    S = getobjsense(form)
    state = OptimizationState{F,S}(
        termination_status, 
        incumbents,
        max_length_ip_primal_sols,
        max_length_lp_primal_sols,
        max_length_lp_dual_sols,
        insert_function_ip_primal_sols,
        insert_function_lp_primal_sols,
        insert_function_lp_dual_sols,
        nothing, nothing, nothing
    )
    return state
end

# function OptimizationState(
#     form::AbstractFormulation, or::OptimizationState
# )
#     newor = OptimizationState(
#         form,
#         termination_status = getterminationstatus(or),
#         ip_primal_bound = get_ip_primal_bound(or),
#         ip_dual_bound = get_ip_dual_bound(or),
#         lp_primal_bound = get_lp_primal_bound(or),
#         lp_dual_bound = get_lp_dual_bound(or)
#     )
#     if or.ip_primal_sols !== nothing
#         newor.ip_primal_sols = copy(or.ip_primal_sols)
#     end
#     if or.lp_primal_sols !== nothing
#         newor.lp_primal_sols = copy(or.lp_primal_sols)
#     end
#     if or.lp_primal_sols !== nothing
#         newor.lp_dual_sols = copy(or.lp_primal_sols)
#     end
#     return newor
# end

getform(state::OptimizationState{F,S}) where {F, S} = F

function OptimizationState(
    form::AbstractFormulation, source_state::OptimizationState, 
    copy_ip_primal_sol::Bool, copy_lp_primal_sol::Bool
)
    state = OptimizationState(
        form,
        termination_status = getterminationstatus(source_state),
        ip_primal_bound = get_ip_primal_bound(source_state),
        lp_primal_bound = PrimalBound(form),
        ip_dual_bound = get_ip_dual_bound(source_state),
        lp_dual_bound = get_lp_dual_bound(source_state)
    )
    if copy_ip_primal_sol && nb_ip_primal_sols(source_state) > 0
        set_ip_primal_sol!(state, copy(get_best_ip_primal_sol(source_state)))
    end
    if copy_lp_primal_sol && nb_lp_primal_sols(source_state) > 0
        set_lp_primal_sol!(state, copy(get_best_lp_primal_sol(source_state)))
    end
    return state
end

getterminationstatus(state::OptimizationState) = state.termination_status
#getsolutionstatus(state::OptimizationState) = state.solution_status

setterminationstatus!(state::OptimizationState, status::TerminationStatus) = state.termination_status = status
#setsolutionstatus!(state::OptimizationState, status::SolutionStatus) = state.solution_status = status

#isfeasible(state::OptimizationState) = state.solution_status == FEASIBLE_SOL
#isinfeasible(state::OptimizationState) = state.solution_status == INFEASIBLE_SOL

getincumbents(state::OptimizationState) = state.incumbents

get_ip_primal_bound(state::OptimizationState) = get_ip_primal_bound(state.incumbents)
get_lp_primal_bound(state::OptimizationState) = get_lp_primal_bound(state.incumbents)
get_ip_dual_bound(state::OptimizationState) = get_ip_dual_bound(state.incumbents)
get_lp_dual_bound(state::OptimizationState) = get_lp_dual_bound(state.incumbents)

update_ip_primal_bound!(state::OptimizationState, val) = update_ip_primal_bound!(state.incumbents, val)
update_ip_dual_bound!(state::OptimizationState, val) = update_ip_dual_bound!(state.incumbents, val)
update_lp_primal_bound!(state::OptimizationState, val) = update_lp_primal_bound!(state.incumbents, val)
update_lp_dual_bound!(state::OptimizationState, val) = update_lp_dual_bound!(state.incumbents, val)

set_ip_primal_bound!(state::OptimizationState, val) = set_ip_primal_bound!(state.incumbents, val)
set_lp_primal_bound!(state::OptimizationState, val) = set_lp_primal_bound!(state.incumbents, val)
set_ip_dual_bound!(state::OptimizationState, val) = set_ip_dual_bound!(state.incumbents, val)
set_lp_dual_bound!(state::OptimizationState, val) = set_lp_dual_bound!(state.incumbents, val)

ip_gap(state::OptimizationState) = ip_gap(state.incumbents)
lp_gap(state::OptimizationState) = lp_gap(state.incumbents)

ip_gap_closed(state::OptimizationState; kw...) = ip_gap_closed(state.incumbents; kw...)
lp_gap_closed(state::OptimizationState; kw...) = lp_gap_closed(state.incumbents; kw...)

function nb_ip_primal_sols(state::OptimizationState)
    return state.ip_primal_sols === nothing ? 0 : length(state.ip_primal_sols)
end

function nb_lp_primal_sols(state::OptimizationState)
    return state.lp_primal_sols === nothing ? 0 : length(state.lp_primal_sols)
end

function nb_lp_dual_sols(state::OptimizationState)
    return state.lp_dual_sols === nothing ? 0 : length(state.lp_dual_sols)
end

get_ip_primal_sols(state::OptimizationState) = state.ip_primal_sols

function get_best_ip_primal_sol(state::OptimizationState)
    nb_ip_primal_sols(state) == 0 && return nothing
    return get_ip_primal_sols(state)[1]
end

get_lp_primal_sols(state::OptimizationState) = state.lp_primal_sols

function get_best_lp_primal_sol(state::OptimizationState)
    nb_lp_primal_sols(state) == 0 && return nothing
    return get_lp_primal_sols(state)[1]
end

get_lp_dual_sols(state::OptimizationState) = state.lp_dual_sols

function get_best_lp_dual_sol(state::OptimizationState)
    nb_lp_dual_sols(state) == 0 && return nothing
    return get_lp_dual_sols(state)[1]
end

function clear_solutions!(state::OptimizationState)
    state.lp_primal_sols = nothing
    state.ip_primal_sols = nothing
    state.lp_dual_sols = nothing
end

function update_ip_primal!(
    dest_state::OptimizationState, orig_state::OptimizationState, set_solution::Bool
)
    update_ip_primal_bound!(dest_state, get_ip_primal_bound(orig_state))    
    if set_solution && nb_ip_primal_sols(orig_state) > 0
        set_ip_primal_sol!(dest_state, get_best_ip_primal_sol(orig_state))
    end
end

function update_all_ip_primal_solutions!(
    dest_state::OptimizationState, orig_state::OptimizationState
)
    # we do it in reverse order in order not to store the same solution several times
    if nb_ip_primal_sols(orig_state) > 0
        for ip_primal_sol in reverse(get_ip_primal_sols(orig_state))
            update_ip_primal_sol!(dest_state, ip_primal_sol)
        end    
    end    
end

function update!(dest_state::OptimizationState, orig_state::OptimizationState)
    setterminationstatus!(dest_state, getterminationstatus(orig_state))
    update_all_ip_primal_solutions!(dest_state, orig_state)
    update_ip_dual_bound!(dest_state, get_ip_dual_bound(orig_state))
    update_lp_dual_bound!(dest_state, get_lp_dual_bound(orig_state))
    set_lp_primal_bound!(dest_state, get_lp_primal_bound(orig_state))
    if nb_lp_primal_sols(orig_state) > 0
        set_lp_primal_sol!(dest_state, get_best_lp_primal_sol(orig_state))
    end        
end

"""
    update_ip_primal_sol!(optstate, sol)

Add the solution `sol` in the solutions list of `optstate` if and only if the 
value of the solution is better than the incumbent. The solution is inserted in the list
by the method defined in `insert_function_ip_primal_sols` field of `OptimizationState`.
If the maximum length of the list is reached, the solution located at the end of the list
is removed.
"""
function update_ip_primal_sol!(state::OptimizationState{F, S}, sol::PrimalSolution{F}) where {F, S}
    state.max_length_ip_primal_sols == 0 && return

    if state.ip_primal_sols === nothing
        state.ip_primal_sols = PrimalSolution{F}[]
    end 
    b = PrimalBound{S}(getvalue(sol))
    if update_ip_primal_bound!(state.incumbents, b)
        state.insert_function_ip_primal_sols(state.ip_primal_sols, state.max_length_ip_primal_sols, sol)
    end
    return
end

"""
    add_ip_primal_sol!(optstate, sol)

Add the solution `sol` in the solutions list of `opstate` and update the incumbent bound if 
the solution is better.
"""
function add_ip_primal_sol!(state::OptimizationState{F, S}, sol::PrimalSolution{F}) where {F, S}
    state.max_length_ip_primal_sols == 0 && return
    
    if state.ip_primal_sols === nothing
        state.ip_primal_sols = PrimalSolution{F}[]
    end
    state.insert_function_ip_primal_sols(state.ip_primal_sols, state.max_length_ip_primal_sols, sol)
    b = PrimalBound{S}(getvalue(sol))
    update_ip_primal_bound!(state.incumbents, b)
    return
end

"""
    set_ip_primal_sol!(optstate, sol)

Add the solution `sol` in the solutions list of `optstate`. The incumbent bound is not 
updated even if the value of the solution is better.
"""
function set_ip_primal_sol!(state::OptimizationState{F, S}, sol::PrimalSolution{F}) where {F, S}
    state.max_length_ip_primal_sols == 0 && return
    
    if state.ip_primal_sols === nothing
        state.ip_primal_sols = PrimalSolution{F}[]
    end
    state.insert_function_ip_primal_sols(state.ip_primal_sols, state.max_length_ip_primal_sols, sol)
    return
end

function update_lp_primal_sol!(state::OptimizationState{F, S}, sol::PrimalSolution{F}) where {F, S}
    state.max_length_lp_primal_sols == 0 && return

    if state.lp_primal_sols === nothing
        state.lp_primal_sols = PrimalSolution{F}[]
    end
    b = PrimalBound{S}(getvalue(sol))
    if update_lp_primal_bound!(state.incumbents, b)
        state.insert_function_lp_primal_sols(state.lp_primal_sols, state.max_length_lp_primal_sols, sol)
    end
    return
end

function add_lp_primal_sol!(state::OptimizationState{F, S}, sol::PrimalSolution{F}) where {F, S}
    state.max_length_lp_primal_sols == 0 && return
    
    if state.lp_primal_sols === nothing
        state.lp_primal_sols = PrimalSolution{F}[]
    end
    
    state.insert_function_lp_primal_sols(state.lp_primal_sols, state.max_length_lp_primal_sols, sol)
    b = PrimalBound{S}(getvalue(sol))
    update_lp_primal_bound!(state.incumbents, b)
    return
end

function set_lp_primal_sol!(state::OptimizationState{F, S}, sol::PrimalSolution{F}) where {F, S}
    state.max_length_lp_primal_sols == 0 && return

    if state.lp_primal_sols === nothing
        state.lp_primal_sols = PrimalSolution{F}[]
    end
    state.insert_function_lp_primal_sols(state.lp_primal_sols, state.max_length_lp_primal_sols, sol)
    return
end

function update_lp_dual_sol!(state::OptimizationState{F, S}, sol::DualSolution{F}) where {F, S}
    state.max_length_lp_dual_sols == 0 && return

    if state.lp_dual_sols === nothing
        state.lp_dual_sols = DualSolution{F}[]
    end
    b = DualBound{S}(getvalue(sol))
    if update_lp_dual_bound!(state.incumbents, b)
        state.insert_function_lp_dual_sols(state.lp_dual_sols, state.max_length_lp_dual_sols, sol)
    end
    return
end

function add_lp_dual_sol!(state::OptimizationState{F, S}, sol::DualSolution{F}) where {F, S}
    state.max_length_lp_dual_sols == 0 && return

    if state.lp_dual_sols === nothing
        state.lp_dual_sols = DualSolution{F}[]
    end
    state.insert_function_lp_dual_sols(state.lp_dual_sols, state.max_length_lp_dual_sols, sol)
    b = DualBound{S}(getvalue(sol))
    update_lp_dual_bound!(state.incumbents, b)
    return
end

function set_lp_dual_sol!(state::OptimizationState{F, S}, sol::DualSolution{F}) where {F, S}
    state.max_length_lp_dual_sols == 0 && return

    if state.lp_dual_sols === nothing
        state.lp_dual_sols = DualSolution{F}[]
    end
    state.insert_function_lp_dual_sols(state.lp_dual_sols, state.max_length_lp_dual_sols, sol)
    return
end


function Base.print(io::IO, form::AbstractFormulation, optstate::OptimizationState)
    println(io, "┌ Optimization state ")
    println(io, "│ Termination status: ", optstate.termination_status)
    println(io, "| Incumbents: ", optstate.incumbents)
    n = nb_ip_primal_sols(optstate)
    println(io, "| IP Primal solutions (",n,")")
    if n > 0
        for sol in optstate.ip_primal_sols
            print(io, form, sol)
        end
    end
    n = nb_lp_primal_sols(optstate)
    println(io, "| LP Primal solutions (",n,"):")
    if n > 0
        for sol in optstate.lp_primal_sols
            print(io, form, sol)
        end
    end
    n = nb_lp_dual_sols(optstate)
    println(io, "| LP Dual solutions (",n,"):")
    if n > 0
        for sol in optstate.lp_dual_sols
            print(io, form, sol)
        end
    end
    println(io, "└")
    return
end