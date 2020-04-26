"""
    SolveLpForm

todo
"""
Base.@kwdef struct SolveLpForm <: AbstractOptimizationAlgorithm 
    get_dual_solution = false
    relax_integrality = false
    set_dual_bound = false
    warm_start = false
    log_level = 0
end

# struct MasterLpRecord <: AbstractAlgorithmResult
#     incumbents::Incumbents
#     proven_infeasible::Bool
# end

# function prepare!(algo::MasterLp, form, node)
#     @logmsg LogLevel(-1) "Prepare MasterLp."
#     return
# end

function run!(algo::SolveLpForm, form::Formulation, input::OptimizationInput)::OptimizationOutput
    optstate = OptimizationState(form)

    if algo.relax_integrality
        relax_integrality!(form)
    end

    if algo.warm_start && nb_lp_primal_sols(optstate) > 0
        setwarmstart!(form, getoptimizer(form), get_best_lp_primal_sol(optstate))
    end

    optimizer_result = optimize!(form)

    setfeasibilitystatus!(optstate, getfeasibilitystatus(optimizer_result))    
    setterminationstatus!(optstate, getterminationstatus(optimizer_result))   

    lp_primal_sol = getbestprimalsol(optimizer_result)
    if lp_primal_sol !== nothing
        add_lp_primal_sol!(optstate, lp_primal_sol)
        if isinteger(lp_primal_sol) && !contains(lp_primal_sol, varid -> isanArtificialDuty(getduty(varid)))
            add_ip_primal_sol!(optstate, lp_primal_sol)
        end
    end

    if algo.get_dual_solution
        lp_dual_sol = getbestdualsol(optimizer_result)
        if lp_dual_sol !== nothing
            if algo.set_dual_bound
                update_lp_dual_sol!(optstate, lp_dual_sol)
            else
                set_lp_dual_sol!(optstate, lp_dual_sol)
            end
        end
    end

    return OptimizationOutput(optstate)
end
