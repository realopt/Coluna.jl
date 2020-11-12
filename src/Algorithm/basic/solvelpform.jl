"""
    Coluna.Algorithm.SolveLpForm(
        get_dual_solution = false,
        relax_integrality = false,
        set_dual_bound = false,
        silent = true
    )

Solve a linear program.
"""
@with_kw struct SolveLpForm <: AbstractOptimizationAlgorithm 
    update_ip_primal_solution = false
    consider_partial_solution = false
    get_dual_solution = false
    relax_integrality = false
    set_dual_bound = false
    silent = true
    log_level = 0
end

# SolveLpForm does not have child algorithms, therefore get_child_algorithms() is not defined

function get_storages_usage(
    algo::SolveLpForm, form::Formulation{Duty}
) where {Duty<:MathProg.AbstractFormDuty}
    # we use storages in the read only mode, as relaxing integrality
    # is reverted before the end of the algorithm, 
    # so the state of the formulation remains the same 
    storages_usage = Tuple{AbstractModel, StorageTypePair, StorageAccessMode}[] 
    push!(storages_usage, (form, StaticVarConstrStoragePair, READ_ONLY))
    if Duty <: MathProg.AbstractMasterDuty
        push!(storages_usage, (form, MasterColumnsStoragePair, READ_ONLY))
        push!(storages_usage, (form, MasterBranchConstrsStoragePair, READ_ONLY))
        push!(storages_usage, (form, MasterCutsStoragePair, READ_ONLY))
    end
    if algo.consider_partial_solution
        push!(storages_usage, (form, PartialSolutionStoragePair, READ_ONLY))
    end
    return storages_usage
end

function optimize_lp_form!(::SolveLpForm, optimizer, ::Formulation, ::OptimizationState) # fallback
    error("Cannot optimize LP formulation with optimizer of type ", typeof(optimizer), ".")
end

function optimize_lp_form!(
    algo::SolveLpForm, optimizer::MoiOptimizer, form::Formulation, result::OptimizationState
)
    MOI.set(form.optimizer.inner, MOI.Silent(), algo.silent)
    optimize_with_moi!(optimizer, form, result)
    return
end

function run!(algo::SolveLpForm, data::ModelData, input::OptimizationInput)::OptimizationOutput
    form = getmodel(data)
    result = OptimizationState(form)

    TO.@timeit Coluna._to "SolveLpForm" begin

    if algo.relax_integrality
        relax_integrality!(form)
    end

    partial_sol = nothing
    partial_sol_val = 0.0
    if algo.consider_partial_solution
        partsolstorage = getstorage(data, PartialSolutionStoragePair)
        partial_sol = get_primal_solution(partsolstorage, form)
        partial_sol_val = getvalue(partial_sol)
    end

    optimizer = getoptimizer(form)
    optimize_lp_form!(algo, optimizer, form, result)
    primal_sols = get_primal_solutions(form, optimizer)

    coeff = getobjsense(form) == MinSense ? 1.0 : -1.0

    if algo.get_dual_solution
        dual_sols = get_dual_solutions(form, optimizer)
        if length(dual_sols) > 0
            lp_dual_sol_pos = argmax(coeff * getvalue.(dual_sols))
            lp_dual_sol = dual_sols[lp_dual_sol_pos]
            set_lp_dual_sol!(result, lp_dual_sol)
            if algo.set_dual_bound
                db = DualBound(form, getvalue(lp_dual_sol) + partial_sol_val)
                set_lp_dual_bound!(result, db)
            end
        end
    end

    if length(primal_sols) > 0
        lp_primal_sol_pos = argmin(coeff * getvalue.(primal_sols))
        lp_primal_sol = primal_sols[lp_primal_sol_pos]
        add_lp_primal_sol!(result, lp_primal_sol)
        pb = PrimalBound(form, getvalue(lp_primal_sol) + partial_sol_val)
        set_lp_primal_bound!(result, pb)
        if algo.update_ip_primal_solution && isinteger(lp_primal_sol) && 
            !contains(lp_primal_sol, varid -> isanArtificialDuty(getduty(varid)))
            if partial_sol !== nothing
                add_ip_primal_sol!(result, cat(lp_primal_sol, partial_sol))
            else
                add_ip_primal_sol!(result, lp_primal_sol)
            end
        end
    end
    end # @timeit
    return OptimizationOutput(result)
end
