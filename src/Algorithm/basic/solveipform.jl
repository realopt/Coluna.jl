"""
    Coluna.Algorithm.SolveIpForm(
        time_limit::Int = 600,
        deactivate_artificial_vars = true,
        enforce_integrality = true,
        silent = true,
        log_level = 0
    )

Solve a mixed integer linear program.
"""
@with_kw struct SolveIpForm <: AbstractOptimizationAlgorithm
    time_limit::Int = 600
    deactivate_artificial_vars = true
    enforce_integrality = true
    get_dual_bound = true
    silent = true
    log_level = 0
end

# SolveIpForm does not have child algorithms, therefore get_child_algorithms() is not defined

function get_storages_usage(
    algo::SolveIpForm, form::Formulation{Duty}
) where {Duty<:MathProg.AbstractFormDuty}
    # we use storages in the read only mode, as all modifications
    # (deactivating artificial vars and enforcing integrality)
    # are reverted before the end of the algorithm,
    # so the state of the formulation remains the same
    storages_usage = Tuple{AbstractModel, StorageTypePair, StorageAccessMode}[] 
    push!(storages_usage, (form, StaticVarConstrStoragePair, READ_ONLY))
    if Duty <: MathProg.AbstractMasterDuty
        push!(storages_usage, (form, PartialSolutionStoragePair, READ_ONLY))
        push!(storages_usage, (form, MasterColumnsStoragePair, READ_ONLY))
        push!(storages_usage, (form, MasterBranchConstrsStoragePair, READ_ONLY))
        push!(storages_usage, (form, MasterCutsStoragePair, READ_ONLY))
    end
    return storages_usage
end

get_storages_usage(algo::SolveIpForm, reform::Reformulation) =
    get_storages_usage(algo, getmaster(reform))

function check_if_optimizer_supports_ip(optimizer::MoiOptimizer)
    return MOI.supports_constraint(optimizer.inner, MOI.SingleVariable, MOI.Integer)
end
check_if_optimizer_supports_ip(optimizer::UserOptimizer) = true
    
function run!(algo::SolveIpForm, data::ModelData, input::OptimizationInput)::OptimizationOutput
    form = getmodel(data)

    result = OptimizationState(
        form, ip_primal_bound = get_ip_primal_bound(getoptstate(input))
    )

    ip_supported = check_if_optimizer_supports_ip(getoptimizer(form))
    if !ip_supported
        @warn "Optimizer of formulation with id =", getuid(form),
              " does not support integer variables. Skip SolveIpForm algorithm."
        setterminationstatus!(result, UNKNOWN_TERMINATION_STATUS)
        return OptimizationOutput(result)
    end

    primal_sols = optimize_ip_form!(algo, getoptimizer(form), form, result)

    partial_sol = nothing
    partial_sol_value = 0.0
    if isa(form, Formulation{MathProg.DwMaster})
        partsolstorage = getstorage(data, PartialSolutionStoragePair)
        partial_sol = get_primal_solution(partsolstorage, form)
        partial_sol_value = getvalue(partial_sol)
    end
    
    if length(primal_sols) > 0
        coeff = getobjsense(form) == MinSense ? 1.0 : -1.0
        bestprimalsol_pos = argmin(coeff * getvalue.(primal_sols))
        bestprimalsol = primal_sols[bestprimalsol_pos]

        if partial_sol !== nothing
            add_ip_primal_sol!(result, cat(partial_sol, bestprimalsol))
        else
            add_ip_primal_sol!(result, bestprimalsol)
        end
        if algo.log_level == 0
            @printf "Found primal solution of %.4f \n" getvalue(get_ip_primal_bound(result))
        end
        @logmsg LogLevel(-3) get_best_ip_primal_sol(result)
    else
        if algo.log_level == 0
            println(
                "No primal solution found. Termination status is ",
                getterminationstatus(result), ". "
            )
        end
    end
    if algo.get_dual_bound && getterminationstatus(result) == OPTIMAL
        dual_bound = getvalue(get_ip_primal_bound(result)) + partial_sol_value
        set_ip_dual_bound!(result, DualBound(form, dual_bound))
    end
    return OptimizationOutput(result)
end

run!(algo::SolveIpForm, data::ReformData, input::OptimizationInput) = 
    run!(algo, getmasterdata(data), input)

function termination_status!(result::OptimizationState, optimizer::MoiOptimizer)
    terminationstatus = MOI.get(getinner(optimizer), MOI.TerminationStatus())
    if terminationstatus != MOI.INFEASIBLE &&
            terminationstatus != MOI.DUAL_INFEASIBLE &&
            terminationstatus != MOI.INFEASIBLE_OR_UNBOUNDED &&
            terminationstatus != MOI.OPTIMIZE_NOT_CALLED &&
            terminationstatus != MOI.TIME_LIMIT

        setterminationstatus!(result, convert_status(terminationstatus))

        if MOI.get(getinner(optimizer), MOI.ResultCount()) <= 0
            msg = """
            Termination status = $(terminationstatus) but no results.
            Please, open an issue at https://github.com/atoptima/Coluna.jl/issues
            """
            error(msg)
        end
    else
        @warn "Solver has no result to show."
        setterminationstatus!(result, INFEASIBLE)
    end
    return
end

function optimize_with_moi!(optimizer::MoiOptimizer, form::Formulation, result::OptimizationState)
    sync_solver!(optimizer, form)
    nbvars = MOI.get(form.optimizer.inner, MOI.NumberOfVariables())
    if nbvars <= 0
        @warn "No variable in the formulation. Coluna does not call the solver."
    else
        MOI.optimize!(getinner(optimizer))
    end
    termination_status!(result, optimizer)
    return
end

function optimize_ip_form!(
    algo::SolveIpForm, optimizer::MoiOptimizer, form::Formulation, result::OptimizationState
)
    MOI.set(optimizer.inner, MOI.TimeLimitSec(), algo.time_limit)
    MOI.set(optimizer.inner, MOI.Silent(), algo.silent)
    # No way to enforce upper bound or lower bound through MOI.
    # Add a constraint c'x <= UB in form ?

    if algo.deactivate_artificial_vars
        deactivate!(form, vcid -> isanArtificialDuty(getduty(vcid)))
    end
    if algo.enforce_integrality
        enforce_integrality!(form)
    end

    optimize_with_moi!(optimizer, form, result)
    primal_sols = get_primal_solutions(form, optimizer)

    if algo.enforce_integrality
        relax_integrality!(form)
    end
    if algo.deactivate_artificial_vars
        activate!(form, vcid -> isanArtificialDuty(getduty(vcid)))
    end
    return primal_sols
end

function optimize_ip_form!(
    ::SolveIpForm, optimizer::UserOptimizer, form::Formulation, result::OptimizationState
)
    @logmsg LogLevel(-2) "Calling user-defined optimization function."

    cbdata = MathProg.PricingCallbackData(form)
    optimizer.user_oracle(cbdata)

    if length(cbdata.primal_solutions) > 0
        setterminationstatus!(result, OPTIMAL)
    else
        setterminationstatus!(result, INFEASIBLE) # TODO : what if no solution found ?
    end
    return cbdata.primal_solutions
end
