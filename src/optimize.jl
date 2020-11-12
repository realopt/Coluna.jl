function _welcome_message()
    welcome = """
    Coluna
    Version 0.3.2 - https://github.com/atoptima/Coluna.jl
    """
    print(welcome)
end

function _adjust_params(params, init_pb)
    if params.global_art_var_cost === nothing
        if init_pb != Inf && init_pb != -Inf
            exp = ceil(log(10, init_pb))
            params.global_art_var_cost = 10^(exp + 1)
        else
            params.global_art_var_cost = 100000.0
            msg = """
            No initial primal bound and no cost for global artificial variables.
            Cost of global artificial variables set to 100000.0
            """
            @warn(msg)
        end
    end
    if params.local_art_var_cost === nothing
        if init_pb != Inf && init_pb != -Inf
            exp = ceil(log(10, init_pb))
            params.local_art_var_cost = 10^exp
        else
            params.local_art_var_cost = 10000.0
            msg = """
            No initial primal bound and no cost for local artificial variables.
            Cost of local artificial variables set to 10000.0
            """
            @warn(msg)
        end
    end
    return
end

"""
Starting point of the solver.
"""
function optimize!(prob::MathProg.Problem, annotations::Annotations, params::Params)
    _welcome_message()

    # Adjust parameters
    ## Retrieve initial bounds on the objective given by the user
    init_pb = get_initial_primal_bound(prob)
    init_db = get_initial_dual_bound(prob)
    _adjust_params(params, init_pb)

    _set_global_params(params)

    # Apply decomposition
    reformulate!(prob, annotations)

    # Coluna ready to start
    _globals_.initial_solve_time = time()
    @logmsg LogLevel(-1) "Coluna ready to start."
    @logmsg LogLevel(-1) _params_

    relax_integrality!(prob.re_formulation.master) # TODO : remove

    TO.@timeit _to "Coluna" begin
        optstate = optimize!(
            prob.re_formulation, params.solver, init_pb, init_db
        )
    end
    #println(_to)
    open("timeroutputs.txt", "a") do io
        println(io, _to)
    end
    @logmsg LogLevel(0) "Terminated"
    @logmsg LogLevel(0) string("Primal bound: ", get_ip_primal_bound(optstate))
    @logmsg LogLevel(0) string("Dual bound: ", get_ip_dual_bound(optstate))
    return optstate
end

"""
Solve a reformulation
"""
function optimize!(
    reform::MathProg.Reformulation, algorithm::Algorithm.AbstractOptimizationAlgorithm,
    initial_primal_bound, initial_dual_bound
)

    master = getmaster(reform)
    initstate = OptimizationState(
        master,
        ip_primal_bound = initial_primal_bound,
        ip_dual_bound = initial_dual_bound,
        lp_dual_bound = initial_dual_bound
    )

    #this will initialize all the storages used by the algorithm and its child algorithms
    reformdata = Algorithm.ReformData(reform)
    Algorithm.initialize_storages!(reformdata, algorithm)

    output = Algorithm.run!(algorithm, reformdata, Algorithm.OptimizationInput(initstate))
    algstate = Algorithm.getoptstate(output)

    Algorithm.check_storage_states_participation(reformdata)

    # we copy optimisation state as we want to project the solution to the compact space
    outstate = OptimizationState(
        master,
        termination_status = getterminationstatus(algstate),
        ip_primal_bound = get_ip_primal_bound(algstate),
        ip_dual_bound = get_ip_dual_bound(algstate),
        lp_primal_bound = get_lp_primal_bound(algstate),
        lp_dual_bound = get_lp_dual_bound(algstate)
    )

    ip_primal_sols = get_ip_primal_sols(algstate)
    if ip_primal_sols !== nothing
        for sol in ip_primal_sols
            add_ip_primal_sol!(outstate, proj_cols_on_rep(sol, master))
        end
    end

    # lp_primal_sol may also be of interest, for example when solving the relaxation
    lp_primal_sol = get_best_lp_primal_sol(algstate)
    if lp_primal_sol !== nothing
        add_lp_primal_sol!(outstate, proj_cols_on_rep(lp_primal_sol, master))
    end

    return outstate
end
