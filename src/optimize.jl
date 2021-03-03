function _welcome_message()
    welcome = """
    Coluna
    Version 0.3.6 | 2021-03-03 | https://github.com/atoptima/Coluna.jl
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
            Warning: No initial primal bound and no cost for global artificial variables.
            Cost of global artificial variables set to 100000.0
            """
            println(msg)
        end
    end
    if params.local_art_var_cost === nothing
        if init_pb != Inf && init_pb != -Inf
            exp = ceil(log(10, init_pb))
            params.local_art_var_cost = 10^exp
        else
            params.local_art_var_cost = 10000.0
            msg = """
            Warning: No initial primal bound and no cost for local artificial variables.
            Cost of local artificial variables set to 10000.0
            """
            println(msg)
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

    env = Env(params)

    # Apply decomposition
    reformulate!(prob, annotations, env)
    
    # Coluna ready to start
    set_optim_start_time!(env)
    @logmsg LogLevel(-1) "Coluna ready to start."
    @logmsg LogLevel(-1) env.params

    TO.@timeit _to "Coluna" begin
        optstate = optimize!(get_optimization_target(prob), env, init_pb, init_db)
    end

    env.kpis.elapsed_optimization_time = elapsed_optim_time(env)

    println(_to)
    TO.reset_timer!(_to)

    println("Terminated")
    println(string("Primal bound: ", get_ip_primal_bound(optstate)))
    println(string("Dual bound: ", get_ip_dual_bound(optstate)))
    return optstate, env.kpis
end

function optimize!(
    reform::MathProg.Reformulation, env::Env, initial_primal_bound, initial_dual_bound
)
    master = getmaster(reform)
    initstate = OptimizationState(
        master,
        ip_primal_bound = initial_primal_bound,
        ip_dual_bound = initial_dual_bound,
        lp_dual_bound = initial_dual_bound
    )

    algorithm = env.params.solver

    #this will initialize all the storages used by the algorithm and its child algorithms
    reformdata = Algorithm.ReformData(reform)
    Algorithm.initialize_storages!(reformdata, algorithm)

    output = Algorithm.run!(algorithm, env, reformdata, Algorithm.OptimizationInput(initstate))
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

function optimize!(
    form::MathProg.Formulation, env::Env, initial_primal_bound, initial_dual_bound
)
    initstate = OptimizationState(
        form,
        ip_primal_bound = initial_primal_bound,
        ip_dual_bound = initial_dual_bound,
        lp_dual_bound = initial_dual_bound
    )
    modeldata = Algorithm.ModelData(form)
    algorithm = env.params.solver
    output = Algorithm.run!(algorithm, env, modeldata, Algorithm.OptimizationInput(initstate))
    return Algorithm.getoptstate(output)
end

"""
Fallback if no solver provided by the user.
"""
function optimize!(::MathProg.Reformulation, ::Nothing, ::Real, ::Real)
    error("""
        No solver to optimize the reformulation. You should provide a solver through Coluna parameters. 
        Please, check the starting guide of Coluna.
    """)
end

function optimize!(::MathProg.Formulation, ::Nothing, ::Real, ::Real)
    error("""
        No solver to optimize the formulation. You should provide a solver through Coluna parameters. 
        Please, check the starting guide of Coluna.
    """)
end