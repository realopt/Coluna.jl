Base.@kwdef struct ColumnGeneration <: AbstractOptimizationAlgorithm
    max_nb_iterations::Int = 1000
    optimality_tol::Float64 = 1e-5
    log_print_frequency::Int = 1
    store_all_ip_primal_sols::Bool = false
end

# Data stored while algorithm is running
mutable struct ColGenRuntimeData
    incumbents::OptimizationState
    has_converged::Bool
    is_feasible::Bool
    ip_primal_sols::Vector{PrimalSolution}
    phase::Int64
end

function ColGenRuntimeData(
    algparams::ColumnGeneration, form::Reformulation, ipprimalbound::PrimalBound
)
    inc = OptimizationState(getmaster(form))
    set_ip_primal_bound!(inc, ipprimalbound)
    return ColGenRuntimeData(inc, false, true, [], 2)
end

function run!(algo::ColumnGeneration, reform::Reformulation, input::NewOptimizationInput)::OptimizationOutput    
    input_result = getinputresult(input)
    data = ColGenRuntimeData(algo, reform, get_ip_primal_bound(input_result))

    cg_main_loop!(algo, data, reform)
    masterform = getmaster(reform)
    if should_do_ph_1(masterform, data)
        set_ph_one(masterform, data)        
        cg_main_loop!(algo, data, reform)
        # TO DO : to implement unsetting phase one !!
        # TO DO : to implement repeating of phase two in the case phase 1 succeded
    end

    if data.is_feasible
        @logmsg LogLevel(-1) "ColumnGeneration terminated with status FEASIBLE."
    else
        data.incumbents = OptimizationState(getmaster(reform))
        setfeasibilitystatus!(data.incumbents, INFEASIBLE)
        @logmsg LogLevel(-1) "ColumnGeneration terminated with status INFEASIBLE."
    end

    if !algo.store_all_ip_primal_sols && nb_ip_primal_sols(data.incumbents) > 0
        for ip_primal_sol in get_ip_primal_sols(data.incumbents)
            push!(data.ip_primal_sols, ip_primal_sol)
        end
    end

    result = OptimizationState(
        masterform, 
        feasibility_status = data.is_feasible ? FEASIBLE : INFEASIBLE,
        termination_status = data.has_converged ? OPTIMAL : OTHER_LIMIT,
        ip_primal_bound = get_ip_primal_bound(data.incumbents),
        ip_dual_bound = get_lp_dual_bound(data.incumbents), # TODO : check if objective function is integer
        lp_dual_bound = get_lp_dual_bound(data.incumbents)
    )

    # add primal sols (data.ip_primal_sols)
    for ip_primal_sol in data.ip_primal_sols
        add_ip_primal_sol!(result, ip_primal_sol)
    end
    
    if nb_lp_primal_sols(data.incumbents) > 0
        for lp_primal_sol in get_lp_primal_sols(data.incumbents)
            add_lp_primal_sol!(result, lp_primal_sol)
        end
    end
    return OptimizationOutput(result)
end

# Internal methods to the column generation
function should_do_ph_1(master::Formulation, data::ColGenRuntimeData)
    ip_gap(data.incumbents) <= 0.00001 && return false
    primal_lp_sol = get_lp_primal_sols(data.incumbents)[1]
    if contains(primal_lp_sol, vid -> isanArtificialDuty(getduty(vid)))
        @logmsg LogLevel(-2) "Artificial variables in lp solution, need to do phase one"
        return true
    else
        @logmsg LogLevel(-2) "No artificial variables in lp solution, will not proceed to do phase one"
        return false
    end
end

function set_ph_one(master::Formulation, data::ColGenRuntimeData)
    for (varid, var) in getvars(master)
        isanArtificialDuty(getduty(varid)) && continue
        setcurcost!(master, varid, 0.0)
    end
    data.phase = 1
    return
end

function update_pricing_problem!(spform::Formulation, dual_sol::DualSolution)
    masterform = getmaster(spform)
    for (varid, var) in getvars(spform)
        iscuractive(spform, varid) || continue
        getduty(varid) <= AbstractDwSpVar || continue
        setcurcost!(spform, var, computereducedcost(masterform, varid, dual_sol))
    end
    return false
end

function update_pricing_target!(spform::Formulation)
    # println("pricing target will only be needed after automating convexity constraints")
end

function insert_cols_in_master!(
    masterform::Formulation, spform::Formulation, sp_solution_ids::Vector{VarId}
) 
    sp_uid = getuid(spform)
    nb_of_gen_col = 0

    for sol_id in sp_solution_ids
        nb_of_gen_col += 1
        name = string("MC_", getsortuid(sol_id)) 

        ## TODO : provide methods or operator to do this in one line
        master_coef_matrix = getcoefmatrix(masterform)
        sp_sol = getprimalsolmatrix(spform)[:,sol_id]
        members = ConstrMembership()
    
        for (sp_var_id, sp_var_val) in sp_sol
            for (master_constrid, sp_var_coef) in master_coef_matrix[:,sp_var_id]
                val = get(members, master_constrid, 0.0)
                members[master_constrid] = val + sp_var_val * sp_var_coef
            end
        end
        ## end TODO
        
        mc = setvar!(
            masterform, name, MasterCol,
            cost = getprimalsolcosts(spform)[sol_id],
            members = members,
            id = sol_id
        )
        @logmsg LogLevel(-2) string("Generated column : ", name)
    end

    return nb_of_gen_col
end

contrib_improves_mlp(::Type{MinSense}, sp_primal_bound::Float64) = (sp_primal_bound < 0.0 - 1e-8)
contrib_improves_mlp(::Type{MaxSense}, sp_primal_bound::Float64) = (sp_primal_bound > 0.0 + 1e-8)
contrib_improves_mlp(sp_primal_bound::PrimalBound{MinSense}) = (sp_primal_bound < 0.0 - 1e-8)
contrib_improves_mlp(sp_primal_bound::PrimalBound{MaxSense}) = (sp_primal_bound > 0.0 + 1e-8)

function compute_pricing_db_contrib(
    spform::Formulation, sp_sol_primal_bound::PrimalBound, sp_lb::Float64,
    sp_ub::Float64
)
    # Since convexity constraints are not automated and there is no stab
    # the pricing_dual_bound_contrib is just the reduced cost * multiplicty
    if contrib_improves_mlp(sp_sol_primal_bound)
        contrib = sp_sol_primal_bound * sp_ub
    else
        contrib = sp_sol_primal_bound * sp_lb
    end
    return contrib
end

function solve_sp_to_gencol!(
    masterform::Formulation, spform::Formulation, dual_sol::DualSolution,
    sp_lb::Float64, sp_ub::Float64
)::Tuple{Bool,Vector{VarId},Vector{VarId},Float64}
    
    recorded_solution_ids = Vector{VarId}()
    sp_solution_ids_to_activate = Vector{VarId}()

    sp_is_feasible = true

    #dual_bound_contrib = 0 # Not used
    #pseudo_dual_bound_contrib = 0 # Not used

    # TODO renable this. Needed at least for the diving
    # if can_not_generate_more_col(princing_prob)
    #     return flag_cannot_generate_more_col
    # end

    # Compute target
    update_pricing_target!(spform)

    # Reset var bounds, var cost, sp minCost
    if update_pricing_problem!(spform, dual_sol) # Never returns true
        #     This code is never executed because update_pricing_prob always returns false
        #     @logmsg LogLevel(-3) "pricing prob is infeasible"
        #     # In case one of the subproblem is infeasible, the master is infeasible
        #     compute_pricing_dual_bound_contrib(alg, pricing_prob)
        #     return flag_is_sp_infeasible
    end

    # if alg.colgen_stabilization != nothing && true #= TODO add conds =#
    #     # switch off the reduced cost estimation when stabilization is applied
    # end

    # Solve sub-problem and insert generated columns in master
    # @logmsg LogLevel(-3) "optimizing pricing prob"
    ipform = SolveIpForm(deactivate_artificial_vars = false, enforce_integrality = false, log_level = 2)
    TO.@timeit Coluna._to "Pricing subproblem" begin
        sp_output = run!(ipform, spform, SolveIpFormInput(ObjValues(spform)))
    end
    sp_result = getresult(sp_output)

    pricing_db_contrib = compute_pricing_db_contrib(
        spform, get_ip_primal_bound(sp_result), sp_lb, sp_ub
    )

    if !isfeasible(sp_result)
        sp_is_feasible = false 
        # @logmsg LogLevel(-3) "pricing prob is infeasible"
        return sp_is_feasible, recorded_solution_ids, PrimalBound(spform)
    end

    if nb_ip_primal_sols(sp_result) > 0
        for sol in get_ip_primal_sols(sp_result)
            if contrib_improves_mlp(getobjsense(spform), getvalue(sol)) # has negative reduced cost
                insertion_status, col_id = setprimalsol!(spform, sol)
                if insertion_status
                    push!(recorded_solution_ids, col_id)
                elseif !insertion_status && !iscuractive(masterform, col_id)
                    push!(sp_solution_ids_to_activate, col_id)
                else
                    msg = """
                    Column already exists as $(getname(masterform, col_id)) and is already active.
                    """
                    @warn string(msg)
                end
            end
        end
    end

    return sp_is_feasible, recorded_solution_ids, sp_solution_ids_to_activate, pricing_db_contrib
end

function solve_sps_to_gencols!(
    reform::Reformulation, dual_sol::DualSolution, 
    sp_lbs::Dict{FormId, Float64}, sp_ubs::Dict{FormId, Float64}
)
    masterform = getmaster(reform)
    nb_new_cols = 0
    dual_bound_contrib = DualBound(masterform, 0.0)
    masterform = getmaster(reform)
    sps = get_dw_pricing_sps(reform)
    recorded_sp_solution_ids = Dict{FormId, Vector{VarId}}()
    sp_solution_to_activate = Dict{FormId, Vector{VarId}}()
    sp_dual_bound_contribs = Dict{FormId, Float64}()

    ### BEGIN LOOP TO BE PARALLELIZED
    for (spuid, spform) in sps
        gen_status, new_sp_solution_ids, sp_solution_ids_to_activate, sp_dual_contrib = solve_sp_to_gencol!(
            masterform, spform, dual_sol, sp_lbs[spuid], sp_ubs[spuid]
        )
        if gen_status # else Sp is infeasible: contrib = Inf
            recorded_sp_solution_ids[spuid] = new_sp_solution_ids
            sp_solution_to_activate[spuid] = sp_solution_ids_to_activate
        end
        sp_dual_bound_contribs[spuid] = sp_dual_contrib #float(contrib)
    end
    ### END LOOP TO BE PARALLELIZED

    nb_new_cols = 0
    for (spuid, spform) in sps
        dual_bound_contrib += sp_dual_bound_contribs[spuid]
        nb_new_cols += insert_cols_in_master!(masterform, spform, recorded_sp_solution_ids[spuid])
        for colid in sp_solution_to_activate[spuid]
            activate!(masterform, colid)
            nb_new_cols += 1
        end
    end    
    return (nb_new_cols, dual_bound_contrib)
end

function compute_master_db_contrib(
    alg::ColGenRuntimeData, restricted_master_sol_value::PrimalBound{S}
) where {S}
    # TODO: will change with stabilization
    return DualBound{S}(restricted_master_sol_value)
end

function calculate_lagrangian_db(
    data::ColGenRuntimeData, restricted_master_sol_value::PrimalBound{S},
    pricing_sp_dual_bound_contrib::DualBound{S}
) where {S}
    lagran_bnd = DualBound{S}(0.0)
    lagran_bnd += compute_master_db_contrib(data, restricted_master_sol_value)
    lagran_bnd += pricing_sp_dual_bound_contrib
    return lagran_bnd
end

function generatecolumns!(
    data::ColGenRuntimeData, reform::Reformulation, master_val, 
    dual_sol, sp_lbs, sp_ubs
)
    nb_new_columns = 0
    while true # TODO Replace this condition when starting implement stabilization
        nb_new_col, sp_db_contrib =  solve_sps_to_gencols!(reform, dual_sol, sp_lbs, sp_ubs)
        nb_new_columns += nb_new_col
        lagran_bnd = calculate_lagrangian_db(data, master_val, sp_db_contrib)
        update_ip_dual_bound!(data.incumbents, lagran_bnd)
        update_lp_dual_bound!(data.incumbents, lagran_bnd)
        if nb_new_col < 0
            # subproblem infeasibility leads to master infeasibility
            return -1
        end
        break # TODO : rm
    end
    return nb_new_columns
end

ph_one_infeasible_db(db::DualBound{MinSense}) = getvalue(db) > (0.0 + 1e-5)
ph_one_infeasible_db(db::DualBound{MaxSense}) = getvalue(db) < (0.0 - 1e-5)

function cg_main_loop!(algo::ColumnGeneration, data::ColGenRuntimeData, reform::Reformulation)
    nb_cg_iterations = 0
    # Phase II loop: Iterate while can generate new columns and
    # termination by bound does not apply
    masterform = getmaster(reform)
    sp_lbs = Dict{FormId, Float64}()
    sp_ubs = Dict{FormId, Float64}()

    # collect multiplicity current bounds for each sp
    for (sp_uid, spform) in get_dw_pricing_sps(reform)
        sp_lbs[sp_uid] = getcurrhs(masterform, reform.dw_pricing_sp_lb[sp_uid])
        sp_ubs[sp_uid] = getcurrhs(masterform, reform.dw_pricing_sp_ub[sp_uid])
    end

    while true

        master_time = @elapsed begin
            master_output = run!(SolveLpForm(), masterform, SolveLpFormInput())
        end
        master_result = getresult(master_output)
        master_val = get_lp_primal_bound(master_result)
        dual_sols = get_lp_dual_sols(master_result)

        if data.phase != 1 && !isfeasible(master_result)
            @warn string("Solver returned that LP restricted master is infeasible or unbounded ",
            "(feasibility status = ", getfeasibilitystatus(master_result),") during phase != 1.")
            data.is_feasible = false
            return
        end

        if nb_lp_primal_sols(master_result) > 0
            set_lp_primal_sol!(data.incumbents, get_best_lp_primal_sol(master_result))
            set_lp_primal_bound!(data.incumbents, get_lp_primal_bound(master_result))
            set_lp_dual_sol!(data.incumbents, get_best_lp_dual_sol(master_result))
        else
            @error string("Solver returned that the LP restricted master is feasible but ",
            "did not return a primal solution. ",
            "Please open an issue (https://github.com/atoptima/Coluna.jl/issues).")
        end

        if nb_ip_primal_sols(master_result) > 0
            # if algo.store_all_ip_primal_sols
            add_ip_primal_sol!(data.incumbents, get_best_ip_primal_sol(master_result))
        end

        # TODO: cleanup restricted master columns        

        nb_cg_iterations += 1

        # generate new columns by solving the subproblems
        sp_time = @elapsed begin
            nb_new_col = generatecolumns!(
                data, reform, master_val, dual_sols[1], sp_lbs, sp_ubs
            )
        end

        if nb_new_col < 0
            @error "Infeasible subproblem."
            data.is_feasible = false        
            return 
        end

        print_intermediate_statistics(
            data, nb_new_col, nb_cg_iterations, master_time, sp_time
        )

        # TODO: update colgen stabilization

        dual_bound = get_ip_dual_bound(data.incumbents)
        primal_bound = get_lp_primal_bound(data.incumbents)     
        ip_primal_bound = get_ip_primal_bound(data.incumbents)

        if diff(dual_bound, ip_primal_bound) < algo.optimality_tol
            data.has_converged = true
            @logmsg LogLevel(1) "Dual bound reached primal bound."
            return 
        end
        if data.phase == 1 && ph_one_infeasible_db(dual_bound)
            data.is_feasible = false
            @logmsg LogLevel(1) "Phase one determines infeasibility."
            return 
        end
        if nb_new_col == 0 || gap(primal_bound, dual_bound) < algo.optimality_tol
            @logmsg LogLevel(1) "Column Generation Algorithm has converged."
            data.has_converged = true
            return 
        end
        if nb_cg_iterations > algo.max_nb_iterations
            @warn "Maximum number of column generation iteration is reached."
            return 
        end
    end
    return 
end

function print_intermediate_statistics(
    algdata::ColGenRuntimeData, nb_new_col::Int, nb_cg_iterations::Int,
    mst_time::Float64, sp_time::Float64
)
    mlp = getvalue(get_lp_primal_bound(algdata.incumbents))
    db = getvalue(get_ip_dual_bound(algdata.incumbents))
    pb = getvalue(get_ip_primal_bound(algdata.incumbents))
    @printf(
        "<it=%3i> <et=%5.2f> <mst=%5.2f> <sp=%5.2f> <cols=%2i> <mlp=%10.4f> <DB=%10.4f> <PB=%.4f>\n",
        nb_cg_iterations, Coluna._elapsed_solve_time(), mst_time, sp_time, nb_new_col, mlp, db, pb
    )
    return
end
