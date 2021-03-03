"""
    Coluna.Algorithm.ColumnGeneration(
        restr_master_solve_alg = SolveLpForm(get_dual_solution = true)
        pricing_prob_solve_alg = SolveIpForm(
            deactivate_artificial_vars = false,
            enforce_integrality = false,
            log_level = 2
        ),
        essential_cut_gen_alg = CutCallbacks(call_robust_facultative = false)
        max_nb_iterations::Int = 1000
        log_print_frequency::Int = 1
        store_all_ip_primal_sols::Bool = false
        redcost_tol::Float = 1e-5
        cleanup_threshold::Int = 10000
        cleanup_ratio::Float = 0.66
        smoothing_stabilization::Float64 = 0.0 # should be in [0, 1]
        opt_atol::Float64 = DEF_OPTIMALITY_ATOL
        opt_rtol::Float64 = DEF_OPTIMALITY_RTOL
    )

Column generation algorithm. It applies `restr_master_solve_alg` to solve the linear
restricted master and `pricing_prob_solve_alg` to solve the subproblems.
"""
@with_kw struct ColumnGeneration <: AbstractOptimizationAlgorithm
    restr_master_solve_alg = SolveLpForm(get_dual_solution = true)
    #TODO : pricing problem solver may be different depending on the
    #       pricing subproblem
    pricing_prob_solve_alg = SolveIpForm(
        deactivate_artificial_vars = false,
        enforce_integrality = false,
        log_level = 2
    )
    essential_cut_gen_alg = CutCallbacks(call_robust_facultative = false)
    max_nb_iterations::Int64 = 1000
    log_print_frequency::Int64 = 1
    store_all_ip_primal_sols::Bool = false
    redcost_tol::Float64 = 1e-5
    solve_subproblems_parallel::Bool = false
    cleanup_threshold::Int64 = 10000
    cleanup_ratio::Float64 = 0.66
    smoothing_stabilization::Float64 = 0.0 # should be in [0, 1]
    opt_atol::Float64 = Coluna.DEF_OPTIMALITY_ATOL
    opt_rtol::Float64 = Coluna.DEF_OPTIMALITY_RTOL
end

stabilization_is_used(algo::ColumnGeneration) = !iszero(algo.smoothing_stabilization)

function get_child_algorithms(algo::ColumnGeneration, reform::Reformulation) 
    child_algs = Tuple{AbstractAlgorithm, AbstractModel}[]
    push!(child_algs, (algo.restr_master_solve_alg, getmaster(reform)))
    push!(child_algs, (algo.essential_cut_gen_alg, getmaster(reform)))
    for (id, spform) in get_dw_pricing_sps(reform)
        push!(child_algs, (algo.pricing_prob_solve_alg, spform))
    end
    return child_algs
end 

function get_storages_usage(algo::ColumnGeneration, reform::Reformulation) 
    storages_usage = Tuple{AbstractModel, StorageTypePair, StorageAccessMode}[] 
    master = getmaster(reform)
    push!(storages_usage, (master, MasterColumnsStoragePair, READ_AND_WRITE))
    push!(storages_usage, (master, PartialSolutionStoragePair, READ_ONLY))
    if stabilization_is_used(algo)
        push!(storages_usage, (master, ColGenStabilizationStoragePair, READ_AND_WRITE))
    end
    return storages_usage
end

struct ReducedCostsVector
    length::Int
    varids::Vector{VarId}
    perencosts::Vector{Float64}
    form::Vector{Formulation}
end

function ReducedCostsVector(varids::Vector{VarId}, form::Vector{Formulation})
    len = length(varids)
    perencosts = zeros(Float64, len)
    p = sortperm(varids)
    permute!(varids, p)
    permute!(form, p)
    for i in 1:len
        perencosts[i] = getcurcost(getmaster(form[i]), varids[i])
    end
    return ReducedCostsVector(len, varids, perencosts, form)
end

function run!(algo::ColumnGeneration, env::Env, data::ReformData, input::OptimizationInput)::OptimizationOutput
    reform = getreform(data)
    master = getmaster(reform)
    optstate = OptimizationState(master, getoptstate(input), false, false)

    restart = true
    stop = false
    while restart && !stop
        set_ph3!(master) # mixed ph1 & ph2
        stop, restart = cg_main_loop!(algo, env, 3, optstate, data)
    end

    restart = true
    while should_do_ph_1(optstate) && restart && !stop
        set_ph1!(master, optstate)
        stop, restart = cg_main_loop!(algo, env, 1, optstate, data)
        restart && break
        if !stop
            set_ph2!(master, optstate) # pure ph2
            stop, restart = cg_main_loop!(algo, env, 2, optstate, data)
        end
    end

    @logmsg LogLevel(-1) string("ColumnGeneration terminated with status ", getterminationstatus(optstate))

    return OptimizationOutput(optstate)
end

function should_do_ph_1(optstate::OptimizationState)
    primal_lp_sol = get_lp_primal_sols(optstate)[1]
    if contains(primal_lp_sol, vid -> isanArtificialDuty(getduty(vid)))
        @logmsg LogLevel(-2) "Artificial variables in lp solution, need to do phase one"
        return true
    end

    @logmsg LogLevel(-2) "No artificial variables in lp solution, will not proceed to do phase one"
    return false
end

function set_ph1!(master::Formulation, optstate::OptimizationState)
    for (varid, var) in getvars(master)
        if !isanArtificialDuty(getduty(varid))
            setcurcost!(master, varid, 0.0)
        end
    end
    set_lp_dual_bound!(optstate, DualBound(master))
    set_ip_dual_bound!(optstate, DualBound(master))
    return
end

function set_ph2!(master::Formulation, optstate::OptimizationState)
    for (varid, var) in getvars(master)
        if isanArtificialDuty(getduty(varid))
            deactivate!(master, varid)
        else
            setcurcost!(master, varid, getperencost(master, var))
        end
    end
    set_lp_dual_bound!(optstate, DualBound(master))
    set_ip_dual_bound!(optstate, DualBound(master))
    set_lp_primal_bound!(optstate, PrimalBound(master))
    return
end

function set_ph3!(master::Formulation)
    for (varid, var) in getvars(master)
        if isanArtificialDuty(getduty(varid))
            activate!(master, varid)
        else
            setcurcost!(master, varid, getperencost(master, var))
        end
    end
    return
end

function update_pricing_target!(spform::Formulation)
    # println("pricing target will only be needed after automating convexity constraints")
    return
end

mutable struct SubprobInfo
    lb_constr_id::ConstrId
    ub_constr_id::ConstrId
    lb::Float64
    ub::Float64
    lb_dual::Float64
    ub_dual::Float64
    bestsol::Union{Nothing, PrimalSolution}
    valid_dual_bound_contrib::Float64
    pseudo_dual_bound_contrib::Float64
    recorded_sol_ids::Vector{VarId}
    sol_ids_to_activate::Vector{VarId}
    isfeasible::Bool
end

function SubprobInfo(reform::Reformulation, spformid::FormId)
    master = getmaster(reform)
    lb_constr_id = get_dw_pricing_sp_lb_constrid(reform, spformid)
    ub_constr_id = get_dw_pricing_sp_ub_constrid(reform, spformid)
    lb = getcurrhs(master, lb_constr_id)
    ub = getcurrhs(master, ub_constr_id)
    return SubprobInfo(
        lb_constr_id, ub_constr_id, lb, ub, 0.0, 0.0, nothing, 0.0, 0.0,
        Vector{VarId}(), Vector{VarId}(), true
    )
end

function clear_before_colgen_iteration!(spinfo::SubprobInfo)
    spinfo.lb_dual = 0.0
    spinfo.ub_dual = 0.0
    spinfo.bestsol = nothing
    spinfo.valid_dual_bound_contrib = 0.0
    spinfo.pseudo_dual_bound_contrib = 0.0
    spinfo.isfeasible = true
    empty!(spinfo.recorded_sol_ids)
    empty!(spinfo.sol_ids_to_activate)
    return
end

set_bestcol_id!(spinfo::SubprobInfo, varid::VarId) = spinfo.bestcol_id = varid

function insert_cols_in_master!(
    masterform::Formulation, spinfo::SubprobInfo, phase::Int64, spform::Formulation,
)
    sp_uid = getuid(spform)
    nb_of_gen_col = 0

    for sol_id in spinfo.recorded_sol_ids
        nb_of_gen_col += 1
        name = string("MC_", getsortuid(sol_id))
        lb = 0.0
        ub = Inf
        kind = Continuous
        duty = MasterCol
        mc = setcol_from_sp_primalsol!(
            masterform, spform, sol_id, name, duty; lb = lb, ub = ub, kind = kind
        )
        if phase == 1
            setcurcost!(masterform, mc, 0.0)
        end
        @logmsg LogLevel(-2) string("Generated column : ", name)
    end

    return nb_of_gen_col
end

function compute_db_contributions!(
    spinfo::SubprobInfo, dualbound::DualBound{MaxSense}, primalbound::PrimalBound{MaxSense}
)
    value = getvalue(dualbound)
    spinfo.valid_dual_bound_contrib = value <= 0 ? value * spinfo.lb : value * spinfo.ub
    value = getvalue(primalbound)
    spinfo.pseudo_dual_bound_contrib = value <= 0 ? value * spinfo.lb : value * spinfo.ub
    return
end

function compute_db_contributions!(
    spinfo::SubprobInfo, dualbound::DualBound{MinSense}, primalbound::PrimalBound{MinSense}
)
    value = getvalue(dualbound)
    spinfo.valid_dual_bound_contrib = value >= 0 ? value * spinfo.lb : value * spinfo.ub
    value = getvalue(primalbound)
    spinfo.pseudo_dual_bound_contrib = value >= 0 ? value * spinfo.lb : value * spinfo.ub
    return
end

function compute_red_cost(
    algo::ColumnGeneration, master::Formulation, spinfo::SubprobInfo,
    spsol::PrimalSolution, lp_dual_sol::DualSolution
)
    red_cost::Float64 = 0.0
    if stabilization_is_used(algo)
        master_coef_matrix = getcoefmatrix(master)
        for (varid, value) in spsol
            red_cost += getcurcost(master, varid) * value
            for (constrid, var_coeff) in @view master_coef_matrix[:,varid]
                red_cost -= value * var_coeff * lp_dual_sol[constrid]
            end
        end
    else
        red_cost = getvalue(spsol)
    end
    #red_cost -= (spinfo.lb * spinfo.lb_dual + spinfo.ub * spinfo.ub_dual)
    red_cost -= (spinfo.lb_dual + spinfo.ub_dual)
    return red_cost
end

function improving_red_cost(redcost::Float64, algo::ColumnGeneration, ::Type{MinSense})
    return (redcost < 0.0 - algo.redcost_tol)
end

function improving_red_cost(redcost::Float64, algo::ColumnGeneration, ::Type{MaxSense})
    return (redcost > 0.0 + algo.redcost_tol)
end

function solve_sp_to_gencol!(
    spinfo::SubprobInfo, algo::ColumnGeneration, env::Env, masterform::Formulation, 
    spdata::ModelData, dualsol::DualSolution
)
    spform = getmodel(spdata)

    # Compute target
    update_pricing_target!(spform)

    output = run!(algo.pricing_prob_solve_alg, env, spdata, OptimizationInput(OptimizationState(spform)))
    sp_optstate = getoptstate(output)
    sp_sol_value = get_ip_primal_bound(sp_optstate)

    compute_db_contributions!(spinfo, get_ip_dual_bound(sp_optstate), sp_sol_value)

    sense = getobjsense(masterform)
    if nb_ip_primal_sols(sp_optstate) > 0
        bestsol = get_best_ip_primal_sol(sp_optstate)
        #@show bestsol compute_red_cost(algo, masterform, spinfo, bestsol, dualsol)
        if bestsol.status == FEASIBLE_SOL
            spinfo.bestsol = bestsol
            spinfo.isfeasible = true
            for sol in get_ip_primal_sols(sp_optstate)
                if improving_red_cost(compute_red_cost(algo, masterform, spinfo, sol, dualsol), algo, sense)
                    insertion_status, col_id = setprimalsol!(spform, sol)
                    if insertion_status
                        push!(spinfo.recorded_sol_ids, col_id)
                    elseif !insertion_status && !iscuractive(masterform, col_id)
                        push!(spinfo.sol_ids_to_activate, col_id)
                    else
                        msg = """
                        Warning: Column already exists as $(getname(masterform, col_id)) and is already active.
                        """
                        println(string(msg))
                    end
                end
            end
        end
    end
    return
end

function updatereducedcosts!(reform::Reformulation, redcostsvec::ReducedCostsVector, dualsol::DualSolution)
    redcosts = deepcopy(redcostsvec.perencosts)
    master = getmaster(reform)
    # sign = getobjsense(master) == MinSense ? -1 : 1
    matrix = getcoefmatrix(master)

    crm = matrix.rowmajor

    constr_key_pos::Int = 1
    next_constr_key_pos::Int = 2

    row_start = 0
    row_end = 0
    row_pos = 0

    terms = Dict{VarId, Float64}(id => 0.0 for id in redcostsvec.varids)

    for dual_pos in 1:length(dualsol.sol.array)
        entry = dualsol.sol.array[dual_pos]
        if entry !== nothing
            constrid, val = entry
            while constr_key_pos <= length(crm.col_keys) && crm.col_keys[constr_key_pos] != constrid
                constr_key_pos += 1
            end
            (constr_key_pos > length(crm.col_keys)) && break
            next_constr_key_pos = constr_key_pos + 1
            while next_constr_key_pos <= length(crm.col_keys) && crm.col_keys[next_constr_key_pos] === nothing
                next_constr_key_pos += 1
            end

            row_start = crm.pcsc.semaphores[constr_key_pos] + 1
            row_end = length(crm.pcsc.pma.array)
            if next_constr_key_pos <= length(crm.col_keys)
                row_end = crm.pcsc.semaphores[next_constr_key_pos] - 1
            end

            for row_pos in row_start:row_end
                entry = crm.pcsc.pma.array[row_pos]
                if entry !== nothing
                    row_varid, coeff = entry
                    if getduty(row_varid) <= AbstractMasterRepDwSpVar
                        terms[row_varid] = get(terms, row_varid, 0.0) + val * coeff
                    end
                end
            end
            constr_key_pos = next_constr_key_pos
        end
    end

    for (i, varid) in enumerate(redcostsvec.varids)
        setcurcost!(redcostsvec.form[i], varid, redcosts[i] - terms[varid])
    end
    return redcosts
end

function solve_sps_to_gencols!(
    spinfos::Dict{FormId, SubprobInfo}, algo::ColumnGeneration, env::Env, phase::Int64, 
    data::ReformData, redcostsvec::ReducedCostsVector, lp_dual_sol::DualSolution, 
    smooth_dual_sol::DualSolution,
)
    reform = getreform(data)
    masterform = getmaster(reform)
    nb_new_cols = 0
    spsdatas = get_dw_pricing_datas(data)

    # update reduced costs
    TO.@timeit Coluna._to "Update reduced costs" begin
        updatereducedcosts!(reform, redcostsvec, smooth_dual_sol)
    end

    ### BEGIN LOOP TO BE PARALLELIZED
    if algo.solve_subproblems_parallel
        spuids = collect(keys(spsdatas))
        Threads.@threads for key in 1:length(spuids)
            spuid = spuids[key]
            spdata = spsdatas[spuid]
            solve_sp_to_gencol!(spinfos[spuid], algo, env, masterform, spdata, lp_dual_sol)
        end
    else
        for (spuid, spdata) in spsdatas
            solve_sp_to_gencol!(spinfos[spuid], algo, env, masterform, spdata, lp_dual_sol)
        end
    end
    ### END LOOP TO BE PARALLELIZED

    for (spuid, spinfo) in spinfos
        !spinfo.isfeasible && return -1
    end

    TO.@timeit Coluna._to "Inserting columns" begin
        nb_new_cols = 0
        for (spuid, spdata) in spsdatas
            spinfo = spinfos[spuid]
            nb_of_gen_cols = insert_cols_in_master!(masterform, spinfo, phase, getmodel(spdata))
            nb_new_cols += nb_of_gen_cols
            for colid in spinfo.sol_ids_to_activate
                activate!(masterform, colid)
                nb_new_cols += 1
            end
            if algo.smoothing_stabilization == 1 && !iszero(spinfo.ub) && spinfo.bestsol === nothing
                @error string("Solutions to all pricing subproblems should be available in order ",
                              " to used automatic dual price smoothing")
            end
        end
    end

    return nb_new_cols
end

can_be_in_basis(algo::ColumnGeneration, ::Type{MinSense}, redcost::Float64) =
    redcost < 0 + algo.redcost_tol

can_be_in_basis(algo::ColumnGeneration, ::Type{MaxSense}, redcost::Float64) =
    redcost > 0 - algo.redcost_tol

function cleanup_columns(algo::ColumnGeneration, iteration::Int64, data::ReformData)

    # we do columns clean up only on every 10th iteration in order not to spend
    # the time retrieving the reduced costs
    # TO DO : master cleanup should be done on every iteration, for this we need
    # to quickly check the number of active master columns
    iteration % 10 != 0 && return

    cols_with_redcost = Vector{Pair{Variable, Float64}}()
    master = getmodel(getmasterdata(data))
    for (id, var) in getvars(master)
        if getduty(id) <= MasterCol && iscuractive(master, var) && isexplicit(master, var)
            push!(cols_with_redcost, var => getreducedcost(master, var))
        end
    end

    num_active_cols = length(cols_with_redcost)
    num_active_cols < algo.cleanup_threshold && return

    # sort active master columns by reduced cost
    reverse_order = getobjsense(master) == MinSense ? true : false
    sort!(cols_with_redcost, by = x -> x.second, rev=reverse_order)

    num_cols_to_keep = floor(Int64, num_active_cols * algo.cleanup_ratio)

    resize!(cols_with_redcost, num_active_cols - num_cols_to_keep)

    num_cols_removed::Int64 = 0
    for (var, redcost) in cols_with_redcost
        # we can remove column only if we are sure is it not in the basis
        # TO DO : we need to get the basis from the LP solver to have this verification
        if !can_be_in_basis(algo, getobjsense(master), redcost)
            deactivate!(master, var)
            num_cols_removed += 1
        end
    end
    @logmsg LogLevel(-1) "Cleaned up $num_cols_removed master columns"
    return
end

ph_one_infeasible_db(algo, db::DualBound{MinSense}) = getvalue(db) > algo.opt_atol
ph_one_infeasible_db(algo, db::DualBound{MaxSense}) = getvalue(db) < - algo.opt_atol

function update_lagrangian_dual_bound!(
    stabstorage::ColGenStabilizationStorage, optstate::OptimizationState{F, S}, algo::ColumnGeneration,
    master::Formulation, puremastervars::Vector{Pair{VarId,Float64}}, dualsol::DualSolution,
    partialsol::PrimalSolution, spinfos::Dict{FormId, SubprobInfo}
) where {F, S}

    sense = getobjsense(master)

    puremastvars_contrib::Float64 = getvalue(partialsol)
    # if smoothing is not active the pure master variables contribution
    # is already included in the value of the dual solution
    if smoothing_is_active(stabstorage)
        master_coef_matrix = getcoefmatrix(master)
        for (varid, mult) in puremastervars
            redcost = getcurcost(master, varid)
            for (constrid, var_coeff) in @view master_coef_matrix[:,varid]
                redcost -= var_coeff * dualsol[constrid]
            end
            mult = improving_red_cost(redcost, algo, sense) ?
                getcurub(master, varid) : getcurlb(master, varid)
            puremastvars_contrib += redcost * mult
        end
    end

    valid_lagr_bound = DualBound{S}(puremastvars_contrib + dualsol.bound)
    for (spuid, spinfo) in spinfos
        valid_lagr_bound += spinfo.valid_dual_bound_contrib
    end

    update_ip_dual_bound!(optstate, valid_lagr_bound)
    update_lp_dual_bound!(optstate, valid_lagr_bound)

    if stabilization_is_used(algo)
        pseudo_lagr_bound = DualBound{S}(puremastvars_contrib + dualsol.bound)
        for (spuid, spinfo) in spinfos
            pseudo_lagr_bound += spinfo.pseudo_dual_bound_contrib
        end
        update_stability_center!(stabstorage, dualsol, valid_lagr_bound, pseudo_lagr_bound)
    end
    return
end

function compute_subgradient_contribution(
    algo::ColumnGeneration, stabstorage::ColGenStabilizationStorage, master::Formulation,
    puremastervars::Vector{Pair{VarId,Float64}}, spinfos::Dict{FormId, SubprobInfo}
)
    sense = getobjsense(master)
    constrids = ConstrId[]
    constrvals = Float64[]

    if subgradient_is_needed(stabstorage, algo.smoothing_stabilization)
        master_coef_matrix = getcoefmatrix(master)

        for (varid, mult) in puremastervars
            for (constrid, var_coeff) in @view master_coef_matrix[:,varid]
                push!(constrids, constrid)
                push!(constrvals, var_coeff * mult)
            end
        end

        for (spuid, spinfo) in spinfos
            iszero(spinfo.ub) && continue
            mult = improving_red_cost(spinfo.bestsol.bound, algo, sense) ? spinfo.ub : spinfo.lb
            for (sp_var_id, sp_var_val) in spinfo.bestsol
                for (master_constrid, sp_var_coef) in @view master_coef_matrix[:,sp_var_id]
                    if !(getduty(master_constrid) <= MasterConvexityConstr)
                        push!(constrids, master_constrid)
                        push!(constrvals, sp_var_coef * sp_var_val * mult)
                    end
                end
            end
        end
    end

    return DualSolution(master, constrids, constrvals, 0.0, UNKNOWN_SOLUTION_STATUS)
end

function move_convexity_constrs_dual_values!(
    spinfos::Dict{FormId, SubprobInfo}, dualsol::DualSolution
)
    newbound = dualsol.bound
    for (spuid, spinfo) in spinfos
        spinfo.lb_dual = dualsol[spinfo.lb_constr_id]
        spinfo.ub_dual = dualsol[spinfo.ub_constr_id]
        dualsol[spinfo.lb_constr_id] = zero(0.0)
        dualsol[spinfo.ub_constr_id] = zero(0.0)
        newbound -= (spinfo.lb_dual * spinfo.lb + spinfo.ub_dual * spinfo.ub)
    end
    constrids = Vector{ConstrId}()
    values = Vector{Float64}()
    for (constrid, value) in dualsol
        if !(getduty(constrid) <= MasterConvexityConstr)
            push!(constrids, constrid)
            push!(values, value)
        end
    end
    return DualSolution(dualsol.model, constrids, values, newbound, FEASIBLE_SOL)
end

function get_pure_master_vars(master::Formulation)
    puremastervars = Vector{Pair{VarId,Float64}}()
    for (varid, var) in getvars(master)
        if isanOriginalRepresentatives(getduty(varid)) &&
            iscuractive(master, var) && isexplicit(master, var)
            push!(puremastervars, varid => 0.0)
        end
    end
    return puremastervars
end

function change_values_sign!(dualsol::DualSolution)
    # note that the bound value remains the same
    for (constrid, value) in dualsol
        dualsol[constrid] = -value
    end
    return
end

# cg_main_loop! returns Tuple{Bool, Bool} :
# - first one is equal to true when colgen algorithm must stop.
# - second one is equal to true when colgen algorithm must restart. 
#   If col gen was at phase 3, it will restart at phase 3. 
#   If it was as phase 1 or 2, it will restart at phase 1.
function cg_main_loop!(
    algo::ColumnGeneration, env::Env, phase::Int, cg_optstate::OptimizationState, 
    data::ReformData
)
    # Phase II loop: Iterate while can generate new columns and
    # termination by bound does not apply
    reform = getreform(data)
    masterform = getmaster(reform)
    spinfos = Dict{FormId, SubprobInfo}()

    # collect multiplicity current bounds for each sp
    dwspvars = Vector{VarId}()
    dwspforms = Vector{Formulation}()
    pure_master_vars = get_pure_master_vars(masterform)

    for (spid, spform) in get_dw_pricing_sps(reform)
        spinfos[spid] = SubprobInfo(reform, spid)

        for (varid, var) in getvars(spform)
            if iscuractive(spform, varid) && getduty(varid) <= AbstractDwSpVar
                push!(dwspvars, varid)
                push!(dwspforms, spform)
            end
        end
    end

    redcostsvec = ReducedCostsVector(dwspvars, dwspforms)
    iteration = 0

    stabstorage = (stabilization_is_used(algo) ? getstorage(getmasterdata(data), ColGenStabilizationStoragePair) 
                                               : ColGenStabilizationStorage(masterform) )

    partsolstorage = getstorage(getmasterdata(data), PartialSolutionStoragePair)
    partial_solution = get_primal_solution(partsolstorage, masterform)

    init_stab_before_colgen_loop!(stabstorage)

    while true
        for (spuid, spinfo) in spinfos
            clear_before_colgen_iteration!(spinfo)
        end

        rm_time = @elapsed begin
            rm_input = OptimizationInput(
                OptimizationState(masterform, ip_primal_bound = get_ip_primal_bound(cg_optstate))
            )
            rm_output = run!(algo.restr_master_solve_alg, env, getmasterdata(data), rm_input)
        end
        rm_optstate = getoptstate(rm_output)

        if phase != 1 && getterminationstatus(rm_optstate) == INFEASIBLE
            @warn string("Solver returned that LP restricted master is infeasible or unbounded ",
            "(termination status = INFEASIBLE) during phase != 1.")
            setterminationstatus!(cg_optstate, INFEASIBLE)
            return true, false
        end

        lp_dual_sol = nothing
        if nb_lp_dual_sols(rm_optstate) > 0
            lp_dual_sol = get_best_lp_dual_sol(rm_optstate)
        else
            error("Solver returned that the LP restricted master is feasible but ",
            "did not return a dual solution. ",
            "Please open an issue (https://github.com/atoptima/Coluna.jl/issues).")
        end
        if getobjsense(masterform) == MaxSense
            # this is needed due to convention that MOI uses for signs of duals in the maximization case
            change_values_sign!(lp_dual_sol)
        end
        lp_dual_sol = move_convexity_constrs_dual_values!(spinfos, lp_dual_sol)
        
        TO.@timeit Coluna._to "Getting primal solution" begin
        if nb_lp_primal_sols(rm_optstate) > 0
            rm_sol = get_best_lp_primal_sol(rm_optstate)
            
            set_lp_primal_sol!(cg_optstate, rm_sol)
            lp_bound = get_lp_primal_bound(rm_optstate) + getvalue(partial_solution)
            set_lp_primal_bound!(cg_optstate, lp_bound)

            if phase != 1 && !contains(rm_sol, varid -> isanArtificialDuty(getduty(varid)))
                proj_sol = proj_cols_on_rep(rm_sol, masterform)
                if isinteger(proj_sol) && isbetter(lp_bound, get_ip_primal_bound(cg_optstate))
                    # Essential cut generation mandatory when colgen finds a feasible solution
                    new_primal_sol = cat(rm_sol, partial_solution)
                    cutcb_input = CutCallbacksInput(new_primal_sol)
                    cutcb_output = run!(
                        algo.essential_cut_gen_alg, env, getmasterdata(data), cutcb_input
                    )
                    if cutcb_output.nb_cuts_added == 0
                        update_ip_primal_sol!(cg_optstate, new_primal_sol)
                    else
                        # because the new cuts may make the master infeasible
                        return false, true
                    end
                end
            end
        else
            @error string("Solver returned that the LP restricted master is feasible but ",
            "did not return a primal solution. ",
            "Please open an issue (https://github.com/atoptima/Coluna.jl/issues).")
        end
        end # @timeit

        TO.@timeit Coluna._to "Cleanup columns" begin
            cleanup_columns(algo, iteration, data)
        end

        iteration += 1

        TO.@timeit Coluna._to "Smoothing update" begin
            smooth_dual_sol = update_stab_after_rm_solve!(stabstorage, algo.smoothing_stabilization, lp_dual_sol)
        end

        nb_new_columns = 0
        sp_time = 0
        while true
            sp_time += @elapsed begin
                nb_new_col = solve_sps_to_gencols!(
                    spinfos, algo, env, phase, data, redcostsvec, lp_dual_sol, 
                    smooth_dual_sol
                )
            end

            if nb_new_col < 0
                @error "Infeasible subproblem."
                setterminationstatus!(cg_optstate, INFEASIBLE)
                return true, false
            end

            nb_new_columns += nb_new_col

            TO.@timeit Coluna._to "Update Lagrangian bound" begin
                update_lagrangian_dual_bound!(
                    stabstorage, cg_optstate, algo, masterform, pure_master_vars, 
                    smooth_dual_sol, partial_solution, spinfos
                )
            end

            if stabilization_is_used(algo)
                TO.@timeit Coluna._to "Smoothing update" begin
                    smooth_dual_sol = update_stab_after_gencols!(
                        stabstorage, algo.smoothing_stabilization, nb_new_col, lp_dual_sol, smooth_dual_sol,
                        compute_subgradient_contribution(algo, stabstorage, masterform, pure_master_vars, spinfos)
                    )
                end
                smooth_dual_sol === nothing && break
            else
                break
            end
        end

        print_colgen_statistics(
            env, phase, iteration, stabstorage.curalpha, cg_optstate, nb_new_columns, 
            rm_time, sp_time
        )

        update_stab_after_colgen_iteration!(stabstorage)

        dual_bound = get_ip_dual_bound(cg_optstate)
        primal_bound = get_lp_primal_bound(cg_optstate)
        ip_primal_bound = get_ip_primal_bound(cg_optstate)

        if ip_gap_closed(cg_optstate, atol = algo.opt_atol, rtol = algo.opt_rtol)
            setterminationstatus!(cg_optstate, OPTIMAL)
            println("Dual bound reached primal bound.")
            return true, false
        end
        if phase == 1 && ph_one_infeasible_db(algo, dual_bound)
            db = - getvalue(DualBound(reform))
            pb = - getvalue(PrimalBound(reform))
            set_lp_dual_bound!(cg_optstate, DualBound(reform, db))
            set_lp_primal_bound!(cg_optstate, PrimalBound(reform, pb))
            setterminationstatus!(cg_optstate, INFEASIBLE)
            println("Phase one determines infeasibility.")
            return true, false
        end
        if lp_gap_closed(cg_optstate, atol = algo.opt_atol, rtol = algo.opt_rtol)
            println("Column generation algorithm has converged.")
            setterminationstatus!(cg_optstate, OPTIMAL)
            return false, false
        end
        if nb_new_columns == 0
            println("No new column generated by the pricing problems.")
            setterminationstatus!(cg_optstate, OTHER_LIMIT)
            return false, false
        end
        if iteration > algo.max_nb_iterations
            setterminationstatus!(cg_optstate, OTHER_LIMIT)
            println("Warning : Maximum number of column generation iteration is reached.")
            return true, false
        end
    end

    return false, false
end

function print_colgen_statistics(
    env::Env, phase::Int64, iteration::Int64, smoothalpha::Float64, 
    optstate::OptimizationState, nb_new_col::Int, mst_time::Float64, sp_time::Float64
)
    mlp = getvalue(get_lp_primal_bound(optstate))
    db = getvalue(get_lp_dual_bound(optstate))
    pb = getvalue(get_ip_primal_bound(optstate))
    phase_string = "  "
    if phase == 1
        phase_string = "# "
    elseif phase == 2
        phase_string = "##"
    end

    @printf(
        "%s<it=%3i> <et=%5.2f> <mst=%5.2f> <sp=%5.2f> <cols=%2i> <al=%5.2f> <DB=%10.4f> <mlp=%10.4f> <PB=%.4f>\n",
        phase_string, iteration, elapsed_optim_time(env), mst_time, sp_time, nb_new_col, smoothalpha, db, mlp, pb
    )
    return
end
