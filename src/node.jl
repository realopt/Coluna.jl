## Defining infos here
@hl mutable struct SetupInfo end

mutable struct TreatOrder
    value::Int
end

@hl mutable struct Node
    params::Params
    children::Vector{Node}
    depth::Int
    # prune_dat_treat_node_start::Bool
    # estimated_sub_tree_size::Int
    # sub_tree_size::Int

    node_inc_lp_dual_bound::Float
    node_inc_ip_dual_bound::Float
    node_inc_lp_primal_bound::Float
    node_inc_ip_primal_bound::Float

    dual_bound_is_updated::Bool
    ip_primal_bound_is_updated::Bool

    node_inc_ip_primal_sol::PrimalSolution
    # partial_solution::PrimalSolution

    # eval_end_time::Int
    treat_order::TreatOrder

    infeasible::Bool
    evaluated::Bool
    treated::Bool

    ### New information recorded when the node was generated
    local_branching_constraints::Vector{MasterBranchConstr}

    ### Information recorded by father
    problem_setup_info::SetupInfo
    # eval_info::EvalInfo
    # children_generation_info::ChildrenGenerationInfo
    # branching_eval_info::BranchingEvaluationInfo #for branching history

    # problem_and_eval_alg_info_saved::Bool
    primal_sol::PrimalSolution # More information than only ::PrimalSolution
    # strong_branch_phase_number::Int
    # strong_branch_node_number::Int

end

function NodeBuilder(problem::ExtendedProblem, dual_bound::Float,
                     problem_setup_info::SetupInfo, primal_sol::PrimalSolution = PrimalSolution())

    return (
        problem.params,
        Node[],
        0,
        dual_bound,
        dual_bound,
        problem.primal_inc_bound,
        problem.primal_inc_bound,
        false,
        false,
        PrimalSolution(),
        TreatOrder(-1),
        false,
        false,
        false,
        MasterBranchConstr[],
        problem_setup_info,
        primal_sol,
    )
end

@hl mutable struct NodeWithParent <: Node
    parent::Node
end

function NodeWithParentBuilder(problem::ExtendedProblem, parent::Node)

    return tuplejoin(NodeBuilder(problem, parent.node_inc_ip_dual_bound,
        parent.problem_setup_info),
        parent)

end

@hl mutable struct DivingNode <: Node end

function DivingNodeBuilder(problem::ExtendedProblem, dual_bound::Float,
                           problem_setup_info::SetupInfo, primal_sol::PrimalSolution)
    return NodeBuilder(problem, dual_bound, problem_setup_info, primal_sol)
end

@hl mutable struct DivingNodeWithParent <: DivingNode
    local_partial_sol::PrimalSolution
    parent::DivingNode
end

function DivingNodeWithParentBuilder(problem::ExtendedProblem, local_partial_sol::PrimalSolution,
                                     parent::DivingNode)
    return tuplejoin(NodeBuilder(problem, parent.node_inc_ip_dual_bound, 
                                 parent.problem_setup_info), local_partial_sol, parent)
end

function DivingNodeWithParentBuilder(problem::ExtendedProblem, local_partial_sol_pair::Tuple{MasterColumn,Float},
                                     parent::DivingNode) 
    local_partial_sol = PrimalSolution()
    master_col = local_partial_sol_pair[1]
    local_partial_sol.cost = master_col.cur_cost_rhs
    local_partial_sol.var_val_map[master_col] = local_partial_sol_pair[2]
    return DivingNodeWithParentBuilder(problem, local_partial_sol, parent)
end

function get_priority(node::Node)
    if node.params.search_strategy == DepthFirst
        return node.depth
    elseif node.params.search_strategy == BestDualBound
        return node.node_inc_lp_dual_bound
    end
end

function is_conquered(node::Node)
    return (node.node_inc_ip_primal_bound - node.node_inc_ip_dual_bound
            <= node.params.mip_tolerance_integrality)
end

function is_to_be_pruned(node::Node, global_primal_bound::Float)
    return (global_primal_bound - node.node_inc_ip_dual_bound
        <= node.params.mip_tolerance_integrality)
end

function set_branch_and_price_order(node::Node, new_value::Int)
    node.treat_order = new_value
end

function exit_treatment(node::Node)
    # Issam: No need for deleting. I prefer deleting the node and storing the info
    # needed for printing the tree in a different light structure (for now)
    # later we can use Nullable for big data such as XXXInfo of node

    node.evaluated = true
    node.treated = true
end

function mark_infeasible_and_exit_treatment(node::Node)
    node.infeasible = true
    node.node_inc_lp_dual_bound = node.node_inc_ip_dual_bound = Inf
    exit_treatment(node)
end

function record_ip_primal_sol_and_update_ip_primal_bound(node::Node,
        sols_and_bounds)

    if node.node_inc_ip_primal_bound > sols_and_bounds.alg_inc_ip_primal_bound
        sol = PrimalSolution(sols_and_bounds.alg_inc_ip_primal_bound,
                             sols_and_bounds.alg_inc_ip_primal_sol_map)
        node.node_inc_ip_primal_sol = sol
        node.node_inc_ip_primal_bound = sols_and_bounds.alg_inc_ip_primal_bound
        node.ip_primal_bound_is_updated = true
    end
end

function update_node_duals(node::Node, sols_and_bounds)
    lp_dual_bound = sols_and_bounds.alg_inc_lp_dual_bound
    ip_dual_bound = sols_and_bounds.alg_inc_ip_dual_bound
    if node.node_inc_lp_dual_bound < lp_dual_bound
        node.node_inc_lp_dual_bound = lp_dual_bound
        node.dual_bound_is_updated = true
    end
    if node.node_inc_ip_dual_bound < ip_dual_bound
        node.node_inc_ip_dual_bound = ip_dual_bound
        node.dual_bound_is_updated = true
    end
end

function update_node_primals(node::Node, sols_and_bounds)
    # sols_and_bounds = node.alg_eval_node.sols_and_bounds
    if sols_and_bounds.is_alg_inc_ip_primal_bound_updated
        record_ip_primal_sol_and_update_ip_primal_bound(node,
            sols_and_bounds)
    end
    node.node_inc_lp_primal_bound = sols_and_bounds.alg_inc_lp_primal_bound
    node.primal_sol = PrimalSolution(node.node_inc_lp_primal_bound,
        sols_and_bounds.alg_inc_lp_primal_sol_map)
end

function update_node_primal_inc(node::Node, ip_bound::Float,
                                sol_map::Dict{Variable, Float})
    if ip_bound < node.node_inc_ip_primal_sol.cost
        new_sol = PrimalSolution(ip_bound, sol_map)
        node.node_inc_ip_primal_sol = new_sol
        node.node_inc_ip_primal_bound = ip_bound
        node.ip_primal_bound_is_updated = true
        if ip_bound < node.node_inc_lp_primal_bound
            node.node_inc_lp_primal_bound = ip_bound
            node.primal_sol = new_sol
        end
    end
end

function update_node_sols(node::Node, sols_and_bounds)
    update_node_primals(node, sols_and_bounds)
    update_node_duals(node, sols_and_bounds)
end

@hl mutable struct AlgLike end

function run(::AlgLike)
    @logmsg LogLevel(0) "Empty algorithm"
    return false
end

function to(alg::AlgLike; args...)
    return alg.extended_problem.timer_output
end

mutable struct TreatAlgs
    alg_setup_node::AlgLike
    alg_preprocess_node::AlgLike
    alg_eval_node::AlgLike
    alg_setdown_node::AlgLike
    alg_vect_primal_heur_node::Vector{AlgLike}
    alg_generate_children_nodes::AlgLike
    TreatAlgs() = new(AlgLike(), AlgLike(), AlgLike(), AlgLike(), AlgLike[], AlgLike())
end

function evaluation(node::Node, treat_algs::TreatAlgs,
                    global_treat_order::TreatOrder,
                    inc_primal_bound::Float)::Bool
    # node.treat_order = TreatOrder(global_treat_order.value)
    node.node_inc_ip_primal_bound = inc_primal_bound
    node.ip_primal_bound_is_updated = false
    node.dual_bound_is_updated = false

    run(treat_algs.alg_setup_node)

    # glpk_prob = treat_algs.alg_generate_children_nodes.extended_problem.master_problem.optimizer.optimizer.inner
    # GLPK.write_lp(glpk_prob, string("mip_ds", node.treat_order.value,".lp")) 

    if run(treat_algs.alg_preprocess_node)
        @logmsg LogLevel(0) string("Preprocess determines infeasibility.")
        run(treat_algs.alg_setdown_node)
        record_node_info(node, treat_algs.alg_setdown_node)
        mark_infeasible_and_exit_treatment(node)
        return true
    end

    # GLPK.write_lp(glpk_prob, string("mip_dp", node.treat_order.value,".lp")) 

    if run(treat_algs.alg_eval_node, inc_primal_bound)
        update_node_sols(node, treat_algs.alg_eval_node.sols_and_bounds)
        run(treat_algs.alg_setdown_node)
        record_node_info(node, treat_algs.alg_setdown_node)
        mark_infeasible_and_exit_treatment(node)
        return true
    end
    node.evaluated = true

    # GLPK.write_lp(glpk_prob, string("mip_de", node.treat_order.value,".lp")) 

    update_node_sols(node, treat_algs.alg_eval_node.sols_and_bounds)

    if is_conquered(node)
        @logmsg LogLevel(-2) string("Node is conquered, no need for branching.")
        run(treat_algs.alg_setdown_node)
        record_node_info(node, treat_algs.alg_setdown_node)
        exit_treatment(node);
        return true
    end

    run(treat_algs.alg_setdown_node)
    record_node_info(node, treat_algs.alg_setdown_node)

    # GLPK.write_lp(glpk_prob, string("mip_dd", node.treat_order.value,".lp")) 

    return true
end

function treat(node::Node, treat_algs::TreatAlgs,
        global_treat_order::TreatOrder, inc_primal_bound::Float)::Bool
    
    node.treat_order = TreatOrder(global_treat_order.value)
    global_treat_order.value += 1

    if !node.evaluated
        evaluation(node, treat_algs, global_treat_order, inc_primal_bound)
    end

    if node.treated
        @logmsg LogLevel(0) "Node is considered as treated after evaluation"
        return true
    end

    for alg in treat_algs.alg_vect_primal_heur_node
        run(alg, global_treat_order, node.primal_sol)
        update_node_primal_inc(node, alg.sols_and_bounds.alg_inc_ip_primal_bound,
                               alg.sols_and_bounds.alg_inc_ip_primal_sol_map)
        println("<", typeof(alg), ">", " <mlp=",
                node.node_inc_lp_primal_bound, "> ",
                "<PB=", node.node_inc_ip_primal_bound, ">")
        if is_conquered(node)
            @logmsg LogLevel(0) string("Node is considered conquered ",
                                       "after primal heuristic ", typeof(alg))
            exit_treatment(node)
            return true
        end
    end

    if !run(treat_algs.alg_generate_children_nodes, node.primal_sol)
        generate_children(node, treat_algs.alg_generate_children_nodes)
    end

    exit_treatment(node)

    return true
end
