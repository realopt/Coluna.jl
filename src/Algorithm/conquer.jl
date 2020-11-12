
"""
    ConquerInput

    Input of a conquer algorithm used by the tree search algorithm.
    Contains the node in the search tree and the collection of storages to restore 
    before running the conquer algorithm. This collection of storages is passed
    in the input so that it is not obtained each time the conquer algorithm runs. 
"""
struct ConquerInput <: AbstractInput 
    node::Node    
    storages_to_restore::StoragesUsageDict
end

getnode(input::ConquerInput) = input.node

restore_states!(input::ConquerInput) = restore_states!(input.node.stateids, input.storages_to_restore) 

"""
    AbstractConquerAlgorithm

    This algorithm type is used by the tree search algorithm to update the incumbents and the formulation.
    For the moment, a conquer algorithm can be run only on reformulation.     
    A conquer algorithm should restore states of storages using function restore_states!(::ConquerInput)
        - each time it runs in the beginning
        - each time after calling a child manager algorithm
"""
abstract type AbstractConquerAlgorithm <: AbstractAlgorithm end

# conquer algorithms are always manager algorithms (they manage storing and restoring storages)
ismanager(algo::AbstractConquerAlgorithm) = true

function run!(algo::AbstractConquerAlgorithm, data::ReformData, input::ConquerInput)
    algotype = typeof(algo)
    error(string("Method run! which takes as parameters ReformData and ConquerInput ", 
                 "is not implemented for algorithm $algotype.")
    )
end    

# this function is needed in strong branching (to have a better screen logging)
isverbose(algo::AbstractConquerAlgorithm) = false

# this function is needed to check whether the best primal solution should be copied to the node optimization state
exploits_primal_solutions(algo::AbstractConquerAlgorithm) = false

# returns the optimization part of the output of the conquer algorithm 
function apply_conquer_alg_to_node!(
    node::Node, algo::AbstractConquerAlgorithm, data::ReformData, storages_to_restore::StoragesUsageDict
)  
    nodestate = getoptstate(node)
    if isverbose(algo)
        @logmsg LogLevel(-1) string("Node IP DB: ", get_ip_dual_bound(nodestate))
        @logmsg LogLevel(-1) string("Tree IP PB: ", get_ip_primal_bound(nodestate))
    end
    if ip_gap_closed(nodestate)
        @info "IP Gap is closed: $(ip_gap(getincumbents(nodestate))). Abort treatment."
    else
        isverbose(algo) && @logmsg LogLevel(-1) string("IP Gap is positive. Need to treat node.")

        run!(algo, data, ConquerInput(node, storages_to_restore))
        store_states!(data, node.stateids)
    end
    node.conquerwasrun = true
    return
end

####################################################################
#                      ParameterisedHeuristic
####################################################################

RestrictedMasterIPHeuristic() = SolveIpForm(get_dual_bound = false)

struct ParameterisedHeuristic
    algorithm::AbstractOptimizationAlgorithm
    root_priority::Float64
    nonroot_priority::Float64
    frequency::Integer
    max_depth::Integer
    name::String 
end

DefaultRestrictedMasterHeuristic() = 
    ParameterisedHeuristic(RestrictedMasterIPHeuristic(), 1.0, 1.0, 1, 1000, "Restricted Master IP")


####################################################################
#                      BendersConquer
####################################################################

@with_kw struct BendersConquer <: AbstractConquerAlgorithm 
    benders::BendersCutGeneration = BendersCutGeneration()
end

isverbose(strategy::BendersConquer) = true

# BendersConquer does not use any storage for the moment, it just calls 
# BendersCutSeparation algorithm, therefore get_storages_usage() is not defined for it

function get_child_algorithms(algo::BendersConquer, reform::Reformulation) 
    return [(algo.benders, reform)]
end

function run!(algo::BendersConquer, data::ReformData, input::ConquerInput)
    restore_states!(input)
    node = getnode(input)    
    nodestate = getoptstate(node)
    output = run!(algo.benders, data, OptimizationInput(nodestate))
    update!(nodestate, getoptstate(output))
    return 
end


####################################################################
#                      ColCutGenConquer
####################################################################

"""
    Coluna.Algorithm.ColCutGenConquer(
        colgen::AbstractOptimizationAlgorithm = ColumnGeneration()
        primal_heuristics::Vector{ParameterisedHeuristic} = [DefaultRestrictedMasterHeuristic()]
        preprocess = PreprocessAlgorithm()
        cutgen = CutCallbacks()
        run_preprocessing::Bool = false
    )

    Column-and-cut-generation based algorithm to find primal and dual bounds for a 
    problem decomposed using Dantzig-Wolfe paradigm. It applies `colgen` for the column 
    generation phase, `cutgen` for the cut generation phase, and it can apply several primal
    heuristics to more efficiently find feasible solutions.
"""
@with_kw struct ColCutGenConquer <: AbstractConquerAlgorithm 
    colgen::AbstractOptimizationAlgorithm = ColumnGeneration()
    primal_heuristics::Vector{ParameterisedHeuristic} = [DefaultRestrictedMasterHeuristic()]
    preprocess = PreprocessAlgorithm()
    cutgen = CutCallbacks()
    max_nb_cut_rounds::Int = 3 # TODO : tailing-off ?
    run_preprocessing::Bool = false
    opt_atol::Float64 = colgen.opt_atol # TODO : force this value in an init() method
    opt_rtol::Float64 = colgen.opt_rtol # TODO : force this value in an init() method
end

isverbose(algo::ColCutGenConquer) = algo.colgen.log_print_frequency > 0

# ColCutGenConquer does not use any storage for the moment, therefore 
# get_storages_usage() is not defined for it

function get_child_algorithms(algo::ColCutGenConquer, reform::Reformulation) 
    child_algos = Tuple{AbstractAlgorithm, AbstractModel}[]
    push!(child_algos, (algo.colgen, reform))
    push!(child_algos, (algo.cutgen, getmaster(reform)))
    algo.run_preprocessing && push!(child_algos, (algo.preprocess, reform))
    for heuristic in algo.primal_heuristics
        push!(child_algos, (heuristic.algorithm, reform))
    end
    return child_algos
end

function run!(algo::ColCutGenConquer, data::ReformData, input::ConquerInput)
    restore_states!(input)
    node = getnode(input)
    nodestate = getoptstate(node)
    reform = getreform(data)
    if algo.run_preprocessing && isinfeasible(run!(algo.preprocess, data, EmptyInput()))
        setterminationstatus!(nodestate, INFEASIBLE)
        return 
    end

    nb_tightening_rounds = 0
    colgen_output = run!(algo.colgen, data, OptimizationInput(nodestate))
    update!(nodestate, getoptstate(colgen_output))

    node_pruned = getterminationstatus(nodestate) == INFEASIBLE ||
        ip_gap_closed(nodestate, atol = algo.opt_atol, rtol = algo.opt_rtol)

    while !node_pruned && nb_tightening_rounds < algo.max_nb_cut_rounds
        sol = get_best_lp_primal_sol(getoptstate(colgen_output))
        if sol !== nothing
            cutcb_input = CutCallbacksInput(sol)
            cutcb_output = run!(CutCallbacks(), getmasterdata(data), cutcb_input)
            cutcb_output.nb_cuts_added == 0 && break
        else
            @warn "Skip cut generation because no best primal solution."
            break
        end

        set_ip_dual_bound!(nodestate, DualBound(reform))
        set_lp_dual_bound!(nodestate, DualBound(reform))
        colgen_output = run!(algo.colgen, data, OptimizationInput(nodestate))
        update!(nodestate, getoptstate(colgen_output))

        node_pruned = getterminationstatus(nodestate) == INFEASIBLE ||
            ip_gap_closed(nodestate, atol = algo.opt_atol, rtol = algo.opt_rtol)

        nb_tightening_rounds += 1
    end

    heuristics_to_run = Tuple{AbstractOptimizationAlgorithm, String, Float64}[]
    for heuristic in algo.primal_heuristics
        #TO DO : get_tree_order of nodes in strong branching is always -1
        if getdepth(node) <= heuristic.max_depth && 
            mod(get_tree_order(node) - 1, heuristic.frequency) == 0
            push!(heuristics_to_run, (
                heuristic.algorithm, heuristic.name,
                isrootnode(node) ? heuristic.root_priority : heuristic.nonroot_priority
            ))
        end
    end
    sort!(heuristics_to_run, by = x -> last(x), rev=true)

    for (heur_algorithm, name, priority) in heuristics_to_run
        node_pruned = ip_gap_closed(
            nodestate, atol = algo.opt_atol, rtol = algo.opt_rtol
            ) 
        node_pruned && break

        @info "Running $name heuristic"
        ismanager(heur_algorithm) && (stateids = store_states!(data))
        heur_output = run!(heur_algorithm, data, OptimizationInput(nodestate))
        update_all_ip_primal_solutions!(nodestate, getoptstate(heur_output))
        ismanager(heur_algorithm) && restore_states!(stateids, input.storages_to_restore)
    end 

    if node_pruned
        setterminationstatus!(nodestate, OPTIMAL)
    else
        setterminationstatus!(nodestate, OTHER_LIMIT)
    end
    return
end

####################################################################
#                      RestrMasterLPConquer
####################################################################

@with_kw struct RestrMasterLPConquer <: AbstractConquerAlgorithm 
    masterlpalgo::SolveLpForm = SolveLpForm(
        update_ip_primal_solution = true, consider_partial_solution = true
        )
end

# RestrMasterLPConquer does not use any storage, therefore get_storages_usage() is not defined for it

function get_child_algorithms(algo::RestrMasterLPConquer, reform::Reformulation) 
    return [(algo.masterlpalgo, getmaster(reform))]
end

function run!(algo::RestrMasterLPConquer, data::ReformData, input::ConquerInput)
    restore_states!(input)
    node = getnode(input)
    nodestate = getoptstate(node)
    output = run!(algo.masterlpalgo, getmasterdata(data), OptimizationInput(nodestate))
    masterlp_state =  getoptstate(output)
    update!(nodestate, masterlp_state)
    if ip_gap_closed(masterlp_state)
        setterminationstatus!(nodestate, OPTIMAL)
    else
        setterminationstatus!(nodestate, OTHER_LIMIT)
    end
    return
end

