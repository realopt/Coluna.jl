using Coluna.Algorithm
using Coluna.MathProg
using Parameters

function user_algorithms_tests()
    conseq_colgen_test()
end

@with_kw struct ConsecutiveColGen <: AbstractOptimizationAlgorithm
    colgen = ColumnGeneration(smoothing_stabilization = 1.0)
    preprocess = PreprocessAlgorithm(preprocess_subproblems = true, printing = true)
    rm_heur = RestrictedMasterIPHeuristic()
    num_calls_to_col_gen = 3  
end

function Coluna.Algorithm.get_child_algorithms(
    algo::ConsecutiveColGen, reform::Reformulation
    )
    return [(algo.colgen, reform), (algo.preprocess, reform), (algo.rm_heur, reform)]
end 

function Coluna.Algorithm.get_storages_usage(
    algo::ConsecutiveColGen, reform::Reformulation
    ) 
    master = Coluna.MathProg.getmaster(reform)
    return [(reform, Coluna.Algorithm.PreprocessingStoragePair, Coluna.Algorithm.READ_AND_WRITE),
            (master, Coluna.Algorithm.PartialSolutionStoragePair, Coluna.Algorithm.READ_AND_WRITE)]
end

function Coluna.Algorithm.run!(
    algo::ConsecutiveColGen, data::ReformData, input::OptimizationInput
    )
    reform = getreform(data)
    master = ClMP.getmaster(reform)
    masterdata = getmasterdata(data)
    optstate = getoptstate(input)

    cg_run_number = 1

    while cg_run_number <= algo.num_calls_to_col_gen 

        cg_output = run!(algo.colgen, data, OptimizationInput(optstate))
        cg_optstate = getoptstate(cg_output)
        update_all_ip_primal_solutions!(optstate, cg_optstate)

        primal_sol = get_best_lp_primal_sol(cg_optstate)

        if cg_run_number == 1
            set_ip_dual_bound!(optstate, get_ip_dual_bound(cg_optstate))
        end

        var_vals = [(var, val) for (var, val) in primal_sol] 
        isempty(var_vals) && break

        sort!(var_vals, by = x -> last(x), rev = true)

        preprocess_storage = getstorage(data, PreprocessingStoragePair)
        partsol_storage = getstorage(masterdata, PartialSolutionStoragePair)
    
        add_to_localpartialsol!(preprocess_storage, first(var_vals[1]), 1.0)
        add_to_solution!(partsol_storage, first(var_vals[1]), 1.0)

        prp_output = run!(algo.preprocess, data, EmptyInput())
        isinfeasible(prp_output) && break
    
        cg_run_number += 1
    end

    heur_output = run!(algo.rm_heur, data, OptimizationInput(optstate))
    update_all_ip_primal_solutions!(optstate, getoptstate(heur_output))

    return OptimizationOutput(optstate)
end

function conseq_colgen_test()
    
    data = CLD.GeneralizedAssignment.data("mediumgapcuts3.txt")

    coluna = JuMP.optimizer_with_attributes(
        CL.Optimizer,
        "params" => CL.Params(solver = ConsecutiveColGen(num_calls_to_col_gen = 3)),
        "default_optimizer" => GLPK.Optimizer
    )

    model, x, dec = CLD.GeneralizedAssignment.model(data, coluna)

    BD.objectiveprimalbound!(model, 2000.0)
    BD.objectivedualbound!(model, 0.0)

    JuMP.optimize!(model)    

    @test JuMP.objective_value(model) < 2000.0
end