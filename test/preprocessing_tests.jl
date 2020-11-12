import Random

function preprocessing_tests()
    # @testset "play gap with preprocessing" begin
    #     play_gap_with_preprocessing_tests()
    # end
    @testset "Preprocessing with random instances" begin
        random_instances_tests()
    end
end

function gen_random_small_gap_instance()
    nb_jobs = Random.rand(7:15)
    nb_machs = Random.rand(2:3)
    data = CLD.GeneralizedAssignment.Data(nb_machs, nb_jobs)
    for m in 1:nb_machs
        data.capacity[m] = Random.rand(100:120)
    end
    avg_weight = sum(data.capacity)/nb_jobs
    for j in 1:nb_jobs, m in 1:nb_machs
        data.cost[j,m] = Random.rand(1:10)
    end
    for j in 1:nb_jobs, m in 1:nb_machs
        data.weight[j,m] = Int(ceil(0.1*Random.rand(6:12)*avg_weight))
    end
    return data
end

function play_gap_with_preprocessing_tests()
    data = CLD.GeneralizedAssignment.data("play2.txt")

    coluna = JuMP.optimizer_with_attributes(
        CL.Optimizer,
        "default_optimizer" => GLPK.Optimizer,
        "params" => CL.Params(
            solver = ClA.TreeSearchAlgorithm(
                conqueralg = ClA.ColCutGenConquer(run_preprocessing = true)
            )
        )
    )

    problem, x, dec = CLD.GeneralizedAssignment.model(data, coluna, true)
    JuMP.optimize!(problem)
    @test abs(JuMP.objective_value(problem) - 75.0) <= 0.00001
    @test JuMP.termination_status(problem) == MOI.OPTIMAL
    @test CLD.GeneralizedAssignment.print_and_check_sol(data, problem, x)
end

function random_instances_tests()
    Random.seed!(3)
    nb_prep_vars = 0
    nb_infeas = 0
    for problem_idx in 1:100
        res = test_random_gap_instance()
        if !res[1]
            nb_infeas += 1
        else
            nb_prep_vars += res[2]
        end
    end
    println("nb_infeas: $(nb_infeas) avg_prep_vars: $(nb_prep_vars/(100 - nb_infeas))")
    return
end

function apply_random_branching_constraint!(problem, x, m, j, leq)
    if leq
        @constraint(problem, random_br, x[m,j] <= 0)
    else
        @constraint(problem, random_br, x[m,j] >= 1)
    end
    return
end

function test_random_gap_instance()
    data = gen_random_small_gap_instance()
    coluna = JuMP.optimizer_with_attributes(
        CL.Optimizer,
        "default_optimizer" => GLPK.Optimizer,
        "params" => CL.Params(
            solver = ClA.TreeSearchAlgorithm(
                conqueralg = ClA.ColCutGenConquer(run_preprocessing = true),
                dividealg = ClA.NoBranching()
            )
        )
    )

    problem, x, dec = CLD.GeneralizedAssignment.model(data, coluna, true)
    # Adding a random branching constraint
    br_j = Random.rand(data.jobs)
    br_m = Random.rand(data.machines)
    leq = Random.rand(Bool)
    apply_random_branching_constraint!(problem, x, br_m, br_j, leq)
    JuMP.optimize!(problem)

    if JuMP.termination_status(problem) == MOI.INFEASIBLE
        coluna = JuMP.optimizer_with_attributes(
            CL.Optimizer,
            "default_optimizer" => GLPK.Optimizer,
            "params" => CL.Params(
                solver = ClA.TreeSearchAlgorithm(
                    conqueralg = ClA.ColCutGenConquer(run_preprocessing = false)
                )
            )
        )
        problem, x, dec = CLD.GeneralizedAssignment.model(data, coluna, true)
        apply_random_branching_constraint!(problem, x, br_m, br_j, leq)
        JuMP.optimize!(problem)
        @test JuMP.termination_status(problem) == MOI.INFEASIBLE
        return (false, 0)
    else
        nb_prep_vars = 0
        coluna_optimizer = problem.moi_backend
        master = CL.getmaster(coluna_optimizer.inner.re_formulation)
        for (moi_index, varid) in coluna_optimizer.varids
            var = CL.getvar(master, varid)
            if CL.getcurlb(master, var) == CL.getcurub(master, var)
                var_name = CL.getname(master, var)
                m = parse(Int, split(split(var_name, ",")[1], "[")[2])
                j = parse(Int, split(split(var_name, ",")[2], "]")[1])
                forbidden_machs = (
                    CL.getcurlb(master, var) == 1 ? [m] : [mach_idx for mach_idx in data.machines if mach_idx != m]
                )
                modified_data = deepcopy(data)
                for mach_idx in forbidden_machs
                    modified_data.weight[j,mach_idx] = modified_data.capacity[mach_idx] + 1
                end
                coluna = JuMP.optimizer_with_attributes(
                    CL.Optimizer,
                    "default_optimizer" => GLPK.Optimizer,
                    "params" => CL.Params(
                        solver = ClA.TreeSearchAlgorithm(
                            conqueralg = ClA.ColCutGenConquer(run_preprocessing = false)
                        )
                    )
                )
                modified_problem, x, dec = CLD.GeneralizedAssignment.model(modified_data, coluna, true)
                apply_random_branching_constraint!(modified_problem, x, br_m, br_j, leq)
                JuMP.optimize!(modified_problem)

                @test JuMP.termination_status(modified_problem) == MOI.INFEASIBLE
                nb_prep_vars += 1
            end
        end
        return (true, nb_prep_vars)
    end
    return
end
