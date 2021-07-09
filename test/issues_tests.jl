# Issue #425
# When the user does not provide decomposition, Coluna should optimize the
# original formulation.
function solve_with_no_decomposition()
    coluna = JuMP.optimizer_with_attributes(
        Coluna.Optimizer,
        "params" => CL.Params(solver=ClA.SolveIpForm()),
        "default_optimizer" => GLPK.Optimizer
    )

    model = BlockModel(coluna, direct_model=true)
    @variable(model, x)
    @constraint(model, x <= 1)
    @objective(model, Max, x)

    optimize!(model)
    @test JuMP.objective_value(model) == 1.0
    @test JuMP.termination_status(model) == MOI.OPTIMAL
end

# Test that empty! empties the Problem
function test_model_empty()
    coluna = JuMP.optimizer_with_attributes(
        Coluna.Optimizer,
        "params" => CL.Params(solver=ClA.SolveIpForm()),
        "default_optimizer" => GLPK.Optimizer
    )

    model = BlockModel(coluna, direct_model=true)
    @variable(model, x)
    @constraint(model, x <= 1)
    @objective(model, Max, x)

    optimize!(model)
    @test JuMP.objective_value(model) == 1.0
    @test JuMP.termination_status(model) == MOI.OPTIMAL

    empty!(model)
    @variable(model, x)
    @constraint(model, x <= 2)
    @objective(model, Max, x)

    optimize!(model)
    @test JuMP.objective_value(model) == 2.0
    @test JuMP.termination_status(model) == MOI.OPTIMAL
end

function decomposition_with_constant_in_objective()
    nb_machines = 4
    nb_jobs = 30
    c = [12.7 22.5 8.9 20.8 13.6 12.4 24.8 19.1 11.5 17.4 24.7 6.8 21.7 14.3 10.5 15.2 14.3 12.6 9.2 20.8 11.7 17.3 9.2 20.3 11.4 6.2 13.8 10.0 20.9 20.6;  19.1 24.8 24.4 23.6 16.1 20.6 15.0 9.5 7.9 11.3 22.6 8.0 21.5 14.7 23.2 19.7 19.5 7.2 6.4 23.2 8.1 13.6 24.6 15.6 22.3 8.8 19.1 18.4 22.9 8.0;  18.6 14.1 22.7 9.9 24.2 24.5 20.8 12.9 17.7 11.9 18.7 10.1 9.1 8.9 7.7 16.6 8.3 15.9 24.3 18.6 21.1 7.5 16.8 20.9 8.9 15.2 15.7 12.7 20.8 10.4;  13.1 16.2 16.8 16.7 9.0 16.9 17.9 12.1 17.5 22.0 19.9 14.6 18.2 19.6 24.2 12.9 11.3 7.5 6.5 11.3 7.8 13.8 20.7 16.8 23.6 19.1 16.8 19.3 12.5 11.0]
    w = [61 70 57 82 51 74 98 64 86 80 69 79 60 76 78 71 50 99 92 83 53 91 68 61 63 97 91 77 68 80; 50 57 61 83 81 79 63 99 82 59 83 91 59 99 91 75 66 100 69 60 87 98 78 62 90 89 67 87 65 100; 91 81 66 63 59 81 87 90 65 55 57 68 92 91 86 74 80 89 95 57 55 96 77 60 55 57 56 67 81 52;  62 79 73 60 75 66 68 99 69 60 56 100 67 68 54 66 50 56 70 56 72 62 85 70 100 57 96 69 65 50]
    Q = [1020 1460 1530 1190]

    coluna = optimizer_with_attributes(
        Coluna.Optimizer,
        "params" => Coluna.Params(
            solver=Coluna.Algorithm.TreeSearchAlgorithm() # default BCP
        ),
        "default_optimizer" => GLPK.Optimizer # GLPK for the master & the subproblems
    )

    @axis(M, 1:nb_machines)
    J = 1:nb_jobs

    model = BlockModel(coluna)
    @variable(model, x[m in M, j in J], Bin)
    @constraint(model, cov[j in J], sum(x[m, j] for m in M) >= 1)
    @constraint(model, knp[m in M], sum(w[m, j] * x[m, j] for j in J) <= Q[m])
    @objective(model, Min, sum(c[m, j] * x[m, j] for m in M, j in J) + 2)
    @dantzig_wolfe_decomposition(model, decomposition, M)
    optimize!(model)
    @test objective_value(model) ≈ 307.5 + 2
end

# Issue #424
# - If you try to solve an empty model with Coluna using a SolveIpForm or SolveLpForm
#   as top solver, the objective value will be 0.
# - If you try to solve an empty model using TreeSearchAlgorithm, then Coluna will
#   throw an error because since there is no decomposition, there is no reformulation
#   and TreeSearchAlgorithm must be run on a reformulation.
function solve_empty_model()
    coluna = JuMP.optimizer_with_attributes(
        Coluna.Optimizer,
        "params" => CL.Params(solver=ClA.SolveIpForm()),
        "default_optimizer" => GLPK.Optimizer
    )
    model = BlockModel(coluna)
    optimize!(model)
    @test JuMP.objective_value(model) == 0

    coluna = JuMP.optimizer_with_attributes(
        Coluna.Optimizer,
        "params" => CL.Params(solver=ClA.SolveLpForm(update_ip_primal_solution=true)),
        "default_optimizer" => GLPK.Optimizer
    )
    model = BlockModel(coluna)
    optimize!(model)
    @test JuMP.objective_value(model) == 0

    coluna = optimizer_with_attributes(
        Coluna.Optimizer,
        "params" => Coluna.Params(
            solver=Coluna.Algorithm.TreeSearchAlgorithm()
        ),
        "default_optimizer" => GLPK.Optimizer
    )
    model = BlockModel(coluna)
    @test_throws ErrorException optimize!(model)
end

function optimize_twice()
    # no reformulation + direct model
    coluna = JuMP.optimizer_with_attributes(
        Coluna.Optimizer,
        "params" => CL.Params(solver=ClA.SolveIpForm()),
        "default_optimizer" => GLPK.Optimizer
    )
    model = BlockModel(coluna, direct_model=true)
    @variable(model, x)
    @constraint(model, x <= 1)
    @objective(model, Max, x)
    optimize!(model)
    @test JuMP.objective_value(model) == 1
    optimize!(model)
    @test JuMP.objective_value(model) == 1

    # no reformulation + no direct model
    model = BlockModel(coluna)
    @variable(model, x)
    @constraint(model, x <= 1)
    @objective(model, Max, x)
    optimize!(model)
    @test JuMP.objective_value(model) == 1
    optimize!(model)
    @test JuMP.objective_value(model) == 1

    # reformulation + direct model
    data = CLD.GeneralizedAssignment.data("play2.txt")
    coluna = JuMP.optimizer_with_attributes(
        Coluna.Optimizer,
        "params" => CL.Params(solver=ClA.TreeSearchAlgorithm()),
        "default_optimizer" => GLPK.Optimizer
    )
    model, x, dec = CLD.GeneralizedAssignment.model(data, coluna)
    BD.objectiveprimalbound!(model, 100)
    BD.objectivedualbound!(model, 0)
    optimize!(model)
    @test JuMP.objective_value(model) ≈ 75.0
    optimize!(model)
    @test JuMP.objective_value(model) ≈ 75.0

    # reformulation + no direct model (`CLD.GeneralizedAssignment.model(data, coluna, false)` threw
    #                                  "MethodError: no method matching Model(; direct_model=false)")
    model = BlockModel(coluna)
    @axis(M, data.machines)
    @variable(model, x[m in M, j in data.jobs], Bin)
    @constraint(model, cov[j in data.jobs], sum(x[m,j] for m in M) >= 1)
    @constraint(model, knp[m in M], sum(data.weight[j,m] * x[m,j] for j in data.jobs) <= data.capacity[m])
    @objective(model, Min, sum(data.cost[j,m] * x[m,j] for m in M, j in data.jobs))
    @dantzig_wolfe_decomposition(model, dec, M)
    subproblems = BlockDecomposition.getsubproblems(dec)
    specify!.(subproblems, lower_multiplicity=0)
    BD.objectiveprimalbound!(model, 100)
    BD.objectivedualbound!(model, 0)
    optimize!(model)
    @test JuMP.objective_value(model) ≈ 75.0
    optimize!(model)
    @test JuMP.objective_value(model) ≈ 75.0
end

function column_generation_solver()
    data = CLD.GeneralizedAssignment.data("play2.txt")
    coluna = JuMP.optimizer_with_attributes(
        Coluna.Optimizer,
        "params" => CL.Params(solver=ClA.TreeSearchAlgorithm(
            maxnumnodes=1,
        )),
        "default_optimizer" => GLPK.Optimizer
    )
    treesearch, x, dec = CLD.GeneralizedAssignment.model_with_penalties(data, coluna)
    optimize!(treesearch)

    coluna = JuMP.optimizer_with_attributes(
        Coluna.Optimizer,
        "params" => Coluna.Params(solver=ClA.ColumnGeneration()),
        "default_optimizer" => GLPK.Optimizer
    )
    colgen, x, dec = CLD.GeneralizedAssignment.model_with_penalties(data, coluna)
    optimize!(colgen)
    
    @test MOI.get(treesearch, MOI.ObjectiveBound()) ≈ MOI.get(colgen, MOI.ObjectiveBound())
end

function branching_file_completion()
    function get_number_of_nodes_in_branching_tree_file(filename::String)
        filepath = string(@__DIR__, "/", filename)
        
        existing_nodes = Set()
        
        open(filepath) do file
            for line in eachline(file)
                for pieceofdata in split(line)
                    regex_match = match(r"n\d+", pieceofdata)
                    if regex_match != nothing
                        regex_match = regex_match.match
                        push!(existing_nodes, parse(Int, regex_match[2:length(regex_match)]))
                    end
                end
            end
        end
        return length(existing_nodes)
    end
    
                 
    data = CLD.GeneralizedAssignment.data("play2.txt")

    coluna = JuMP.optimizer_with_attributes(
        Coluna.Optimizer,
        "params" => CL.Params(solver=ClA.TreeSearchAlgorithm(
            branchingtreefile="playgap.dot"
        )),
        "default_optimizer" => GLPK.Optimizer
    )

    model, x, dec = CLD.GeneralizedAssignment.model(data, coluna)
    BD.objectiveprimalbound!(model, 100)
    BD.objectivedualbound!(model, 0)

    JuMP.optimize!(model)

    @test MOI.get(model, MOI.NodeCount()) == get_number_of_nodes_in_branching_tree_file("playgap.dot")
    @test JuMP.objective_value(model) ≈ 75.0
    @test JuMP.termination_status(model) == MOI.OPTIMAL
end

# Issue https://github.com/atoptima/Coluna.jl/issues/550
function continuous_vars_in_sp()
    # Simple min cost flow problem
    #   
    # n1 ------ f[1] -------> n3
    #  \                     ^
    #   \                   /
    #    -f[2]-> n2 --f[3]--
    # 
    #  n1: demand = -10.1
    #  n2: demand = 0
    #  n3: demand = 10.1 
    #  f[1]: cost = 0, capacity = 8.5, in mip model only integer flow allowed
    #  f[2]: cost = 50, (capacity = 5 can be activated by removing comment at constraint, line 93)
    #  f[3]: cost = 50
    #
    #  Correct solution for non-integer f[1]
    #           f[1] = 8.5, f[2] = f[3] = 1.6, cost = 8.5*0 + 1.6*2*50 = 160
    #  Correct solution for integer f[1]
    #           f[1] = 8, f[2] = f[3] = 2.1, cost = 8.5*0 + 2.1*2*50 = 210
    #
    function solve_flow_model(f1_integer, coluna)
        @axis(M, 1:1)
        model = BlockDecomposition.BlockModel(coluna, direct_model=true)
        @variable(model, f[1:3, m in M])
        if f1_integer
            JuMP.set_integer(f[1, 1])
        end
        @constraint(model, n1[m in M], f[1,m] + f[2,m] == 10.1)
        @constraint(model, n2[m in M], f[2,m] == f[3,m])
        @constraint(model, n3[m in M], f[1,m] + f[3,m] == 10.1)
        @constraint(model, cap1, sum(f[1,m] for m in M) <= 8.5)
        #@JuMP.constraint(model, cap2, sum(f[2,m] for m in M) <= 5)
        @objective(model, Min, 50 * f[2,1] + 50 * f[3,1])

        @dantzig_wolfe_decomposition(model, decomposition, M)

        subproblems = BlockDecomposition.getsubproblems(decomposition)
        BlockDecomposition.specify!.(subproblems, lower_multiplicity=1, upper_multiplicity=1)

        optimize!(model)

        if f1_integer
            @test termination_status(model) == MOI.OPTIMAL
            @test primal_status(model) == MOI.FEASIBLE_POINT
            @test objective_value(model) ≈ 210
            @test value(f[1,1]) ≈ 8
            @test value(f[2,1]) ≈ 2.1
            @test value(f[3,1]) ≈ 2.1
        else
            @test termination_status(model) == MOI.OPTIMAL
            @test primal_status(model) == MOI.FEASIBLE_POINT
            @test objective_value(model) ≈ 160
            @test value(f[1,1]) ≈ 8.5
            @test value(f[2,1]) ≈ 1.6
            @test value(f[3,1]) ≈ 1.6
        end
    end
    
    coluna = JuMP.optimizer_with_attributes(
        Coluna.Optimizer,
        "params" => Coluna.Params(
            solver=Coluna.Algorithm.TreeSearchAlgorithm(),
        ),
        "default_optimizer" => GLPK.Optimizer 
    );

    solve_flow_model(false, coluna)
    solve_flow_model(true, coluna)
    return
end

# Issue https://github.com/atoptima/Coluna.jl/issues/553
function anonymous_constraint_model()
    coluna = optimizer_with_attributes(
        Coluna.Optimizer,
        "params" => Coluna.Params(
            solver = Coluna.Algorithm.TreeSearchAlgorithm()
        ),
        "default_optimizer" => GLPK.Optimizer
    )
    
    model = BlockModel(coluna)
    
    @axis(S, 1:2)
    
    @variable(model, x[i in S], Bin)
    t = @constraint(model, sum(x[i] for i in S) <= 1)
    @objective(model, Min, -sum(x))
    
    @dantzig_wolfe_decomposition(model, decomposition, S)
    
    optimize!(model)

    @show BlockDecomposition.annotation(model, t)
    @test objectivevalue(model) ≈ -1
end

function test_issues_fixed()
    # @testset "no_decomposition" begin
    #     solve_with_no_decomposition()
    # end

    # @testset "moi_empty" begin
    #     test_model_empty()
    # end

    # @testset "decomposition_with_constant_in_objective" begin
    #     decomposition_with_constant_in_objective()
    # end

    # @testset "solve_empty_model" begin
    #     solve_empty_model()
    # end
    
    # @testset "optimize_twice" begin
    #     optimize_twice()
    # end

    # @testset "column_generation_solver" begin
    #     column_generation_solver()
    # end
    
    # @testset "branching_file_completion" begin
    #     branching_file_completion()
    # end

    # @testset "continuous_vars_in_sp" begin
    #     continuous_vars_in_sp()
    # end

    @testset "anonymous_constraint_model" begin
        anonymous_constraint_model()
        exit()
    end
    return
end

test_issues_fixed()