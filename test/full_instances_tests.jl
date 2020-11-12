using LightGraphs

function full_instances_tests()
    generalized_assignment_tests()
    capacitated_lot_sizing_tests()
    lot_sizing_tests()
    #facility_location_tests()
    cutting_stock_tests()
    cvrp_tests()
end

function generalized_assignment_tests()
    @testset "play gap" begin
        data = CLD.GeneralizedAssignment.data("play2.txt")

        coluna = JuMP.optimizer_with_attributes(
            Coluna.Optimizer,
            "params" => CL.Params(solver = ClA.TreeSearchAlgorithm(
                branchingtreefile = "playgap.dot"
            )),
            "default_optimizer" => GLPK.Optimizer
        )

        model, x, dec = CLD.GeneralizedAssignment.model(data, coluna)
        BD.objectiveprimalbound!(model, 100)
        BD.objectivedualbound!(model, 0)

        JuMP.optimize!(model)

        @test JuMP.objective_value(model) ≈ 75.0
        @test JuMP.termination_status(model) == MOI.OPTIMAL
        @test CLD.GeneralizedAssignment.print_and_check_sol(data, model, x)
        @test MOI.get(model, MOI.NumberOfVariables()) == length(x)
        @test MOI.get(model, MOI.SolverName()) == "Coluna"
    end

    @testset "gap - JuMP/MOI modeling" begin
        data = CLD.GeneralizedAssignment.data("smallgap3.txt")

        coluna = JuMP.optimizer_with_attributes(
            Coluna.Optimizer,
            "params" => CL.Params(solver = ClA.TreeSearchAlgorithm()),
            "default_optimizer" => GLPK.Optimizer
        )

        model, x, dec = CLD.GeneralizedAssignment.model(data, coluna)
        BD.objectiveprimalbound!(model, 500.0)
        BD.objectivedualbound!(model, 0.0)

        JuMP.optimize!(model)
        @test JuMP.objective_value(model) ≈ 438.0
        @test JuMP.termination_status(model) == MOI.OPTIMAL
        @test CLD.GeneralizedAssignment.print_and_check_sol(data, model, x)
    end

    @testset "gap - strong branching" begin
        data = CLD.GeneralizedAssignment.data("mediumgapcuts3.txt")

        conquer_with_small_cleanup_threshold = ClA.ColCutGenConquer(
            colgen = ClA.ColumnGeneration(cleanup_threshold = 150, smoothing_stabilization = 1.0)
        )

        branching = ClA.StrongBranching()
        push!(branching.phases, ClA.BranchingPhase(5, ClA.RestrMasterLPConquer()))
        push!(branching.phases, ClA.BranchingPhase(1, conquer_with_small_cleanup_threshold))
        push!(branching.rules, ClA.PrioritisedBranchingRule(ClA.VarBranchingRule(), 1.0, 1.0))

        coluna = JuMP.optimizer_with_attributes(
            CL.Optimizer,
            "params" => CL.Params(
                solver = ClA.TreeSearchAlgorithm(
                    conqueralg = conquer_with_small_cleanup_threshold,
                    dividealg = branching,
                    maxnumnodes = 300
                )
            ),
            "default_optimizer" => GLPK.Optimizer
        )

        model, x, dec = CLD.GeneralizedAssignment.model(data, coluna)

        BD.objectiveprimalbound!(model, 2000.0)
        BD.objectivedualbound!(model, 0.0)

        JuMP.optimize!(model)

        @test JuMP.objective_value(model) ≈ 1553.0
        @test JuMP.termination_status(model) == MOI.OPTIMAL
        @test CLD.GeneralizedAssignment.print_and_check_sol(data, model, x)
    end

    @testset "gap - ColGen max nb iterations" begin
        data = CLD.GeneralizedAssignment.data("smallgap3.txt")

        coluna = JuMP.optimizer_with_attributes(
            CL.Optimizer,
            "params" => CL.Params(
                solver = ClA.TreeSearchAlgorithm(
                    conqueralg = ClA.ColCutGenConquer(
                        colgen = ClA.ColumnGeneration(max_nb_iterations = 8)
                    )
                )
            ),
            "default_optimizer" => GLPK.Optimizer
        )

        problem, x, dec = CLD.GeneralizedAssignment.model(data, coluna)

        JuMP.optimize!(problem)
        @test abs(JuMP.objective_value(problem) - 438.0) <= 0.00001
        @test JuMP.termination_status(problem) == MOI.OPTIMAL # Problem with final dual bound ?
        @test CLD.GeneralizedAssignment.print_and_check_sol(data, problem, x)
    end

    @testset "gap with penalties - pure master variables" begin
        data = CLD.GeneralizedAssignment.data("smallgap3.txt")

        coluna = JuMP.optimizer_with_attributes(
            Coluna.Optimizer,
            "params" => CL.Params(solver = ClA.TreeSearchAlgorithm()),
            "default_optimizer" => GLPK.Optimizer
        )

        problem, x, y, dec = CLD.GeneralizedAssignment.model_with_penalties(data, coluna)
        JuMP.optimize!(problem)
        @test JuMP.termination_status(problem) == MOI.OPTIMAL
        @test abs(JuMP.objective_value(problem) - 416.4) <= 0.00001
    end

    @testset "gap with maximisation objective function" begin
        data = CLD.GeneralizedAssignment.data("smallgap3.txt")

        coluna = JuMP.optimizer_with_attributes(
            Coluna.Optimizer,
            "params" => CL.Params(solver = ClA.TreeSearchAlgorithm()),
            "default_optimizer" => GLPK.Optimizer
        )

        problem, x, dec = CLD.GeneralizedAssignment.model_max(data, coluna)
        JuMP.optimize!(problem)
        @test JuMP.termination_status(problem) == MOI.OPTIMAL
        @test abs(JuMP.objective_value(problem) - 580.0) <= 0.00001
    end

    @testset "gap with infeasible master" begin
        data = CLD.GeneralizedAssignment.data("master_infeas.txt")

        coluna = JuMP.optimizer_with_attributes(
            Coluna.Optimizer,
            "params" => CL.Params(solver = ClA.TreeSearchAlgorithm()),
            "default_optimizer" => GLPK.Optimizer
        )

        problem, x, dec = CLD.GeneralizedAssignment.model(data, coluna)

        JuMP.optimize!(problem)
        @test JuMP.termination_status(problem) == MOI.INFEASIBLE
    end

    @testset "gap with infeasible subproblem" begin
        data = CLD.GeneralizedAssignment.data("sp_infeas.txt")

        coluna = JuMP.optimizer_with_attributes(
            Coluna.Optimizer,
            "params" => CL.Params(solver = ClA.TreeSearchAlgorithm()),
            "default_optimizer" => GLPK.Optimizer
        )

        problem, x, dec = CLD.GeneralizedAssignment.model(data, coluna)

        JuMP.optimize!(problem)
        @test JuMP.termination_status(problem) == MOI.INFEASIBLE
    end

    @testset "gap with all phases in col.gen" begin
        data = CLD.GeneralizedAssignment.data("mediumgapcuts1.txt")
        for m in data.machines
            data.capacity[m] = floor(Int, data.capacity[m] * 0.5)
        end

        coluna = JuMP.optimizer_with_attributes(
            Coluna.Optimizer,
            "params" => CL.Params(solver = ClA.TreeSearchAlgorithm(
                conqueralg = ClA.ColCutGenConquer(
                    colgen = ClA.ColumnGeneration(opt_rtol = 1e-4, smoothing_stabilization = 0.5)
                )
            )),
            "default_optimizer" => GLPK.Optimizer
        )

        problem, x, y, dec = CLD.GeneralizedAssignment.model_with_penalty(data, coluna)

        JuMP.optimize!(problem)
        @test abs(JuMP.objective_value(problem) - 31895.0) <= 0.00001
    end

    @testset "gap with max. obj., pure mast. vars., and stabilization" begin
        data = CLD.GeneralizedAssignment.data("gapC-5-100.txt")

        coluna = JuMP.optimizer_with_attributes(
            CL.Optimizer,
            "params" => CL.Params(
                solver = ClA.TreeSearchAlgorithm(
                    conqueralg = ClA.ColCutGenConquer(
                        colgen = ClA.ColumnGeneration(smoothing_stabilization = 1.0)
                    ),
                    maxnumnodes = 300
                )
            ),
            "default_optimizer" => GLPK.Optimizer
        )

        model, x, y, dec = CLD.GeneralizedAssignment.max_model_with_subcontracts(data, coluna)

        JuMP.optimize!(model)

        @test JuMP.objective_value(model) ≈ 3520.1
        @test JuMP.termination_status(model) == MOI.OPTIMAL
    end

    @testset "play gap" begin
        data = CLD.GeneralizedAssignment.data("play2.txt")

        coluna = JuMP.optimizer_with_attributes(
            Coluna.Optimizer,
            "params" => CL.Params(solver = ClA.TreeSearchAlgorithm()),
            "default_optimizer" => GLPK.Optimizer
        )

        problem, x, dec = CLD.GeneralizedAssignment.model(data, coluna)
        JuMP.optimize!(problem)
        @test abs(JuMP.objective_value(problem) - 75.0) <= 0.00001
        @test JuMP.termination_status(problem) == MOI.OPTIMAL
        @test CLD.GeneralizedAssignment.print_and_check_sol(data, problem, x)
    end

    @testset "play gap with no solver" begin
        data = CLD.GeneralizedAssignment.data("play2.txt")

        coluna = JuMP.optimizer_with_attributes(
            Coluna.Optimizer,
            "params" => CL.Params(solver = ClA.TreeSearchAlgorithm())
        )

        problem, x, dec = CLD.GeneralizedAssignment.model(data, coluna)
        try
            JuMP.optimize!(problem)
        catch e
            @test repr(e) == "ErrorException(\"Cannot optimize LP formulation with optimizer of type Coluna.MathProg.NoOptimizer.\")"
        end
    end
    return
end

function lot_sizing_tests()
    @testset "play single mode multi items lot sizing" begin
        data = CLD.SingleModeMultiItemsLotSizing.data("lotSizing-3-20-2.txt")

        coluna = JuMP.optimizer_with_attributes(
            Coluna.Optimizer,
            "params" => CL.Params(
                solver = ClA.BendersCutGeneration()
            ),
            "default_optimizer" => GLPK.Optimizer
        )

        problem, x, y, dec = CLD.SingleModeMultiItemsLotSizing.model(data, coluna)
        JuMP.optimize!(problem)
        @test 600 - 1e-6 <= objective_value(problem) <= 600 + 1e-6
    end
    return
end

function capacitated_lot_sizing_tests()
    @testset "clsp small instance" begin
        data = CLD.CapacitatedLotSizing.readData("testSmall")

        coluna = JuMP.optimizer_with_attributes(
            Coluna.Optimizer,
            "params" => CL.Params(solver = ClA.TreeSearchAlgorithm()),
            "default_optimizer" => GLPK.Optimizer
        )

        model, x, y, s, dec = CLD.CapacitatedLotSizing.model(data, coluna)
        JuMP.optimize!(model)

        @test JuMP.termination_status(model) == MOI.OPTIMAL
    end
end

function facility_location_tests()
    @testset "play facility location test " begin
        data = CLD.FacilityLocation.data("play.txt")

        coluna = JuMP.optimizer_with_attributes(
            Coluna.Optimizer,
            "params" => CL.Params(
                solver = ClA.BendersCutGeneration()
            ),
            "default_optimizer" => GLPK.Optimizer
        )

        problem, x, y, dec = CLD.FacilityLocation.model(data, coluna)
        JuMP.optimize!(problem)
    end
    return
end

function cutting_stock_tests()
    @testset "play cutting stock" begin
        data = CLD.CuttingStock.data("randomInstances/inst10-10")

        coluna = JuMP.optimizer_with_attributes(
            Coluna.Optimizer,
            "params" => CL.Params(solver = ClA.TreeSearchAlgorithm()),
            "default_optimizer" => GLPK.Optimizer
        )

        problem, x, y, dec = CLD.CuttingStock.model(data, coluna)
        JuMP.optimize!(problem)
        @test 4 - 1e-6 <= objective_value(problem) <= 4 + 1e-6
    end
    return
end

function cvrp_tests()
    @testset "play cvrp" begin
        data = CLD.CapacitatedVehicleRouting.data("A-n16-k3.vrp")

        coluna = JuMP.optimizer_with_attributes(
            Coluna.Optimizer,
            "params" => CL.Params(solver = ClA.TreeSearchAlgorithm(
                maxnumnodes = 10000,
                branchingtreefile = "cvrp.dot"
            )),
            "default_optimizer" => GLPK.Optimizer
        )

        model, x, dec = CLD.CapacitatedVehicleRouting.model(data, coluna)
        JuMP.optimize!(model)
        @test objective_value(model) ≈ 504
    end
    return
end
