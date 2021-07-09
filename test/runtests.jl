using Base.CoreLogging: Error
using Coluna

using Test, GLPK, ColunaDemos, JuMP, BlockDecomposition
using Random

using MathOptInterface, MathOptInterface.Utilities

using Base.CoreLogging, Logging
global_logger(ConsoleLogger(stderr, LogLevel(0)))

const MOI = MathOptInterface
const MOIU = MOI.Utilities
const MOIT = MOI.Test
const MOIB = MOI.Bridges

const CL = Coluna
const CLD = ColunaDemos
const BD = BlockDecomposition

const ClF = Coluna.MathProg # Must be deleted
const ClMP = Coluna.MathProg
const ClA = Coluna.Algorithm

# include("unit/unit_tests.jl")
# include("MathOptInterface/MOI_wrapper.jl")
# include("interfaces/model.jl")
include("issues_tests.jl")
include("show_functions_tests.jl")
include("full_instances_tests.jl")
include("user_algorithms_tests.jl")
include("preprocessing_tests.jl")
include("pricing_callback_tests.jl")
include("bound_callback_tests.jl")
include("optimizer_with_attributes_test.jl")
include("subproblem_solvers_tests.jl")
include("custom_var_cuts_tests.jl")
include("sol_disaggregation_tests.jl")

rng = MersenneTwister(1234123)

unit_tests()

@testset "Full instances " begin
    full_instances_tests()
end

@testset "User algorithms" begin
    user_algorithms_tests()
end

# @testset "Preprocessing " begin
#     preprocessing_tests()
# end

@testset "pricing callback" begin
    pricing_callback_tests()
end

@testset "bound callback" begin
    bound_callback_tests()
end

@testset "Base.show functions " begin
    backup_stdout = stdout
    (rd_out, wr_out) = redirect_stdout()
    show_functions_tests()
    close(wr_out)
    close(rd_out)
    redirect_stdout(backup_stdout)
end

@testset "Optimizer with Attributes" begin
    optimizer_with_attributes_test()
end

@testset "Subproblem Solvers" begin
    subproblem_solvers_test()
end

@testset "Custom Variables and Cuts" begin
    custom_var_cuts_test()
end

@testset "Solution Disaggregation" begin
    sol_disaggregation_tests()
end
