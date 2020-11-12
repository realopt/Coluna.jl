module Algorithm

import DataStructures
import MathOptInterface
import TimerOutputs
import DynamicSparseArrays

using ..Coluna, ..ColunaBase, ..MathProg

using Logging, Parameters, Printf, Statistics

const TO = TimerOutputs
const DS = DataStructures
const MOI = MathOptInterface

import Base: push!

# Import to extend methods to OptimizationState
import ..MathProg:
    get_ip_primal_bound, get_ip_dual_bound,
    get_lp_primal_bound, get_lp_dual_bound, update_ip_primal_bound!, update_ip_dual_bound!,
    update_lp_primal_bound!, update_lp_dual_bound!, set_ip_primal_bound!,
    set_ip_dual_bound!, set_lp_primal_bound!, set_lp_dual_bound!, ip_gap, lp_gap, ip_gap_closed, lp_gap_closed

# Utilities to build algorithms
include("utilities/optimizationstate.jl")

include("storage.jl")
include("data.jl")
include("formstorages.jl")

# Abstract algorithm
include("interface.jl")

# Basic algorithms
include("basic/solveipform.jl")
include("basic/solvelpform.jl")
include("basic/cutcallback.jl")

# Child algorithms used by conquer algorithms
include("colgenstabilization.jl")
include("colgen.jl")
include("benders.jl")
include("preprocessing.jl")

# Algorithms and structures used by the tree search algorithm
include("node.jl")
include("conquer.jl")
include("divide.jl")

# Here include divide algorithms
include("branching/branchinggroup.jl")
include("branching/branchingrule.jl")
include("branching/varbranching.jl")
include("branching/branchingalgo.jl")

include("treesearch.jl")

# Algorithm should export only methods usefull to define & parametrize algorithms, and
# data structures from utilities.
# Other Coluna's submodules should be independent to Algorithm

# Utilities
export getterminationstatus, setterminationstatus!, 
    nb_ip_primal_sols, nb_lp_primal_sols, nb_lp_dual_sols,
    get_ip_primal_sols, get_lp_primal_sols, get_lp_dual_sols, get_best_ip_primal_sol,
    get_best_lp_primal_sol, get_best_lp_dual_sol, update_ip_primal_sol!,
    update_lp_primal_sol!, update_lp_dual_sol!, add_ip_primal_sol!, add_lp_primal_sol!,
    add_lp_dual_sol!, set_ip_primal_sol!, set_lp_primal_sol!, set_lp_dual_sol!,
    get_ip_dual_bound, set_ip_dual_bound!, update_all_ip_primal_solutions!, getreform,
    getmasterdata, getoptstate, run!, isinfeasible
    
# Algorithm's types
export AbstractOptimizationAlgorithm, TreeSearchAlgorithm, ColCutGenConquer, ColumnGeneration,
       BendersConquer, BendersCutGeneration, SolveIpForm, SolveLpForm, ExactBranchingPhase,
       OnlyRestrictedMasterBranchingPhase, PreprocessAlgorithm, RestrictedMasterIPHeuristic,
       OptimizationInput, OptimizationOutput, OptimizationState, ModelData, ReformData,
       EmptyInput       

# Storages 
export PartialSolutionStoragePair, PreprocessingStoragePair       

# Storage functions 
export getstorage, add_to_solution!, add_to_localpartialsol!

end
