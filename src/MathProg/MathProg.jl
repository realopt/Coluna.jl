module MathProg

import BlockDecomposition
import Distributed
import MathOptInterface
import TimerOutputs

import ..Coluna # for NestedEnum (types.jl:210)
using ..Coluna: iterate # to be deleted
using ..ColunaBase

import Base: haskey, length, iterate, diff

using DynamicSparseArrays
using Logging
using Printf

global const BD = BlockDecomposition
global const MOI = MathOptInterface
global const MOIU = MathOptInterface.Utilities
global const TO = TimerOutputs
const Cl = Coluna

include("counters.jl")
include("types.jl")
include("vcids.jl")
include("variable.jl")
include("constraint.jl")
include("duties.jl")
include("varconstr.jl") # to rm
include("manager.jl")
include("bounds.jl")
include("solutions.jl")
include("incumbents.jl") # to rm
include("optimizationresults.jl")
include("buffer.jl")
include("formulation.jl")
include("new_varconstr.jl") 
include("optimizerwrappers.jl")
include("clone.jl")
include("reformulation.jl")
include("projection.jl")
include("problem.jl")
include("decomposition.jl")
include("MOIinterface.jl")

# TODO : clean up
# Types
export  MaxSense, MinSense, MoiOptimizer,
         Id, ConstrSense, VarSense,
        FormId, FormulationPhase, Annotations,
       Original, Counter, UserOptimizer, MoiObjective

# Methods
export no_optimizer_builder, set_original_formulation!,
       getid, store!, getuid,
       register_objective_sense!, nbprimalsols, getdualbound,
       getprimalbound,
       enforce_integrality!, relax_integrality!, 
       getobjsense, getoptimizer, getbestprimalsol,
       add_primal_sol!, setdualbound!,
       getprimalsols, getdualsols,
       computereducedcost,
       update!,
       convert_status, getduty, getbestdualsol, 
       computereducedrhs, 
       unsafe_getbestprimalsol,
        find_owner_formulation,
        get_dw_pricing_sps, 
       getsortuid,
       get_benders_sep_sps,
       contains, setprimalbound!

# Below this line, clean up has been done :
export reformulate!, optimize!, setwarmstart!

# Methods related to Problem
export Problem, set_initial_dual_bound!, set_initial_primal_bound!,
       get_initial_dual_bound, get_initial_primal_bound

# Methods related to formulations
export AbstractFormulation, Reformulation, Formulation, getmaster, getreformulation,
       getvar, getvars, getconstr, getconstrs, getelem,
       getcoefmatrix,
       getprimalsolmatrix,
       getprimalsolcosts,
       getdualsolmatrix,
       getdualsolrhss,
       setvar!, setconstr!,
       setprimalsol!, setdualsol!,
       setcol_from_sp_primalsol!, setcut_from_sp_dualsol! # TODO : merge with setvar! & setconstr!

# Methods related to duties
export isanArtificialDuty, isaStaticDuty, isaDynamicDuty, isanOriginalRepresentatives

# Types and methods related to variables and constraints
export Variable, Constraint, VarId, ConstrId, VarMembership, ConstrMembership, 
    getperenecost, getcurcost, setcurcost!, getperenelb, getcurlb, setcurlb!, 
    getpereneub, getcurub, setcurub!, getperenerhs, getcurrhs, setcurrhs!, getperenesense,
    getcursense, setcursense!, getperenekind, getcurkind, setcurkind!, getpereneincval, 
    getcurincval, setcurincval!, ispereneactive, iscuractive, activate!, deactivate!,
    ispereneexplicit, iscurexplicit, setiscurexplicit!, getname, reset!

# Types & methods related to solutions & bounds
# Note : we should export only get methods for MoiResult (the solution is built in MathProg)
export PrimalBound, DualBound, PrimalSolution, DualSolution, ObjValues, TerminationStatus,
    FeasibilityStatus, MoiResult, OPTIMAL, TIME_LIMIT, NODE_LIMIT, OTHER_LIMIT, 
    EMPTY_RESULT, NOT_YET_DETERMINED, INFEASIBLE, UNKNOWN_FEASIBILITY, FEASIBLE,
    getterminationstatus, getfeasibilitystatus, setterminationstatus!, 
    setfeasibilitystatus!, isfeasible, get_ip_primal_bound, get_lp_primal_bound,
    get_ip_dual_bound, get_lp_dual_bound, update_ip_primal_bound!, update_lp_primal_bound!,
    update_ip_dual_bound!, update_lp_dual_bound!, set_ip_primal_bound!, 
    set_lp_primal_bound!, set_ip_dual_bound!, set_lp_dual_bound!, ip_gap

# methods related to projections
export projection_is_possible, proj_cols_on_rep

# convert methods
export convert_coluna_sense_to_moi,
       convert_moi_sense_to_coluna,
       convert_moi_rhs_to_coluna,
       convert_moi_kind_to_coluna

end
