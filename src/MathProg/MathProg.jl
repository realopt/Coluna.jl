module MathProg

import BlockDecomposition
import Distributed
import MathOptInterface
import TimerOutputs

import ..Coluna # for NestedEnum (types.jl:210)
using ..ColunaBase

import Base: haskey, length, iterate, diff

using DynamicSparseArrays, Logging, Printf

if VERSION >= v"1.5"
    import Base: contains
end

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
include("manager.jl")
include("bounds.jl")
include("solutions.jl")
include("buffer.jl")
include("formulation.jl")
include("varconstr.jl")
include("optimizerwrappers.jl")
include("clone.jl")
include("reformulation.jl")
include("projection.jl")
include("problem.jl")
include("MOIinterface.jl")

# TODO : clean up
# Types
export  MaxSense, MinSense, MoiOptimizer,
        Id, ConstrSense, VarSense,
        FormId, FormulationPhase, Annotations,
        Counter, UserOptimizer, MoiObjective

# Methods
export no_optimizer_builder, set_original_formulation!,
       getid, getuid,
       enforce_integrality!, relax_integrality!,
       getobjsense, getoptimizer,
       setdualbound!,
       computereducedcost,
       update!,
       convert_status, getduty,
       computereducedrhs,
       find_owner_formulation,
       getsortuid,
       contains, setprimalbound!, get_original_formulation,
       getoriginformuid, getspsol, sync_solver!, getinner,
       get_primal_solutions, get_dual_solutions

# Below this line, clean up has been done :

# Methods related to Problem
export Problem, set_initial_dual_bound!, set_initial_primal_bound!,
       get_initial_dual_bound, get_initial_primal_bound

# Methods related to Reformulation
export Reformulation, getmaster, add_dw_pricing_sp!, add_benders_sep_sp!, get_dw_pricing_sps,
    set_reformulation!, get_benders_sep_sps, get_dw_pricing_sp_ub_constrid,
    get_dw_pricing_sp_lb_constrid, setmaster!

# Methods related to formulations
export AbstractFormulation, Formulation, getreformulation, getvar, getvars, getconstr,
    getconstrs, getelem, getcoefmatrix, getprimalsolmatrix, getprimalsolcosts,
    getdualsolmatrix, getdualsolrhss, setvar!, setconstr!, setprimalsol!, setdualsol!,
    set_robust_constr_generator!, get_robust_constr_generators,
    setcol_from_sp_primalsol!, setcut_from_sp_dualsol!, # TODO : merge with setvar! & setconstr
    set_objective_sense!, clonevar!, cloneconstr!, clonecoeffs!, initialize_optimizer!

# Duties of formulations
export Original, DwMaster, BendersMaster, DwSp, BendersSp

# Methods related to duties
export isanArtificialDuty, isaStaticDuty, isaDynamicDuty, isanOriginalRepresentatives

# Types and methods related to variables and constraints
export Variable, Constraint, VarId, ConstrId, VarMembership, ConstrMembership,
    getperencost, getcurcost, setcurcost!, getperenlb, getcurlb, setcurlb!,
    getperenub, getcurub, setcurub!, getperenrhs, getcurrhs, setcurrhs!, getperensense,
    getcursense, setcursense!, getperenkind, getcurkind, setcurkind!, getperenincval,
    getcurincval, setcurincval!, isperenactive, iscuractive, activate!, deactivate!,
    isexplicit, getname, reset!, getreducedcost

# Types & methods related to solutions & bounds
export PrimalBound, DualBound, PrimalSolution, DualSolution, ObjValues,
    get_ip_primal_bound, get_lp_primal_bound,
    get_ip_dual_bound, get_lp_dual_bound, update_ip_primal_bound!, update_lp_primal_bound!,
    update_ip_dual_bound!, update_lp_dual_bound!, set_ip_primal_bound!,
    set_lp_primal_bound!, set_ip_dual_bound!, set_lp_dual_bound!, ip_gap, lp_gap, ip_gap_closed, 
    lp_gap_closed

# Methods related to projections
export projection_is_possible, proj_cols_on_rep

end
