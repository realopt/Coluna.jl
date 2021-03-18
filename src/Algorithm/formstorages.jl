"""
    VarState

    Used in formulation records
"""

struct VarState
    cost::Float64
    lb::Float64
    ub::Float64
end

function apply_data!(form::Formulation, var::Variable, var_state::VarState)
    # Bounds
    if getcurlb(form, var) != var_state.lb || getcurub(form, var) != var_state.ub
        @logmsg LogLevel(-2) string("Reseting bounds of variable ", getname(form, var))
        setcurlb!(form, var, var_state.lb)
        setcurub!(form, var, var_state.ub)
        @logmsg LogLevel(-3) string("New lower bound is ", getcurlb(form, var))
        @logmsg LogLevel(-3) string("New upper bound is ", getcurub(form, var))
    end
    # Cost
    if getcurcost(form, var) != var_state.cost
        @logmsg LogLevel(-2) string("Reseting cost of variable ", getname(form, var))
        setcurcost!(form, var, var_state.cost)
        @logmsg LogLevel(-3) string("New cost is ", getcurcost(form, var))
    end
    return
end

"""
    ConstrState

    Used in formulation records
"""

struct ConstrState
    rhs::Float64
end

function apply_data!(form::Formulation, constr::Constraint, constr_state::ConstrState)
    # Rhs
    if getcurrhs(form, constr) != constr_state.rhs
        @logmsg LogLevel(-2) string("Reseting rhs of constraint ", getname(form, constr))
        setcurrhs!(form, constr, constr_state.rhs)
        @logmsg LogLevel(-3) string("New rhs is ", getcurrhs(form, constr))
    end
    return
end

"""
    FormulationRecord

    Formulation record is empty and it is used to implicitely keep 
    the data which is changed inside the model 
    (for example, dynamic variables and constraints of a formulaiton) 
    in order to store it to the record state and restore it afterwards. 
"""

struct FormulationRecord <: AbstractRecord end

FormulationRecord(form::Formulation) = FormulationRecord()

"""
    MasterBranchConstrsRecordPair

    Record pair for master branching constraints. 
    Consists of FormulationRecord and MasterBranchConstrsRecordState.    
"""

mutable struct MasterBranchConstrsRecordState <: AbstractRecordState
    constrs::Dict{ConstrId, ConstrState}
end

function Base.show(io::IO, state::MasterBranchConstrsRecordState)
    print(io, "[")
    for (id, constr) in state.constrs
        print(io, " ", MathProg.getuid(id))
    end
    print(io, "]")
end

function MasterBranchConstrsRecordState(form::Formulation, record::FormulationRecord)
    @logmsg LogLevel(-2) "Storing branching constraints"
    state = MasterBranchConstrsRecordState(Dict{ConstrId, ConstrState}())
    for (id, constr) in getconstrs(form)
        if getduty(id) <= AbstractMasterBranchingConstr && 
           iscuractive(form, constr) && isexplicit(form, constr)
            
            constrstate = ConstrState(getcurrhs(form, constr))
            state.constrs[id] = constrstate
        end
    end
    return state
end

function restorefromstate!(
    form::Formulation, record::FormulationRecord, state::MasterBranchConstrsRecordState
)
    @logmsg LogLevel(-2) "Restoring branching constraints"
    for (id, constr) in getconstrs(form)
        if getduty(id) <= AbstractMasterBranchingConstr && isexplicit(form, constr)
            @logmsg LogLevel(-4) "Checking " getname(form, constr)
            if haskey(state.constrs, id) 
                if !iscuractive(form, constr) 
                    @logmsg LogLevel(-2) string("Activating branching constraint", getname(form, constr))
                    activate!(form, constr)
                else    
                    @logmsg LogLevel(-2) string("Leaving branching constraint", getname(form, constr))
                end
                @logmsg LogLevel(-4) "Updating data"
                apply_data!(form, constr, state.constrs[id])
            else
                if iscuractive(form, constr) 
                    @logmsg LogLevel(-2) string("Deactivating branching constraint", getname(form, constr))
                    deactivate!(form, constr)
                end
            end    
        end
    end
end

const MasterBranchConstrsRecordPair = (FormulationRecord => MasterBranchConstrsRecordState)

"""
    MasterColumnsRecordPair

    Record pair for branching constraints of a formulation. 
    Consists of EmptyRecord and MasterColumnsState.    
"""

mutable struct MasterColumnsState <: AbstractRecordState
    cols::Dict{VarId, VarState}
end

function Base.show(io::IO, state::MasterColumnsState)
    print(io, "[")
    for (id, val) in state.cols
        print(io, " ", MathProg.getuid(id))
    end
    print(io, "]")
end

function MasterColumnsState(form::Formulation, record::FormulationRecord)
    @logmsg LogLevel(-2) "Storing master columns"
    state = MasterColumnsState(Dict{VarId, ConstrState}())
    for (id, var) in getvars(form)
        if getduty(id) <= MasterCol && 
           iscuractive(form, var) && isexplicit(form, var)
            
            varstate = VarState(getcurcost(form, var), getcurlb(form, var), getcurub(form, var))
            state.cols[id] = varstate
        end
    end
    return state
end

function restorefromstate!(
    form::Formulation, record::FormulationRecord, state::MasterColumnsState
)
    @logmsg LogLevel(-2) "Restoring master columns"
    for (id, var) in getvars(form)
        if getduty(id) <= MasterCol && isexplicit(form, var)
            @logmsg LogLevel(-4) "Checking " getname(form, var)
            if haskey(state.cols, id) 
                if !iscuractive(form, var) 
                    @logmsg LogLevel(-4) string("Activating column", getname(form, var))
                    activate!(form, var)
                end
                @logmsg LogLevel(-4) "Updating data"
                apply_data!(form, var, state.cols[id])
            else
                if iscuractive(form, var) 
                    @logmsg LogLevel(-4) string("Deactivating column", getname(form, var))
                    deactivate!(form, var)
                end
            end    
        end
    end
end

const MasterColumnsRecordPair = (FormulationRecord => MasterColumnsState)

"""
    MasterCutsRecordPair

    Record pair for cutting planes of a formulation. 
    Consists of EmptyRecord and MasterCutsState.    
"""

mutable struct MasterCutsState <: AbstractRecordState
    cuts::Dict{ConstrId, ConstrState}
end

function Base.show(io::IO, state::MasterCutsState)
    print(io, "[")
    for (id, constr) in state.cuts
        print(io, " ", MathProg.getuid(id))
    end
    print(io, "]")
end

function MasterCutsState(form::Formulation, record::FormulationRecord)
    @logmsg LogLevel(-2) "Storing master cuts"
    state = MasterCutsState(Dict{ConstrId, ConstrState}())
    for (id, constr) in getconstrs(form)
        if getduty(id) <= AbstractMasterCutConstr && 
           iscuractive(form, constr) && isexplicit(form, constr)
            
            constrstate = ConstrState(getcurrhs(form, constr))
            state.cuts[id] = constrstate
        end
    end
    return state
end

function restorefromstate!(
    form::Formulation, record::FormulationRecord, state::MasterCutsState
)
    @logmsg LogLevel(-2) "Storing master cuts"
    for (id, constr) in getconstrs(form)
        if getduty(id) <= AbstractMasterCutConstr && isexplicit(form, constr)
            @logmsg LogLevel(-4) "Checking " getname(form, constr)
            if haskey(state.cuts, id) 
                if !iscuractive(form, constr) 
                    @logmsg LogLevel(-4) string("Activating cut", getname(form, constr))
                    activate!(form, constr)
                end
                @logmsg LogLevel(-4) "Updating data"
                apply_data!(form, constr, state.cuts[id])
            else
                if iscuractive(form, constr) 
                    @logmsg LogLevel(-4) string("Deactivating cut", getname(form, constr))
                    deactivate!(form, constr)
                end
            end    
        end
    end
end

const MasterCutsRecordPair = (FormulationRecord => MasterCutsState)

"""
    StaticVarConstrRecordPair

    Record pair for static variables and constraints of a formulation.
    Consists of EmptyRecord and StaticVarConstrRecordState.    
"""

mutable struct StaticVarConstrRecordState <: AbstractRecordState
    constrs::Dict{ConstrId, ConstrState}
    vars::Dict{VarId, VarState}
end

#TO DO: we need to keep here only the difference with the initial data

function Base.show(io::IO, state::StaticVarConstrRecordState)
    print(io, "[vars:")
    for (id, var) in state.vars
        print(io, " ", MathProg.getuid(id))
    end
    print(io, ", constrs:")
    for (id, constr) in state.constrs
        print(io, " ", MathProg.getuid(id))
    end
    print(io, "]")
end

function StaticVarConstrRecordState(form::Formulation, record::FormulationRecord)
    @logmsg LogLevel(-2) string("Storing static vars and consts")
    state = StaticVarConstrRecordState(Dict{ConstrId, ConstrState}(), Dict{VarId, VarState}())
    for (id, constr) in getconstrs(form)
        if isaStaticDuty(getduty(id)) && iscuractive(form, constr) && isexplicit(form, constr)            
            constrstate = ConstrState(getcurrhs(form, constr))
            state.constrs[id] = constrstate
        end
    end
    for (id, var) in getvars(form)
        if isaStaticDuty(getduty(id)) && iscuractive(form, var) && isexplicit(form, var)            
            varstate = VarState(getcurcost(form, var), getcurlb(form, var), getcurub(form, var))
            state.vars[id] = varstate
        end
    end
    return state
end

function restorefromstate!(
    form::Formulation, record::FormulationRecord, state::StaticVarConstrRecordState
)
    @logmsg LogLevel(-2) "Restoring static vars and consts"
    for (id, constr) in getconstrs(form)
        if isaStaticDuty(getduty(id)) && isexplicit(form, constr)
            @logmsg LogLevel(-4) "Checking " getname(form, constr)
            if haskey(state.constrs, id) 
                if !iscuractive(form, constr) 
                    @logmsg LogLevel(-2) string("Activating constraint", getname(form, constr))
                    activate!(form, constr)
                end
                @logmsg LogLevel(-4) "Updating data"
                apply_data!(form, constr, state.constrs[id])
            else
                if iscuractive(form, constr) 
                    @logmsg LogLevel(-2) string("Deactivating constraint", getname(form, constr))
                    deactivate!(form, constr)
                end
            end    
        end
    end
    for (id, var) in getvars(form)
        if isaStaticDuty(getduty(id)) && isexplicit(form, var)
            @logmsg LogLevel(-4) "Checking " getname(form, var)
            if haskey(state.vars, id) 
                if !iscuractive(form, var) 
                    @logmsg LogLevel(-4) string("Activating variable", getname(form, var))
                    activate!(form, var)
                end
                @logmsg LogLevel(-4) "Updating data"
                apply_data!(form, var, state.vars[id])
            else
                if iscuractive(form, var) 
                    @logmsg LogLevel(-4) string("Deactivating variable", getname(form, var))
                    deactivate!(form, var)
                end
            end    
        end
    end
end

const StaticVarConstrRecordPair = (FormulationRecord => StaticVarConstrRecordState)

"""
    PartialSolutionRecordPair

    Record pair for partial solution of a formulation.
    Consists of PartialSolutionRecord and PartialSolutionRecordState.    
"""

# TO DO : to replace dictionaries by PrimalSolution
# issues to see : 1) PrimalSolution is parametric; 2) we need a solution concatenation functionality

mutable struct PartialSolutionRecord <: AbstractRecord
    solution::Dict{VarId, Float64}
end

function add_to_solution!(record::PartialSolutionRecord, varid::VarId, value::Float64)
    cur_value = get(record.solution, varid, 0.0)
    record.solution[varid] = cur_value + value
    return
end

function get_primal_solution(record::PartialSolutionRecord, form::Formulation)
    varids = collect(keys(record.solution))
    vals = collect(values(record.solution))
    solcost = 0.0
    for (varid, value) in record.solution
        solcost += getcurcost(form, varid) * value
    end
    return PrimalSolution(form, varids, vals, solcost, UNKNOWN_FEASIBILITY)
end    


function PartialSolutionRecord(form::Formulation) 
    return PartialSolutionRecord(Dict{VarId, Float64}())
end

# the record state is the same as the record here
mutable struct PartialSolutionRecordState <: AbstractRecordState
    solution::Dict{VarId, Float64}
end

function PartialSolutionRecordState(form::Formulation, record::PartialSolutionRecord)
    @logmsg LogLevel(-2) "Storing partial solution"
    return PartialSolutionRecordState(copy(record.solution))
end

function restorefromstate!(
    form::Formulation, record::PartialSolutionRecord, state::PartialSolutionRecordState
)
    @logmsg LogLevel(-2) "Restoring partial solution"
    record.solution = copy(state.solution)
end

const PartialSolutionRecordPair = (PartialSolutionRecord => PartialSolutionRecordState)
