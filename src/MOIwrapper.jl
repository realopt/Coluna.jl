const CleverDicts = MOI.Utilities.CleverDicts

const SupportedObjFunc = Union{MOI.ScalarAffineFunction{Float64}, MOI.SingleVariable}

const SupportedVarSets = Union{
    MOI.ZeroOne, MOI.Integer, MOI.LessThan{Float64}, MOI.EqualTo{Float64}, 
    MOI.GreaterThan{Float64}
}

const SupportedConstrFunc = Union{MOI.ScalarAffineFunction{Float64}}

const SupportedConstrSets = Union{
    MOI.EqualTo{Float64}, MOI.GreaterThan{Float64}, MOI.LessThan{Float64}
}

@enum(ObjectiveType, SINGLE_VARIABLE, SCALAR_AFFINE)
mutable struct Optimizer <: MOI.AbstractOptimizer
    env::Env
    inner::Problem
    objective_type::ObjectiveType
    annotations::Annotations
    #varmap::Dict{MOI.VariableIndex,VarId} # For the user to get VariablePrimal
    vars::CleverDicts.CleverDict{MOI.VariableIndex, Variable}
    #varids::CleverDicts.CleverDict{MOI.VariableIndex, VarId}
    moi_varids::Dict{VarId, MOI.VariableIndex}
    names_to_vars::Dict{String, MOI.VariableIndex}
    constrs::Dict{MOI.ConstraintIndex, Constraint}
    constrs_on_single_var_to_vars::Dict{MOI.ConstraintIndex, VarId}
    constrs_on_single_var_to_names::Dict{MOI.ConstraintIndex, String}
    names_to_constrs::Dict{String, MOI.ConstraintIndex}
    result::OptimizationState
    disagg_result::Union{Nothing, OptimizationState}
    default_optimizer_builder::Union{Nothing, Function}

    feasibility_sense::Bool # Coluna supports only Max or Min.

    function Optimizer()
        model = new()
        model.env = Env(Params())
        model.inner = Problem(model.env)
        model.annotations = Annotations()
        model.vars = CleverDicts.CleverDict{MOI.VariableIndex, Variable}()
        #model.varids = CleverDicts.CleverDict{MOI.VariableIndex, VarId}() # TODO : check if necessary to have two dicts for variables
        model.moi_varids = Dict{VarId, MOI.VariableIndex}()
        model.names_to_vars = Dict{String, MOI.VariableIndex}()
        model.constrs = Dict{MOI.ConstraintIndex, Union{Constraint, Nothing}}()
        model.constrs_on_single_var_to_vars = Dict{MOI.ConstraintIndex, VarId}()
        model.constrs_on_single_var_to_names = Dict{MOI.ConstraintIndex, String}()
        model.names_to_constrs = Dict{String, MOI.ConstraintIndex}()
        model.result = OptimizationState(get_optimization_target(model.inner))
        model.disagg_result = nothing
        model.default_optimizer_builder = nothing
        model.feasibility_sense = false
        return model
    end
end

MOI.Utilities.supports_default_copy_to(::Coluna.Optimizer, ::Bool) = true
MOI.supports(::Optimizer, ::MOI.VariableName, ::Type{MOI.VariableIndex}) = true
MOI.supports(::Optimizer, ::MOI.ConstraintName, ::Type{<:MOI.ConstraintIndex}) = true
MOI.supports_constraint(::Optimizer, ::Type{<:SupportedConstrFunc}, ::Type{<:SupportedConstrSets}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.SingleVariable}, ::Type{<: SupportedVarSets}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{<:SupportedObjFunc}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true
MOI.supports(::Optimizer, ::MOI.ConstraintPrimalStart) = false
MOI.supports(::Optimizer, ::MOI.ConstraintDualStart) = false

# Parameters
function MOI.set(model::Optimizer, param::MOI.RawParameter, val)
    if param.name == "params"
        model.env.params = val
    elseif param.name == "default_optimizer"
        optimizer_builder = () -> MoiOptimizer(MOI._instantiate_and_check(val))
        model.default_optimizer_builder = optimizer_builder
        set_default_optimizer_builder!(model.inner, optimizer_builder)
    else
        @warn("Unknown parameter $(param.name).")
    end
    return
end

function _get_orig_varid(optimizer::Optimizer, x::MOI.VariableIndex)
    if haskey(optimizer.vars, x)
        return optimizer.env.varids[x]
    end
    throw(MOI.InvalidIndex(x))
    return origid
end

function _get_orig_varid_in_form(
    optimizer::Optimizer, form::Formulation, x::MOI.VariableIndex
)
    origid = _get_orig_varid(optimizer, x)
    return getid(getvar(form, origid))
end

MOI.get(optimizer::Coluna.Optimizer, ::MOI.SolverName) = "Coluna"

function MOI.optimize!(optimizer::Optimizer)
    optimizer.result, optimizer.disagg_result = optimize!(
        optimizer.env, optimizer.inner, optimizer.annotations
    )
    return
end

function MOI.copy_to(dest::Coluna.Optimizer, src::MOI.ModelLike; kwargs...)
    return MOI.Utilities.automatic_copy_to(dest, src; kwargs...)
end

############################################################################################
# Add variables
############################################################################################
function MOI.add_variable(model::Coluna.Optimizer)
    orig_form = get_original_formulation(model.inner)
    var = setvar!(orig_form, "v", OriginalVar)
    index = CleverDicts.add_item(model.vars, var)
    model.moi_varids[getid(var)] = index
    index2 = CleverDicts.add_item(model.env.varids, getid(var))
    @assert index == index2
    return index
end

############################################################################################
# Add constraint
############################################################################################
function _constraint_on_variable!(var::Variable, ::MOI.Integer)
    # set perene data
    var.perendata.kind = Integ
    var.curdata.kind = Integ
    return
end

function _constraint_on_variable!(var::Variable, ::MOI.ZeroOne)
    # set perene data
    var.perendata.kind = Binary
    var.curdata.kind = Binary
    var.perendata.lb = max(0.0, var.perendata.lb)
    var.curdata.lb = max(0.0, var.curdata.lb)
    var.perendata.ub = min(1.0, var.perendata.ub)
    var.curdata.ub = min(1.0, var.curdata.ub)
    return
end

function _constraint_on_variable!(var::Variable, set::MOI.GreaterThan{Float64})
    # set perene data
    var.perendata.lb = max(set.lower, var.perendata.lb)
    var.curdata.lb = max(set.lower, var.perendata.lb)
    return
end

function _constraint_on_variable!(var::Variable, set::MOI.LessThan{Float64})
    # set perene data
    var.perendata.ub = min(set.upper, var.perendata.ub)
    var.curdata.ub = min(set.upper, var.curdata.ub)
    return
end

function _constraint_on_variable!(var::Variable, set::MOI.EqualTo{Float64})
    # set perene data
    var.perendata.lb = max(set.value, var.perendata.lb)
    var.curdata.lb = max(set.value, var.curdata.lb)
    var.perendata.ub = min(set.value, var.perendata.ub)
    var.curdata.ub = min(set.value, var.curdata.ub)
    return
end

function _constraint_on_variable!(var::Variable, set::MOI.Interval{Float64})
    # set perene data
    var.perendata.lb = max(set.lower, var.perendata.lb)
    var.curdata.lb = max(set.lower, var.curdata.lb)
    var.perendata.ub = min(set.upper, var.perendata.ub)
    var.curdata.ub = min(set.upper, var.curdata.ub)
    return
end

function MOI.add_constraint(
    model::Coluna.Optimizer, func::MOI.SingleVariable, set::S
) where {S<:SupportedVarSets}
    orig_form = get_original_formulation(model.inner)
    var = model.vars[func.variable]
    _constraint_on_variable!(var, set)
    constrid = MOI.ConstraintIndex{MOI.SingleVariable, S}(func.variable.value)
    model.constrs_on_single_var_to_names[constrid] = ""
    model.constrs_on_single_var_to_vars[constrid] = getid(var)
    return constrid
end

function MOI.add_constraint(
    model::Coluna.Optimizer, func::MOI.ScalarAffineFunction{Float64}, set::S
) where {S<:SupportedConstrSets}
    orig_form = get_original_formulation(model.inner)
    members = Dict{VarId, Float64}()
    for term in func.terms
        var = model.vars[term.variable_index]
        members[getid(var)] = term.coefficient
    end
    constr = setconstr!(
        orig_form, "c", OriginalConstr;
        rhs = MathProg.convert_moi_rhs_to_coluna(set),
        kind = Essential,
        sense = MathProg.convert_moi_sense_to_coluna(set),
        inc_val = 10.0,
        members = members
    )
    constrid =  MOI.ConstraintIndex{typeof(func), typeof(set)}(length(model.constrs))
    model.constrs[constrid] = constr
    return constrid
end

############################################################################################
# Get variables
############################################################################################
function MOI.get(model::Coluna.Optimizer, ::Type{MOI.VariableIndex}, name::String)
    return get(model.names_to_vars, name, nothing)
end

function MOI.get(model::Coluna.Optimizer, ::MOI.ListOfVariableIndices)
    indices = Vector{MathOptInterface.VariableIndex}()
    for (key,value) in model.moi_varids
        push!(indices, value)
    end
    return sort!(indices)
end

############################################################################################
# Get constraints
############################################################################################
function _moi_bounds_type(lb, ub)
    lb == ub && return MOI.EqualTo{Float64}
    lb == -Inf && ub < Inf && return MOI.LessThan{Float64}
    lb > -Inf && ub == Inf && return MOI.GreaterThan{Float64}
    lb > -Inf && ub < -Inf && return MOI.Interval{Float64}
    return nothing
end

function MOI.get(
    model::Coluna.Optimizer, C::Type{MOI.ConstraintIndex{F,S}}, name::String
) where {F,S}
    index = get(model.names_to_constrs, name, nothing)
    typeof(index) == C && return index
    return nothing
end

function MOI.get(model::Coluna.Optimizer, ::MOI.ListOfConstraints)
    orig_form = get_original_formulation(model.inner)
    constraints = Set{Tuple{DataType, DataType}}()
    for (id, var) in model.vars
        # Bounds
        lb = getperenlb(orig_form, var)
        ub = getperenub(orig_form, var)
        bound_type = _moi_bounds_type(lb, ub)
        if bound_type !== nothing
            push!(constraints, (MOI.SingleVariable, bound_type))
        end
        # Kind
        var_kind = MathProg.convert_coluna_kind_to_moi(getperenkind(orig_form, var))
        if var_kind !== nothing
            push!(constraints, (MOI.SingleVariable, var_kind))
        end
    end
    for (id, constr) in model.constrs
        constr_sense = MathProg.convert_coluna_sense_to_moi(getperensense(orig_form, constr))
        push!(constraints, (MOI.ScalarAffineFunction{Float64}, constr_sense))
    end
    return collect(constraints)
end

_add_constraint!(indices::Vector, index) = nothing
function _add_constraint!(
    indices::Vector{MOI.ConstraintIndex{F,S}}, index::MOI.ConstraintIndex{F,S}
) where {F,S}
    push!(indices, index)
    return
end

function MOI.get(
    model::Coluna.Optimizer, ::MOI.ListOfConstraintIndices{F, S}
) where {F<:MOI.ScalarAffineFunction{Float64}, S}
    indices = MOI.ConstraintIndex{F,S}[]
    for (id, constr) in model.constrs
        _add_constraint!(indices, id)
    end
    return sort!(indices, by = x -> x.value)
end

function MOI.get(
    model::Coluna.Optimizer, ::MOI.ListOfConstraintIndices{F, S}
) where {F<:MOI.SingleVariable, S}
    orig_form = get_original_formulation(model.inner)
    indices = MOI.ConstraintIndex{F,S}[]
    for (id, var) in model.vars
        if S == MathProg.convert_coluna_kind_to_moi(getperenkind(orig_form, var))
            push!(indices, MOI.ConstraintIndex{F,S}(id.value))
        end
    end
    return sort!(indices, by = x -> x.value)
end

function MOI.get(
    model::Coluna.Optimizer, ::MOI.ConstraintFunction, index::MOI.ConstraintIndex{F,S}
) where {F<:MOI.ScalarAffineFunction{Float64}, S}
    orig_form = get_original_formulation(model.inner)
    constrid = getid(model.constrs[index])
    terms = MOI.ScalarAffineTerm{Float64}[]
    for (varid, coeff) in @view getcoefmatrix(orig_form)[constrid, :]
        push!(terms, MOI.ScalarAffineTerm(coeff, model.moi_varids[varid]))
    end
    return MOI.ScalarAffineFunction(terms, 0.0)
end

function MOI.get(
    model::Coluna.Optimizer, ::MOI.ConstraintFunction, index::MOI.ConstraintIndex{F,S}
) where {F<:MOI.SingleVariable, S}
    return MOI.SingleVariable(MOI.VariableIndex(index.value))
end

function MOI.get(
    model::Coluna.Optimizer, ::MOI.ConstraintSet, index::MOI.ConstraintIndex{F,S}
) where {F<:MOI.ScalarAffineFunction{Float64},S}
    orig_form = get_original_formulation(model.inner)
    rhs = getperenrhs(orig_form, model.constrs[index])
    return S(rhs)
end

function MOI.get(
    model::Coluna.Optimizer, ::MOI.ConstraintSet,
    index::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}}
)
    orig_form = get_original_formulation(model.inner)
    lb = getperenlb(orig_form, model.vars[MOI.VariableIndex(index.value)])
    return MOI.GreaterThan(lb)
end

function MOI.get(
    model::Coluna.Optimizer, ::MOI.ConstraintSet,
    index::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}}
)
    orig_form = get_original_formulation(model.inner)
    ub = getperenub(orig_form, model.vars[MOI.VariableIndex(index.value)])
    return MOI.LessThan(ub)
end

function MOI.get(
    model::Coluna.Optimizer, ::MOI.ConstraintSet,
    index::MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}}
)
    orig_form = get_original_formulation(model.inner)
    lb = getperenlb(orig_form, model.vars[MOI.VariableIndex(index.value)])
    ub = getperenub(orig_form, model.vars[MOI.VariableIndex(index.value)])
    @assert lb == ub
    return MOI.EqualTo(lb)
end

function MOI.get(
    model::Coluna.Optimizer, ::MOI.ConstraintSet,
    index::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}
)
    orig_form = get_original_formulation(model.inner)
    lb = getperenlb(orig_form, model.vars[MOI.VariableIndex(index.value)])
    ub = getperenub(orig_form, model.vars[MOI.VariableIndex(index.value)])
    return MOI.Interval(lb, ub)
end

function MOI.get(
    model::Coluna.Optimizer, ::MOI.ConstraintSet,
    index::MOI.ConstraintIndex{MOI.SingleVariable, MOI.ZeroOne}
)
    return MOI.ZeroOne()
end

function MOI.get(
    model::Coluna.Optimizer, ::MOI.ConstraintSet,
    index::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Integer}
)
    return MOI.Integer()
end

function MOI.get(model::Coluna.Optimizer, ::Type{MOI.ConstraintIndex}, name::String)
    return get(model.names_to_constrs, name, nothing)
end

############################################################################################
# Attributes of variables
############################################################################################
function MOI.set(
    model::Coluna.Optimizer, ::BD.VariableDecomposition, varid::MOI.VariableIndex,
    annotation::BD.Annotation
)
    store!(model.annotations, annotation, model.vars[varid])
    return
end

function MOI.set(
    model::Coluna.Optimizer, ::MOI.VariableName, varid::MOI.VariableIndex, name::String
)
    var = model.vars[varid]
    # TODO : rm set perene name
    var.name = name
    model.names_to_vars[name] = varid
    return
end

function MOI.set(
    model::Coluna.Optimizer, ::BD.VarBranchingPriority, varid::MOI.VariableIndex, branching_priority::Int
)
    var = model.vars[varid]
    var.branching_priority = Float64(branching_priority)
    return
end

function MOI.get(model::Coluna.Optimizer, ::MOI.VariableName, index::MOI.VariableIndex)
    orig_form = get_original_formulation(model.inner)
    return getname(orig_form, model.vars[index])
end

function MOI.get(model::Coluna.Optimizer, ::BD.VarBranchingPriority, varid::MOI.VariableIndex)
    var = model.vars[varid]
    return var.branching_priority
end

function MOI.get(model::Optimizer, ::MOI.ListOfVariableAttributesSet)
    return MOI.AbstractVariableAttribute[MOI.VariableName()]
end

############################################################################################
# Attributes of constraints
############################################################################################
function MOI.set(
    model::Coluna.Optimizer, ::BD.ConstraintDecomposition, constrid::MOI.ConstraintIndex,
    annotation::BD.Annotation
)
    constr = get(model.constrs, constrid, nothing)
    if constr !== nothing
        store!(model.annotations, annotation, constr)
    end
    return
end

function MOI.set(
    model::Coluna.Optimizer, ::MOI.ConstraintName, constrid::MOI.ConstraintIndex{F,S}, name::String
) where {F<:MOI.ScalarAffineFunction,S}
    MOI.throw_if_not_valid(model, constrid)
    constr = model.constrs[constrid]
    # TODO : rm set perene name
    constr.name = name
    model.names_to_constrs[name] = constrid
    return
end

function MOI.set(
    model::Coluna.Optimizer, ::MOI.ConstraintName, constrid::MOI.ConstraintIndex{F,S}, name::String
) where {F<:MOI.SingleVariable,S}
    MOI.throw_if_not_valid(model, constrid)
    model.constrs_on_single_var_to_names[constrid] = name
    model.names_to_constrs[name] = constrid
    return
end

function MOI.get(model::Coluna.Optimizer, ::MOI.ConstraintName, constrid::MOI.ConstraintIndex)
    MOI.throw_if_not_valid(model, constrid)
    orig_form = get_original_formulation(model.inner)
    constr = get(model.constrs, constrid, nothing)
    if constr !== nothing
        return getname(orig_form, constr)
    end
    return ""
end

function MOI.get(model::Optimizer, ::MOI.ListOfConstraintAttributesSet)
    return MOI.AbstractConstraintAttribute[MOI.ConstraintName()]
end

############################################################################################
# Objective
############################################################################################
function MOI.set(model::Coluna.Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    orig_form = get_original_formulation(model.inner)
    if sense == MOI.MIN_SENSE
        model.feasibility_sense = false
        set_objective_sense!(orig_form, true) # Min
    elseif sense == MOI.MAX_SENSE
        model.feasibility_sense = false
        set_objective_sense!(orig_form, false) # Max
    else
        model.feasibility_sense = true
        set_objective_sense!(orig_form, true) # Min
    end
    return
end

function MOI.get(model::Coluna.Optimizer, ::MOI.ObjectiveSense)
    sense = getobjsense(get_original_formulation(model.inner))
    model.feasibility_sense && return MOI.FEASIBILITY_SENSE
    sense == MaxSense && return MOI.MAX_SENSE
    return MOI.MIN_SENSE
end

function MOI.get(model::Coluna.Optimizer, ::MOI.ObjectiveFunctionType)
    if model.objective_type == SINGLE_VARIABLE
        return MOI.SingleVariable
    end
    @assert model.objective_type == SCALAR_AFFINE
    return MOI.ScalarAffineFunction{Float64}
end

function MOI.set(
    model::Coluna.Optimizer, ::MOI.ObjectiveFunction{F}, func::F
) where {F<:MOI.ScalarAffineFunction{Float64}}
    model.objective_type = SCALAR_AFFINE
    for term in func.terms
        var = model.vars[term.variable_index]
        cost = term.coefficient
        # TODO : rm set peren cost
        var.perendata.cost = cost
        var.curdata.cost = cost
    end
    if func.constant != 0
        orig_form = get_original_formulation(model.inner)
        setobjconst!(orig_form, func.constant)
    end
    return
end

function MOI.set(
    model::Coluna.Optimizer, ::MOI.ObjectiveFunction{MOI.SingleVariable},
    func::MOI.SingleVariable
)
    model.objective_type = SINGLE_VARIABLE
    var = model.vars[func.variable]
    # TODO : rm set perene cost
    var.perendata.cost = 1.0
    var.curdata.cost = 1.0
    return
end

function MOI.get(
    model::Coluna.Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}
)
    @assert model.objective_type == SCALAR_AFFINE
    orig_form = get_original_formulation(model.inner)
    terms = MOI.ScalarAffineTerm{Float64}[]
    for (id, var) in model.vars
        cost = getperencost(orig_form, var)
        iszero(cost) && continue
        push!(terms, MOI.ScalarAffineTerm(cost, id))
    end
    return MOI.ScalarAffineFunction(terms, getobjconst(orig_form))
end

function MOI.get(
    model::Coluna.Optimizer, ::MOI.ObjectiveFunction{MOI.SingleVariable}
)
    @assert model.objective_type == SINGLE_VARIABLE
    orig_form = get_original_formulation(model.inner)
    for (id, var) in model.vars
        cost = getperencost(orig_form, var)
        if cost != 0
            return MOI.SingleVariable(id)
        end
    end
    error("Could not find the variable with cost != 0.")
end

############################################################################################
# Attributes of model
############################################################################################
function MOI.set(model::Coluna.Optimizer, ::BD.DecompositionTree, tree::BD.Tree)
    model.annotations.tree = tree
    return
end

function MOI.set(model::Coluna.Optimizer, ::BD.ObjectiveDualBound, db)
    set_initial_dual_bound!(model.inner, db)
    return
end

function MOI.set(model::Coluna.Optimizer, ::BD.ObjectivePrimalBound, pb)
    set_initial_primal_bound!(model.inner, pb)
    return
end

function MOI.empty!(model::Coluna.Optimizer)
    model.inner = Problem(model.env)
    model.annotations = Annotations()
    model.vars = CleverDicts.CleverDict{MOI.VariableIndex, Variable}()
    model.env.varids = CleverDicts.CleverDict{MOI.VariableIndex, VarId}()
    model.moi_varids = Dict{VarId, MOI.VariableIndex}()
    model.constrs = Dict{MOI.ConstraintIndex, Constraint}()
    if model.default_optimizer_builder !== nothing
        set_default_optimizer_builder!(model.inner, model.default_optimizer_builder)
    end
    model.result = OptimizationState(get_optimization_target(model.inner))
    model.disagg_result = nothing
    return
end

mutable struct ColumnInfo <: BD.AbstractColumnInfo
    optimizer::Coluna.Optimizer
    column_var_id::VarId
    column_val::Float64
end

function BD.getsolutions(model::Coluna.Optimizer, k)
    ip_primal_sol = model.disagg_result.ip_primal_sols[k]
    sp_columns_info = Vector{ColumnInfo}()
    for (varid, val) in ip_primal_sol
        if getduty(varid) <= MasterCol
            push!(sp_columns_info, ColumnInfo(model, varid, val))
        end
    end
    return sp_columns_info
end

BD.value(info::ColumnInfo) = info.column_val

function BD.value(info::ColumnInfo, index::MOI.VariableIndex)
    varid = info.optimizer.env.varids[index]
    origin_form_uid = getoriginformuid(info.column_var_id)
    spform = get_dw_pricing_sps(info.optimizer.inner.re_formulation)[origin_form_uid]
    return info.column_val * getprimalsolmatrix(spform)[varid, info.column_var_id]
end

function MOI.get(model::Coluna.Optimizer, ::MOI.NumberOfVariables)
    orig_form = get_original_formulation(model.inner)
    return length(getvars(orig_form))
end

function MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{F, S}) where {F, S}
    return length(MOI.get(model, MOI.ListOfConstraintIndices{F, S}()))
end

function MOI.get(model::Optimizer, ::MOI.ListOfModelAttributesSet)
    attributes = Any[MOI.ObjectiveSense()]
    typ = MOI.get(model, MOI.ObjectiveFunctionType())
    if typ !== nothing
        push!(attributes, MOI.ObjectiveFunction{typ}())
    end
    return attributes
end

# ######################
# ### Get functions ####
# ######################

function MOI.is_empty(optimizer::Optimizer)
    return optimizer.inner === nothing || optimizer.inner.re_formulation === nothing
end

function MOI.is_valid(
    optimizer::Optimizer, index::MOI.ConstraintIndex{F,S}
) where {F<:MOI.SingleVariable,S}
    return haskey(optimizer.constrs_on_single_var_to_names, index)
end

function MOI.is_valid(
    optimizer::Optimizer, index::MOI.ConstraintIndex{F,S}
) where {F<:MOI.ScalarAffineFunction,S}
    return haskey(optimizer.constrs, index)
end

function MOI.is_valid(optimizer::Optimizer, index::MOI.VariableIndex)
    return haskey(optimizer.vars, index)
end

function MOI.get(optimizer::Optimizer, ::MOI.ObjectiveBound)
    return getvalue(get_ip_dual_bound(optimizer.result))
end

function MOI.get(optimizer::Optimizer, ::MOI.ObjectiveValue)
    return getvalue(get_ip_primal_bound(optimizer.result))
end

function MOI.get(optimizer::Optimizer, ::MOI.RelativeGap)
    return ip_gap(optimizer.result)
end

function MOI.get(optimizer::Optimizer, ::MOI.VariablePrimal, ref::MOI.VariableIndex)
    id = getid(optimizer.vars[ref]) # This gets a coluna Id{Variable}
    best_primal_sol = get_best_ip_primal_sol(optimizer.result)
    if best_primal_sol === nothing
        @warn "Coluna did not find a primal feasible solution."
        return NaN
    end
    return get(best_primal_sol, id, 0.0)
end

function MOI.get(optimizer::Optimizer, ::MOI.VariablePrimal, refs::Vector{MOI.VariableIndex})
    best_primal_sol = get_best_ip_primal_sol(optimizer.result)
    if best_primal_sol === nothing
        @warn "Coluna did not find a primal feasible solution."
        return [NaN for ref in refs]
    end
    return [get(best_primal_sol, getid(optimizer.vars[ref]), 0.0) for ref in refs]
end

function MOI.get(optimizer::Optimizer, ::MOI.TerminationStatus)
    return convert_status(getterminationstatus(optimizer.result))
end

function MOI.get(optimizer::Optimizer, ::MOI.PrimalStatus)
    primal_sol = get_best_ip_primal_sol(optimizer.result)
    primal_sol === nothing && return MOI.NO_SOLUTION
    return convert_status(getstatus(primal_sol))
end

function MOI.get(optimizer::Optimizer, ::MOI.DualStatus)    
    dual_sol = get_best_lp_dual_sol(optimizer.result)
    dual_sol === nothing && return MOI.NO_SOLUTION
    return convert_status(getstatus(dual_sol))
end

function MOI.get(optimizer::Optimizer, ::MOI.RawStatusString)
    return string(getterminationstatus(optimizer.result))
end

function MOI.get(optimizer::Optimizer, ::MOI.ResultCount)
    return length(get_ip_primal_sols(optimizer.result))
end

function MOI.get(
    optimizer::Optimizer, ::MOI.ConstraintPrimal, index::MOI.ConstraintIndex{F,S}
) where {F<:MOI.SingleVariable,S}
    varid = get(optimizer.constrs_on_single_var_to_vars, index, nothing)
    if varid === nothing
        @warn "Could not find constraint with id $(index)."
        return NaN
    end
    best_primal_sol = get_best_ip_primal_sol(optimizer.result)
    return get(best_primal_sol, varid, 0.0)
end

function MOI.get(optimizer::Optimizer, ::MOI.ConstraintPrimal, index::MOI.ConstraintIndex)
    constrid = get(optimizer.constrs, index, nothing)
    if constrid === nothing
        @warn "Could not find constraint with id $(index)."
        return NaN
    end
    best_primal_sol = get_best_ip_primal_sol(optimizer.result)
    return constraint_primal(best_primal_sol, getid(constrid))
end

MOI.get(optimizer::Optimizer, ::MOI.NodeCount) = optimizer.env.kpis.node_count
MOI.get(optimizer::Optimizer, ::MOI.SolveTime) = optimizer.env.kpis.elapsed_optimization_time

# function MOI.get(optimizer::Optimizer, ::MOI.ConstraintDual, index::MOI.ConstraintIndex)
#     return 0.0
# end

# function MOI.get(optimizer::Optimizer, ::MOI.SolveTime)
#     return 0.0
# end