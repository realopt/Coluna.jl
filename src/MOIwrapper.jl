const CleverDicts = MOI.Utilities.CleverDicts

const SupportedObjFunc = Union{MOI.ScalarAffineFunction{Float64},
    MOI.SingleVariable}

const SupportedVarSets = Union{MOI.ZeroOne, MOI.Integer, MOI.LessThan{Float64},
    MOI.EqualTo{Float64}, MOI.GreaterThan{Float64}}

const SupportedConstrFunc = Union{MOI.ScalarAffineFunction{Float64}}

const SupportedConstrSets = Union{MOI.EqualTo{Float64}, MOI.GreaterThan{Float64},
    MOI.LessThan{Float64}}

@enum(
    ObjectiveType,
    SINGLE_VARIABLE,
    SCALAR_AFFINE
)
mutable struct Optimizer <: MOI.AbstractOptimizer
    inner::Problem
    objective_type::ObjectiveType
    params::Params
    annotations::Annotations
    #varmap::Dict{MOI.VariableIndex,VarId} # For the user to get VariablePrimal
    vars::CleverDicts.CleverDict{MOI.VariableIndex, Variable}
    varids::CleverDicts.CleverDict{MOI.VariableIndex, VarId}
    moi_varids::Dict{VarId, MOI.VariableIndex}
    constrs::Dict{MOI.ConstraintIndex, Constraint}
    result::Union{Nothing,OptimizationState}

    function Optimizer()
        model = new()
        model.inner = Problem()
        model.params = Params()
        model.annotations = Annotations()
        model.vars = CleverDicts.CleverDict{MOI.VariableIndex, Variable}()
        model.varids = CleverDicts.CleverDict{MOI.VariableIndex, VarId}() # TODO : check if necessary to have two dicts for variables
        model.moi_varids = Dict{VarId, MOI.VariableIndex}()
        model.constrs = Dict{MOI.ConstraintIndex, Constraint}()
        return model
    end
end

MOI.Utilities.supports_default_copy_to(::Coluna.Optimizer, ::Bool) = true
MOI.supports(::Optimizer, ::MOI.VariableName, ::Type{MOI.VariableIndex}) = true
MOI.supports(::Optimizer, ::MOI.ConstraintName, ::Type{<:MOI.ConstraintIndex}) = true
MOI.supports_constraint(::Optimizer, ::Type{<:SupportedConstrFunc}, ::Type{<:SupportedConstrSets}) = true
MOI.supports_constraint(::Optimizer, ::Type{MOI.SingleVariable}, ::Type{<: SupportedVarSets}) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{<:SupportedObjFunc}) = true

# Parameters
function MOI.set(model::Optimizer, param::MOI.RawParameter, val)
    if param.name == "params"
        model.params = val
    elseif param.name == "default_optimizer"
        optimizer_builder = () -> MoiOptimizer(val())
        model.inner.default_optimizer_builder = optimizer_builder
    else
        @warn("Unknown parameter $(param.name).")
    end
    return
end

function _get_orig_varid(optimizer::Optimizer, x::MOI.VariableIndex)
    if haskey(optimizer.vars, x)
        return optimizer.varids[x]
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
    optimizer.result = optimize!(
        optimizer.inner, optimizer.annotations, optimizer.params
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
    index2 = CleverDicts.add_item(model.varids, getid(var))
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
    var.perendata.lb = 0.0
    var.curdata.lb = 0.0
    var.perendata.ub = 1.0
    var.curdata.ub = 1.0
    return
end

function _constraint_on_variable!(var::Variable, set::MOI.GreaterThan{Float64})
    # set perene data
    var.perendata.lb = set.lower
    var.curdata.lb = set.lower
    return
end

function _constraint_on_variable!(var::Variable, set::MOI.LessThan{Float64})
    # set perene data
    var.perendata.ub = set.upper
    var.curdata.ub = set.upper
    return
end

function _constraint_on_variable!(var::Variable, set::MOI.EqualTo{Float64})
    # set perene data
    var.perendata.lb = set.value
    var.curdata.lb = set.value
    var.perendata.ub = set.value
    var.curdata.ub = set.value
    return
end

function _constraint_on_variable!(var::Variable, set::MOI.Interval{Float64})
    # set perene data
    var.perendata.lb = set.lower
    var.curdata.lb = set.lower
    var.perendata.ub = set.upper
    var.curdata.ub = set.upper
    return
end

function MOI.add_constraint(
    model::Coluna.Optimizer, func::MOI.SingleVariable, set::S
) where {S<:SupportedVarSets}
    orig_form = get_original_formulation(model.inner)
    var = model.vars[func.variable]
    _constraint_on_variable!(var, set)
    return MOI.ConstraintIndex{MOI.SingleVariable, S}(func.variable.value)
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
# Get constraints
############################################################################################
function _moi_bounds_type(lb, ub)
    lb == ub && return MOI.EqualTo{Float64}
    lb == -Inf && ub < Inf && return MOI.LessThan{Float64}
    lb > -Inf && ub == Inf && return MOI.GreaterThan{Float64}
    lb > -Inf && ub < -Inf && return MOI.Interval{Float64}
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
    return
end

function MOI.get(model::Coluna.Optimizer, ::MOI.VariableName, index::MOI.VariableIndex)
    orig_form = get_original_formulation(model.inner)
    return getname(orig_form, model.vars[index])
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
        store!(model.annotations, annotation, model.constrs[constrid])
    end
    return
end

function MOI.set(
    model::Coluna.Optimizer, ::MOI.ConstraintName, constrid::MOI.ConstraintIndex, name::String
)
    constr = model.constrs[constrid]
    # TODO : rm set perene name
    constr.name = name
    return
end

function MOI.get(model::Coluna.Optimizer, ::MOI.ConstraintName, index::MOI.ConstraintIndex)
    orig_form = get_original_formulation(model.inner)
    constr = get(model.constrs, index, nothing)
    if constr !== nothing
        return getname(orig_form, constr)
    end
    return ""
end

############################################################################################
# Objective
############################################################################################
function MOI.set(model::Coluna.Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    orig_form = get_original_formulation(model.inner)
    min_sense = (sense == MOI.MIN_SENSE)
    set_objective_sense!(orig_form, min_sense)
    return
end

function MOI.get(model::Coluna.Optimizer, ::MOI.ObjectiveSense)
    sense = getobjsense(get_original_formulation(model.inner))
    sense == MinSense && return MOI.MIN_SENSE
    return MOI.MAX_SENSE
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
    # TODO : missing constant
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
    # TODO : missing constant
    return MOI.ScalarAffineFunction(terms, 0.0)
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

function MOI.empty!(optimizer::Optimizer)
    optimizer.inner.re_formulation = nothing
end

function MOI.get(model::Coluna.Optimizer, ::MOI.NumberOfVariables)
    orig_form = get_original_formulation(model.inner)
    return length(getvars(orig_form))
end

# ######################
# ### Get functions ####
# ######################

function MOI.is_empty(optimizer::Optimizer)
    return optimizer.inner === nothing || optimizer.inner.re_formulation === nothing
end

function MOI.get(optimizer::Optimizer, ::MOI.ObjectiveBound)
    return getvalue(get_ip_dual_bound(optimizer.result))
end

function MOI.get(optimizer::Optimizer, ::MOI.ObjectiveValue)
    return getvalue(get_ip_primal_bound(optimizer.result))
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

function MOI.get(optimizer::Optimizer, object::MOI.TerminationStatus)
    return MathProg.convert_status(getterminationstatus(optimizer.result))
end
