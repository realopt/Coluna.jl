getid(vc::AbstractVarConstr) = vc.id
getmoirecord(vc::AbstractVarConstr) = vc.moirecord
getoriginformuid(vc::AbstractVarConstr) = getoriginformuid(getid(vc))

# Variables
## Cost
"""
    getperencost(formulation, variable)
    getperencost(formulation, varid)

Return the cost as defined by the user of a variable in a formulation.
"""
getperencost(form::Formulation, varid::VarId) = getperencost(form, getvar(form, varid))
getperencost(form::Formulation, var::Variable) = var.perendata.cost

"""
    getcurcost(formulation, variable)
    getcurcost(formulation, varid)

Return the current cost of the variable in the formulation.
"""
getcurcost(form::Formulation, varid::VarId) = getcurcost(form, getvar(form, varid))
getcurcost(form::Formulation, var::Variable) = var.curdata.cost

"""
    setcurcost!(formulation, varid, cost::Float64)
    setcurcost!(formulation, variable, cost::Float64)

Set the current cost of variable in the formulation.
If the variable is active and explicit, this change is buffered before application to the 
subsolver.
"""
function setcurcost!(form::Formulation, var::Variable, cost::Float64)
    var.curdata.cost = cost
    if isexplicit(form, var) && iscuractive(form, var)
        change_cost!(form.buffer, getid(var))
    end
    return
end

function setcurcost!(form::Formulation, varid::VarId, cost::Float64)
    return setcurcost!(form, getvar(form, varid), cost)
end

## Lower bound
"""
    getperenlb(formulation, varid)
    getperenlb(formulation, var)

Return the lower bound as defined by the user of a variable in a formulation.
"""
getperenlb(form::Formulation, varid::VarId) = getperenlb(form, getvar(form, varid))
getperenlb(form::Formulation, var::Variable) = var.perendata.lb

"""
    getcurlb(formulation, varid)
    getcurlb(formulation, var)

Return the current lower bound of a variable in a formulation.
"""
getcurlb(form::Formulation, varid::VarId) = getcurlb(form, getvar(form, varid))
getcurlb(form::Formulation, var::Variable) = var.curdata.lb

"""
    setcurlb!(formulation, varid, lb::Float64)
    setcurlb!(formulation, var, lb::Float64)

Set the current lower bound of a variable in a formulation.
If the variable is active and explicit, change is buffered before application to the
subsolver.
"""
function setcurlb!(form::Formulation, var::Variable, lb::Float64)
    var.curdata.lb = lb
    if isexplicit(form, var) && iscuractive(form, var)
        change_bound!(form.buffer, getid(var))
    end
    return
end
setcurlb!(form::Formulation, varid::VarId, lb::Float64) =  setcurlb!(form, getvar(form, varid), lb)

## Upper bound
"""
    getperenub(formulation, varid)
    getperenub(formulation, var)

Return the upper bound as defined by the user of a variable in a formulation.
"""
getperenub(form::Formulation, varid::VarId) = getperenub(form, getvar(form, varid))
getperenub(form::Formulation, var::Variable) = var.perendata.ub

"""
    getcurub(formulation, varid)
    getcurub(formulation, var)

Return the current upper bound of a variable in a formulation.
"""
getcurub(form::Formulation, varid::VarId) = getcurub(form, getvar(form, varid))
getcurub(form::Formulation, var::Variable) = var.curdata.ub

"""
    setcurub!(formulation, varid, ub::Float64)
    setcurub!(formulation, var, ub::Float64)

Set the current upper bound of a variable in a formulation.
If the variable is active and explicit, change is buffered before application to the
subsolver.
"""
function setcurub!(form::Formulation, var::Variable, ub::Float64)
    var.curdata.ub = ub
    if isexplicit(form, var) && iscuractive(form, var)
        change_bound!(form.buffer, getid(var))
    end
    return
end
setcurub!(form::Formulation, varid::VarId, ub::Float64) = setcurub!(form, getvar(form, varid), ub)


# Constraint
## rhs
"""
    getperenrhs(formulation, constraint)
    getperenrhs(formulation, constrid)

Return the right-hand side as defined by the user of a constraint in a formulation.
"""
getperenrhs(form::Formulation, constrid::ConstrId) = getperenrhs(form, getconstr(form, constrid))
getperenrhs(form::Formulation, constr::Constraint) = constr.perendata.rhs

"""
    getcurrhs(formulation, constraint)
    getcurrhs(formulation, constrid)

Return the current right-hand side of a constraint in a formulation.
"""
getcurrhs(form::Formulation, constrid::ConstrId) = getcurrhs(form, getconstr(form, constrid))
getcurrhs(form::Formulation, constr::Constraint) = constr.curdata.rhs

"""
    setcurrhs(formulation, constraint, rhs::Float64)
    setcurrhs(formulation, constrid, rhs::Float64)

Set the current right-hand side of a constraint in a formulation. 
If the constraint is active and explicit, this change is buffered before application to the
subsolver.
"""
function setcurrhs!(form::Formulation, constr::Constraint, rhs::Float64) 
    constr.curdata.rhs = rhs
    if isexplicit(form, constr) && iscuractive(form, constr)
        change_rhs!(form.buffer, getid(constr))
    end
    return
end
setcurrhs!(form::Formulation, constrid::ConstrId, rhs::Float64) = setcurrhs!(form, getconstr(form, constrid), rhs)


# Constraint
## incval
"""
    getcurincval(formulation, constraint)
    getcurincval(formulation, constrid)

Return the current incumbent value of a constraint in a formulation.
"""
getcurincval(form::Formulation, constrid::ConstrId) = getcurincval(form, getconstr(form, constrid))
getcurincval(form::Formulation, constr::Constraint) = constr.curdata.inc_val

"""
    setcurincval!(formulation, constraint, inc_val::Float64)
    setcurincval!(formulation, constrid, inc_val::Float64)

Set the current incumbent value of a constraint in a formulation. 
"""
function setcurincval!(form::Formulation, constr::Constraint, inc_val::Float64) 
    constr.curdata.inc_val = inc_val
    return
end
setcurincval!(form::Formulation, constrid::ConstrId, inc_val::Float64) =
    setcurincval!(form, getconstr(form, constrid), inc_val)

# Variable & Constraints
## kind
"""
    getperenkind(formulation, varconstr)
    getperenkind(formulation, varconstrid)

Return the kind as defined by the user of a variable or a constraint in a formulation.

Kinds of variable (`enum VarKind`) are `Continuous`, `Binary`, or `Integ`.
    
Kinds of a constraint (`enum ConstrKind`) are : 
 - `Essential` when the constraint structures the problem
 - `Facultative` when the constraint does not structure the problem
 - `SubSystem` (to do)

The kind of a constraint cannot change.
"""
getperenkind(form::Formulation, varid::VarId) = getperenkind(form, getvar(form, varid))
getperenkind(form::Formulation, var::Variable) = var.perendata.kind
getperenkind(form::Formulation, constrid::ConstrId) = getperenkind(form, getconstr(form, constrid))
getperenkind(form::Formulation, constr::Constraint) = constr.perendata.kind

"""
    getcurkind(formulation, variable)
    getcurkind(formulation, varid)

Return the current kind of a variable in a formulation.
"""
getcurkind(form::Formulation, varid::VarId) = getcurkind(form, getvar(form, varid))
getcurkind(form::Formulation, var::Variable) = var.curdata.kind

"""
    setcurkind!(formulation, variable, kind::VarKind)
    setcurkind!(formulation, varid, kind::VarKind)

Set the current kind of a variable in a formulation.
If the variable is active and explicit, this change is buffered before
application to the subsolver
"""
function setcurkind!(form::Formulation, var::Variable, kind::VarKind)
    var.curdata.kind = kind
    if isexplicit(form, var) && iscuractive(form, var)
        change_kind!(form.buffer, getid(var))
    end
    return
end
setcurkind!(form::Formulation, varid::VarId, kind::VarKind) = setcurkind!(form, getvar(form, varid), kind) 


## sense
function _senseofvar(lb::Float64, ub::Float64)
    lb >= 0 && return Positive
    ub <= 0 && return Negative
    return Free 
end

"""
    getperensense(formulation, varconstr)
    getperensense(formulation, varconstrid)

Return the sense as defined by the user of a variable or a constraint in a formulation.

Senses or a variable are (`enum VarSense`)  `Positive`, `Negative`, and `Free`.
Senses or a constraint are (`enum ConstrSense`) `Greater`, `Less`, and `Equal`.

The perennial sense of a variable depends on its perennial bounds.
"""
getperensense(form::Formulation, varid::VarId) = getperensense(form, getvar(form, varid))
getperensense(form::Formulation, var::Variable) = _senseofvar(getperenlb(form, var), getperenub(form, var))
getperensense(form::Formulation, constrid::ConstrId) = getperensense(form, getconstr(form, constrid))
getperensense(form::Formulation, constr::Constraint) = constr.perendata.sense

"""
    getcursense(formulation, varconstr)
    getcursense(formulation, varconstrid)

Return the current sense of a variable or a constraint in a formulation.
The current sense of a variable depends on its current bounds.
"""
getcursense(form::Formulation, varid::VarId) = getcursense(form, getconstr(form, varid))
getcursense(form::Formulation, var::Variable) = _senseofvar(getcurlb(form, var), getcurub(form, var))
getcursense(form::Formulation, constrid::ConstrId) = getcursense(form, getconstr(form, constrid))
getcursense(form::Formulation, constr::Constraint) = constr.curdata.sense

"""
    setcursense!(formulation, constr, sense::ConstrSense)
    setcursense!(formulation, constrid, sense::ConstrSense)

Set the current sense of a constraint in a formulation.

This method is not applicable to variables because the sense of a variable depends on its
bounds.
"""
function setcursense!(form::Formulation, constr::Constraint, sense::ConstrSense)
    constr.curdata.sense = sense
    #if isexplicit(form, constr) 
        #change_sense!(form.buffer, getvar(form, varid))
    #end
    return
end

function setcursense!(form::Formulation, constrid::ConstrId, sense::ConstrSense)
    return setcursense!(form, getconstr(form, constrid), sense)
end

## inc_val
"""
    getperenincval(formulation, varconstrid)
    getperenincval(formulation, varconstr)

Return the incumbent value as defined by the user of a variable or a constraint in a formulation. 
The incumbent value is ?
"""
getperenincval(form::Formulation, varid::VarId) = getperenincval(form, getvar(form, varid))
getperenincval(form::Formulation, var::Variable) = var.perendata.inc_val
getperenincval(form::Formulation, constrid::ConstrId) = getperenincval(form, getconstr(form, constrid))
getperenincval(form::Formulation, constr::Constraint) = constr.perendata.inc_val

"""
    getcurincval(formulation, varconstrid)
    getcurincval(formulation, varconstr)

Return the current incumbent value of a variable or a constraint in a formulation.
"""
getcurincval(form::Formulation, varid::VarId) = getcurincval(form, getvar(form, varid))
getcurincval(form::Formulation, var::Variable) = var.curdata.inc_val
getcurincval(form::Formulation, constrid::ConstrId) = getcurincval(form, getconstr(form, constrid))
getcurincval(form::Formulation, constr::Constraint) = constr.curdata.inc_val

"""
    setcurincval!(formulation, varconstrid, value::Real)

Set the current incumbent value of a variable or a constraint in a formulation.
"""
function setcurincval!(form::Formulation, var::Variable, inc_val::Real)
    var.curdata.inc_val = inc_val
    return
end
setcurincval!(form::Formulation, varid::VarId, inc_val::Real) = setcurincval!(form, getvar(form, varid), inc_val)
function setcurincval!(form::Formulation, constr::Constraint, inc_val::Real)
    constr.curdata.inc_val = inc_val
    return
end
setcurincval!(form::Formulation, constrid::ConstrId, inc_val::Real) = setcurincval!(form, getconstr(form, constrid), inc_val)

## active
"""
    isperenactive(formulation, varconstrid)
    isperenactive(formulation, varconstr)

Return `true` if the variable or the constraint is active in the formulation; `false` otherwise.
A variable (or a constraint) is active if it is used in the formulation. You can fake the 
deletion of the variable by deativate it. This allows you to keep the variable if you want 
to reactivate it later.
"""
isperenactive(form::Formulation, varid::VarId) = isperenactive(form, getvar(form, varid))
isperenactive(form::Formulation, var::Variable) = var.perendata.is_active
isperenactive(form::Formulation, constrid::ConstrId) = isperenactive(form, getconstr(form, constrid))
isperenactive(form::Formulation, constr::Constraint) = constr.perendata.is_active

"""
    iscuractive(formulation, varconstrid)
    iscuractive(formulation, varconstr)

Return `true` if the variable or the constraint is currently active; `false` otherwise.
"""
iscuractive(form::Formulation, varid::VarId) = iscuractive(form, getvar(form, varid))
iscuractive(form::Formulation, var::Variable) = var.curdata.is_active
iscuractive(form::Formulation, constrid::ConstrId) = iscuractive(form, getconstr(form, constrid))
iscuractive(form::Formulation, constr::Constraint) = constr.curdata.is_active

function _activate!(form::Formulation, varconstr::AbstractVarConstr)
    if isexplicit(form, varconstr) && !iscuractive(form, varconstr)
        add!(form.buffer, getid(varconstr))
    end
    varconstr.curdata.is_active = true
    return
end

"""
    activate!(formulation, varconstrid)
    activate!(formulation, varconstr)

Activate a variable or a constraint in a formulation.

    activate!(formulation, function)

It is also possible to activate variables and constraints of a formulation such that 
`function(varconstrid)` returns `true`.
"""
function activate!(form::Formulation, constr::Constraint)
    _activate!(form, constr)
    for varid in constr.art_var_ids
        _activate!(form, getvar(form, varid))
    end
    return
end

activate!(form::Formulation, var::Variable) = _activate!(form, var)

function activate!(form::Formulation, varconstrid::Id{VC}) where {VC <: AbstractVarConstr}
    return activate!(form, getelem(form, varconstrid))
end

function activate!(form::Formulation, f::Function)
    for (varid, var) in getvars(form)
        if !iscuractive(form, varid) && f(varid)
            activate!(form, var)
        end
    end
    for (constrid, constr) in getconstrs(form)
        if !iscuractive(form, constrid) && f(constrid)
            activate!(form, constr)
        end
    end
    return
end

function _deactivate!(form::Formulation, varconstr::AbstractVarConstr)
    if isexplicit(form, varconstr) && iscuractive(form, varconstr)
        remove!(form.buffer, getid(varconstr))
    end
    varconstr.curdata.is_active = false
    return
end

"""
    deactivate!(formulation, varconstrid)
    deactivate!(formulation, varconstr)

Deactivate a variable or a constraint in a formulation.

    deactivate!(formulation, function)

It is also possible to deactivate variables and constraints such that 
`function(varconstrid)` returns `true`.
"""
function deactivate!(form::Formulation, constr::Constraint)
    _deactivate!(form, constr)
    for varid in constr.art_var_ids
        _deactivate!(form, getvar(form, varid))
    end
    return
end
deactivate!(form::Formulation, var::Variable) = _deactivate!(form, var)
deactivate!(form::Formulation, varconstrid::Id{VC}) where {VC<:AbstractVarConstr} = deactivate!(form, getelem(form, varconstrid))

function deactivate!(form::Formulation, f::Function)
    for (varid, var) in getvars(form)
        if iscuractive(form, var) && f(varid)
            deactivate!(form, var)
        end
    end
    for (constrid, constr) in getconstrs(form)
        if iscuractive(form, constr) && f(constrid)
            deactivate!(form, constr)
        end
    end
    return
end

## explicit
"""
    isexplicit(formulation, varconstr)
    isexplicit(formulation, varconstrid)

Return `true` if a variable or a constraint is explicit in a formulation; `false` otherwise.
"""
isexplicit(form::Formulation, varid::VarId) = isexplicit(form, getvar(form, varid))
isexplicit(form::Formulation, var::Variable) = var.perendata.is_explicit
isexplicit(form::Formulation, constrid::ConstrId) = isexplicit(form, getconstr(form, constrid))
isexplicit(form::Formulation, constr::Constraint) = constr.perendata.is_explicit

## name
"""
    getname(formulation, varconstr)
    getname(formulation, varconstrid)

Return the name of a variable or a constraint in a formulation.
"""
getname(form::Formulation, varid::VarId) = getvar(form, varid).name
getname(form::Formulation, var::Variable) = var.name
getname(form::Formulation, constrid::ConstrId) = getconstr(form, constrid).name
getname(form::Formulation, constr::Constraint) = constr.name

## branching_priority
"""
    getbranchingpriority(formulation, var)
    getbranchingpriority(formulation, varid)

Return the branching priority of a variable
"""
getbranchingpriority(form::Formulation, varid::VarId) = getvar(form, varid).branching_priority
getbranchingpriority(form::Formulation, var::Variable) = var.branching_priority

# Reset
"""
    reset!(form, var)
    reset!(form, varid)
    reset!(form, constr)
    reset!(form, constraint)

doc todo
"""
function reset!(form::Formulation, var::Variable)
    setcurcost!(form, var, getperencost(form, var))
    setcurlb!(form, var, getperenlb(form, var))
    setcurub!(form, var, getperenub(form, var))
    setcurkind!(form, var, getperenkind(form, var))
    setcurincval!(form, var, getperenincval(form, var))
    if isperenactive(form, var)
        activate!(form, var)
    else
        deactivate!(form, var)
    end
    return
end
reset!(form::Formulation, varid::VarId)  = reset!(form, getvar(form, varid)) 

function reset!(form::Formulation, constr::Constraint)
    setcurrhs!(form, constr, getperenrhs(form, constr))
    setcursense!(form, constr, getperensense(form, constr))
    setcurincval!(form, constr , getperenincval(form, constr))
    if isperenactive(form, constr)
        activate!(form, constr)
    else
        deactivate!(form, constr)
    end
    return
end
reset!(form::Formulation, constrid::ConstrId) = reset!(form, getconstr(form,constrid))
