@enum(TerminationStatus, OPTIMAL, TIME_LIMIT, NODE_LIMIT, OTHER_LIMIT, EMPTY_RESULT, NOT_YET_DETERMINED)
@enum(FeasibilityStatus, FEASIBLE, INFEASIBLE, UNKNOWN_FEASIBILITY)

function convert_status(moi_status::MOI.TerminationStatusCode)
    moi_status == MOI.OPTIMAL && return OPTIMAL
    moi_status == MOI.TIME_LIMIT && return TIME_LIMIT
    moi_status == MOI.NODE_LIMIT && return NODE_LIMIT
    moi_status == MOI.OTHER_LIMIT && return OTHER_LIMIT
    return NOT_YET_DETERMINED
end

function convert_status(coluna_status::TerminationStatus)
    coluna_status == OPTIMAL && return MOI.OPTIMAL
    coluna_status == TIME_LIMIT && return MOI.TIME_LIMIT
    coluna_status == NODE_LIMIT && return MOI.NODE_LIMIT
    coluna_status == OTHER_LIMIT && return MOI.OTHER_LIMIT
    return MOI.OTHER_LIMIT
end

### START : TO BE DELETED
#### Old Optimization Result (usefull only to return the result of the call to
### a solver through MOI). Will be removed very soon
mutable struct MoiResult{F<:AbstractFormulation,S<:Coluna.AbstractSense}
    termination_status::TerminationStatus
    feasibility_status::FeasibilityStatus
    primal_bound::PrimalBound{S}
    dual_bound::DualBound{S}
    primal_sols::Vector{PrimalSolution{F}}
    dual_sols::Vector{DualSolution{F}}
end

function MoiResult(model::M) where {M<:AbstractModel}
    S = getobjsense(model)
    return MoiResult{M,S}(
        NOT_YET_DETERMINED, UNKNOWN_FEASIBILITY, PrimalBound(model),
        DualBound(model),
        PrimalSolution{M}[], DualSolution{M}[]
    )
end

function MoiResult(
    model::M, ts::TerminationStatus, fs::FeasibilityStatus; pb = nothing,
    db = nothing, primal_sols = nothing, dual_sols = nothing
) where {M<:AbstractModel}
    S = getobjsense(model)
    return OptimizationState{M,S}(
        ts, fs,
        pb === nothing ? PrimalBound(model) : pb,
        db === nothing ? DualBound(model) : db,
        primal_sols === nothing ? PrimalSolution{M}[] : primal_sols,
        dual_sols === nothing ? DualSolution{M}[] : dual_sols
    )
end

getterminationstatus(res::MoiResult) = res.termination_status
getfeasibilitystatus(res::MoiResult) = res.feasibility_status
isfeasible(res::MoiResult) = res.feasibility_status == FEASIBLE
getprimalbound(res::MoiResult) = res.primal_bound
getdualbound(res::MoiResult) = res.dual_bound
getprimalsols(res::MoiResult) = res.primal_sols
getdualsols(res::MoiResult) = res.dual_sols
nbprimalsols(res::MoiResult) = length(res.primal_sols)
nbdualsols(res::MoiResult) = length(res.dual_sols)

# For documentation : Only unsafe methods must be used to retrieve best
# solutions in the core of Coluna.
unsafe_getbestprimalsol(res::MoiResult) = res.primal_sols[1]
unsafe_getbestdualsol(res::MoiResult) = res.dual_sols[1]
getbestprimalsol(res::MoiResult) = get(res.primal_sols, 1, nothing)
getbestdualsol(res::MoiResult) = get(res.dual_sols, 1, nothing)

setprimalbound!(res::MoiResult, b::PrimalBound) = res.primal_bound = b
setdualbound!(res::MoiResult, b::DualBound) = res.dual_bound = b
setterminationstatus!(res::MoiResult, status::TerminationStatus) = res.termination_status = status
setfeasibilitystatus!(res::MoiResult, status::FeasibilityStatus) = res.feasibility_status = status
result_gap(res::MoiResult) = gap(getprimalbound(res), getdualbound(res))

function add_primal_sol!(res::MoiResult{M,S}, solution::PrimalSolution{M}) where {M,S}
    push!(res.primal_sols, solution)
    pb = PrimalBound{S}(getvalue(solution))
    if isbetter(pb, getprimalbound(res))
        setprimalbound!(res, pb)
    end
    sort!(res.primal_sols; by = x->valueinminsense(PrimalBound{S}(getvalue(x))))
    return
end

function add_dual_sol!(res::MoiResult{M,S}, solution::DualSolution{M}) where {M,S}
    push!(res.dual_sols, solution)
    db = DualBound{S}(getvalue(solution))
    if isbetter(db, getdualbound(res))
        setdualbound!(res, db)
    end
    #sort!(res.dual_sols; by = x->valueinminsense(DualBound{S}(getvalue(x)))))
    return
end

function determine_statuses(res::MoiResult, fully_explored::Bool)
    gap_is_zero = result_gap(res) <= 0.00001
    found_sols = length(getprimalsols(res)) >= 1
    # We assume that gap cannot be zero if no solution was found
    gap_is_zero && @assert found_sols
    found_sols && setfeasibilitystatus!(res, FEASIBLE)
    gap_is_zero && setterminationstatus!(res, OPTIMAL)
    if !found_sols # Implies that gap is not zero
        setterminationstatus!(res, EMPTY_RESULT)
        # Determine if we can prove that is was infeasible
        if fully_explored
            setfeasibilitystatus!(res, INFEASIBLE)
        else
            setfeasibilitystatus!(res, UNKNOWN_FEASIBILITY)
        end
    elseif !gap_is_zero
        setterminationstatus!(res, OTHER_LIMIT)
    end
    return
end

### END : TO BE DELETED


"""
    NoOptimizer <: AbstractOptimizer

Wrapper to indicate that no optimizer is assigned to a `Formulation`
"""
struct NoOptimizer <: AbstractOptimizer end

no_optimizer_builder(args...) = NoOptimizer()

"""
    UserOptimizer <: AbstractOptimizer

Wrapper that is used when the `optimize!(f::Formulation)` function should call an user-defined callback.
"""
mutable struct UserOptimizer <: AbstractOptimizer
    user_oracle::Function
end

mutable struct PricingCallbackData
    form::Formulation
    result::Union{Nothing, MoiResult}
end

function optimize!(form::Formulation, optimizer::UserOptimizer)
    @logmsg LogLevel(-2) "Calling user-defined optimization function."
    cbdata = PricingCallbackData(form, nothing)
    optimizer.user_oracle(cbdata)
    return cbdata.result
end

"""
    MoiOptimizer <: AbstractOptimizer

Wrapper that is used when the optimizer of a formulation 
is an `MOI.AbstractOptimizer`, thus inheriting MOI functionalities.
"""
struct MoiOptimizer <: AbstractOptimizer
    inner::MOI.ModelLike
end

getinner(optimizer::MoiOptimizer) = optimizer.inner

function retrieve_result(form::Formulation, optimizer::MoiOptimizer)
    result = MoiResult(form)
    terminationstatus = MOI.get(getinner(optimizer), MOI.TerminationStatus())
    if terminationstatus != MOI.INFEASIBLE &&
            terminationstatus != MOI.DUAL_INFEASIBLE &&
            terminationstatus != MOI.INFEASIBLE_OR_UNBOUNDED &&
            terminationstatus != MOI.OPTIMIZE_NOT_CALLED
        fill_primal_result!(form, optimizer, result)
        fill_dual_result!(form, optimizer, result)
        if MOI.get(getinner(optimizer), MOI.ResultCount()) >= 1 
            setfeasibilitystatus!(result, FEASIBLE)
            setterminationstatus!(result, convert_status(terminationstatus))
        else
            msg = """
            Termination status = $(terminationstatus) but no results.
            Please, open an issue at https://github.com/atoptima/Coluna.jl/issues
            """
            error(msg)
        end
    else
        @warn "Solver has no result to show."
        setfeasibilitystatus!(result, INFEASIBLE)
        setterminationstatus!(result, EMPTY_RESULT)
    end
    return result
end

function setwarmstart!(form::Formulation, optimizer::MoiOptimizer, sol::PrimalSolution)
    for (varid, val) in sol
        moirec = getmoirecord(getvar(form, varid))
        moi_index = getindex(moirec)
        MOI.set(optimizer.inner, MOI.VariablePrimalStart(), moi_index, val)
    end
    return
end

function optimize!(form::Formulation, optimizer::MoiOptimizer)
    @logmsg LogLevel(-4) "MOI formulation before synch: "
    @logmsg LogLevel(-4) getoptimizer(form)
    sync_solver!(getoptimizer(form), form)
    @logmsg LogLevel(-3) "MOI formulation after synch: "
    @logmsg LogLevel(-3) getoptimizer(form)
    nbvars = MOI.get(form.optimizer.inner, MOI.NumberOfVariables())
    if nbvars <= 0
        @warn "No variable in the formulation. Coluna does not call the solver."
        return retrieve_result(form, optimizer)
    end
    call_moi_optimize_with_silence(form.optimizer)
    status = MOI.get(form.optimizer.inner, MOI.TerminationStatus())
    @logmsg LogLevel(-2) string("Optimization finished with status: ", status)
    return retrieve_result(form, optimizer)
end

function sync_solver!(optimizer::MoiOptimizer, f::Formulation)
    @logmsg LogLevel(-1) string("Synching formulation ", getuid(f))
    buffer = f.buffer
    matrix = getcoefmatrix(f)

    # Remove constrs
    @logmsg LogLevel(-2) string("Removing constraints")
    remove_from_optimizer!(buffer.constr_buffer.removed, f)

    # Remove vars
    @logmsg LogLevel(-2) string("Removing variables")
    remove_from_optimizer!(buffer.var_buffer.removed, f)

    # Add vars
    for id in buffer.var_buffer.added
        v = getvar(f, id)
        @logmsg LogLevel(-4) string("Adding variable ", getname(f, v))
        add_to_optimizer!(f, v)
    end

    # Add constrs
    for constr_id in buffer.constr_buffer.added
        constr = getconstr(f, constr_id)
        @logmsg LogLevel(-4) string("Adding constraint ", getname(f, constr))
        add_to_optimizer!(f, constr, (f, constr) -> iscuractive(f, constr) && iscurexplicit(f, constr))  
    end

    # Update variable costs
    for id in buffer.changed_cost
        (id in buffer.var_buffer.added || id in buffer.var_buffer.removed) && continue
        update_cost_in_optimizer!(f, getvar(f, id))
    end

    # Update variable bounds
    for id in buffer.changed_bound
        (id in buffer.var_buffer.added || id in buffer.var_buffer.removed) && continue
        @logmsg LogLevel(-4) "Changing bounds of variable " getname(f, id)
        @logmsg LogLevel(-5) string("New lower bound is ", getcurlb(f, id))
        @logmsg LogLevel(-5) string("New upper bound is ", getcurub(f, id))
        update_bounds_in_optimizer!(f, getvar(f, id))
    end

    # Update variable kind
    for id in buffer.changed_var_kind
        (id in buffer.var_buffer.added || id in buffer.var_buffer.removed) && continue
        @logmsg LogLevel(-2) "Changing kind of variable " getname(f, id)
        @logmsg LogLevel(-3) string("New kind is ", getcurkind(f, id))
        enforce_kind_in_optimizer!(f, getvar(f,id))
    end

    # Update constraint rhs
    for id in buffer.changed_rhs
        (id in buffer.constr_buffer.added || id in buffer.constr_buffer.removed) && continue
        @logmsg LogLevel(-2) "Changing rhs of constraint " getname(f, id)
        @logmsg LogLevel(-3) string("New rhs is ", getcurrhs(f, id))
        update_constr_rhs_in_optimizer!(f, getconstr(f, id))
    end

    # Update matrix
    # First check if should update members of just-added vars
    matrix = getcoefmatrix(f)
    for id in buffer.var_buffer.added
        for (constrid, coeff) in  matrix[:,id]
            iscuractive(f, constrid) || continue
            iscurexplicit(f, constrid) || continue
            constrid ∉ buffer.constr_buffer.added || continue
            c = getconstr(f, constrid)
            update_constr_member_in_optimizer!(optimizer, c, getvar(f, id), coeff)
        end
    end

    # Then updated the rest of the matrix coeffs
    for ((c_id, v_id), coeff) in buffer.reset_coeffs
        # Ignore modifications involving vc's that were removed
        (c_id in buffer.constr_buffer.removed || v_id in buffer.var_buffer.removed) && continue
        c = getconstr(f, c_id)
        v = getvar(f, v_id)
        @logmsg LogLevel(-2) string("Setting matrix coefficient: (", getname(f, c), ",", getname(f, v), ") = ", coeff)
        update_constr_member_in_optimizer!(optimizer, c, v, coeff)
    end
    _reset_buffer!(f)
    return
end

# Fallbacks
optimize!(f::Formulation, ::S) where {S<:AbstractOptimizer} = error(
    string("Function `optimize!` is not defined for object of type ", S)
)

# Initialization of optimizers
function _initialize_optimizer!(optimizer::MoiOptimizer, form::Formulation)
    f = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 0.0)
    MOI.set(form.optimizer.inner, MoiObjective(), f)
    set_obj_sense!(form.optimizer, getobjsense(form))
    return
end

_initialize_optimizer!(optimizer, form::Formulation) = return