function convert_status(moi_status::MOI.TerminationStatusCode)
    moi_status == MOI.OPTIMAL && return OPTIMAL
    moi_status == MOI.INFEASIBLE && return INFEASIBLE
    moi_status == MOI.TIME_LIMIT && return TIME_LIMIT
    moi_status == MOI.NODE_LIMIT && return NODE_LIMIT
    moi_status == MOI.OTHER_LIMIT && return OTHER_LIMIT
    return UNCOVERED_TERMINATION_STATUS
end

function convert_status(coluna_status::TerminationStatus)
    coluna_status == OPTIMAL && return MOI.OPTIMAL
    coluna_status == INFEASIBLE && return MOI.INFEASIBLE
    coluna_status == TIME_LIMIT && return MOI.TIME_LIMIT
    coluna_status == NODE_LIMIT && return MOI.NODE_LIMIT
    coluna_status == OTHER_LIMIT && return MOI.OTHER_LIMIT
    return MOI.OTHER_LIMIT
end

"""
    NoOptimizer <: AbstractOptimizer

Wrapper to indicate that no optimizer is assigned to a `Formulation`
"""
struct NoOptimizer <: AbstractOptimizer end

no_optimizer_builder(args...) = NoOptimizer()

"""
    UserOptimizer <: AbstractOptimizer

TODO
"""
mutable struct UserOptimizer <: AbstractOptimizer
    user_oracle::Function
end

struct PricingCallbackData
    form::Formulation
    primal_solutions::Vector{PrimalSolution}
end

function PricingCallbackData(form::F) where {F<:Formulation} 
    return PricingCallbackData(form, PrimalSolution{F}[])
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
        @logmsg LogLevel(-2) string("Adding constraint ", getname(f, constr))
        add_to_optimizer!(f, constr, (f, constr) -> iscuractive(f, constr) && isexplicit(f, constr))
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
        @logmsg LogLevel(-3) "Changing kind of variable " getname(f, id)
        @logmsg LogLevel(-4) string("New kind is ", getcurkind(f, id))
        enforce_kind_in_optimizer!(f, getvar(f,id))
    end

    # Update constraint rhs
    for id in buffer.changed_rhs
        (id in buffer.constr_buffer.added || id in buffer.constr_buffer.removed) && continue
        @logmsg LogLevel(-3) "Changing rhs of constraint " getname(f, id)
        @logmsg LogLevel(-4) string("New rhs is ", getcurrhs(f, id))
        update_constr_rhs_in_optimizer!(f, getconstr(f, id))
    end

    # Update matrix
    # First check if should update members of just-added vars
    matrix = getcoefmatrix(f)
    for id in buffer.var_buffer.added
        for (constrid, coeff) in @view matrix[:,id]
            iscuractive(f, constrid) || continue
            isexplicit(f, constrid) || continue
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

# Initialization of optimizers
function _initialize_optimizer!(optimizer::MoiOptimizer, form::Formulation)
    f = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 0.0)
    MOI.set(form.optimizer.inner, MoiObjective(), f)
    set_obj_sense!(form.optimizer, getobjsense(form))
    return
end

_initialize_optimizer!(optimizer, form::Formulation) = return

