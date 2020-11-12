const PrimalBound{S} = Bound{Primal, S}
const DualBound{S} = Bound{Dual, S}

function PrimalBound(form::AbstractFormulation)
    Se = getobjsense(form)
    return Bound{Primal,Se}()
end

function PrimalBound(form::AbstractFormulation, val::Float64)
    Se = getobjsense(form)
    return Bound{Primal,Se}(val)
end

function PrimalBound(form::AbstractFormulation, pb::PrimalBound{S}) where {S}
    Se = getobjsense(form)
    if Se != S
        msg = """
        Cannot create primal bound.
        Sense of the formulation is $Se and sense of the bound is $S.
        """
        error(msg)
    end
    return Bound{Primal,Se}(getvalue(pb))
end

function DualBound(form::AbstractFormulation)
    Se = getobjsense(form)
    return Bound{Dual,Se}()
end

function DualBound(form::AbstractFormulation, val::Float64)
    Se = getobjsense(form)
    return Bound{Dual,Se}(val)
end

function DualBound(form::AbstractFormulation, db::DualBound{S}) where {S}
    Se = getobjsense(form)
    if Se != S
        msg = """
        Cannot create primal bound.
        Sense of the formulation is $Se and sense of the bound is $S.
        """
        error(msg)
    end
    return Bound{Dual,Se}(getvalue(db))
end

valueinminsense(b::PrimalBound{MinSense}) = b.value
valueinminsense(b::DualBound{MinSense}) = b.value
valueinminsense(b::PrimalBound{MaxSense}) = -b.value
valueinminsense(b::DualBound{MaxSense}) = -b.value

# ObjValues
mutable struct ObjValues{S}
    lp_primal_bound::PrimalBound{S}
    lp_dual_bound::DualBound{S}
    ip_primal_bound::PrimalBound{S}
    ip_dual_bound::DualBound{S}
end

"""
    ObjValues(form)

todo
"""
function ObjValues(
    form::M; 
    ip_primal_bound = nothing,
    ip_dual_bound = nothing,
    lp_primal_bound = nothing,
    lp_dual_bound = nothing
) where {M<:AbstractFormulation}
    S = getobjsense(form)
    ov = ObjValues{S}(
        PrimalBound(form), DualBound(form), PrimalBound(form), DualBound(form)
    )
    if ip_primal_bound !== nothing
        set_ip_primal_bound!(ov, PrimalBound(form, ip_primal_bound))
    end
    if ip_dual_bound !== nothing
        set_ip_dual_bound!(ov, DualBound(ip_dual_bound))
    end
    if lp_primal_bound !== nothing
        set_lp_primal_bound!(ov, PrimalBound(lp_primal_bound))
    end
    if lp_dual_bound !== nothing
        set_lp_dual_bound!(ov, DualBound(lp_dual_bound))
    end
    return ov
end

## Getters bounds
"Return the best primal bound of the mixed-integer program."
get_ip_primal_bound(ov::ObjValues) = ov.ip_primal_bound

"Return the best dual bound of the mixed-integer program."
get_ip_dual_bound(ov::ObjValues) = ov.ip_dual_bound

"Return the best primal bound of the linear program."
get_lp_primal_bound(ov::ObjValues) = ov.lp_primal_bound

"Return the best dual bound of the linear program."
get_lp_dual_bound(ov::ObjValues) = ov.lp_dual_bound

## Gaps
"""
Return the gap between the best primal and dual bounds of the integer program.
Should not be used to check convergence
"""
ip_gap(ov::ObjValues) = gap(get_ip_primal_bound(ov), get_ip_dual_bound(ov))

"Return the gap between the best primal and dual bounds of the linear program."
lp_gap(ov::ObjValues) = gap(get_lp_primal_bound(ov), get_lp_dual_bound(ov))

"Return true if the gap between the best primal and dual bounds of the integer program is closed given optimality tolerances."
function ip_gap_closed(
    ov::ObjValues; atol = Coluna.DEF_OPTIMALITY_ATOL, rtol = Coluna.DEF_OPTIMALITY_RTOL
)
    closed = ip_gap(ov) <= 0
    closed = closed || isapprox(
        get_ip_primal_bound(ov).value, get_ip_dual_bound(ov).value, atol = atol, rtol = rtol
    )
    return closed
end

"Return true if the gap between the best primal and dual bounds of the linear program is closed given optimality tolerances."
function lp_gap_closed(
    ov::ObjValues; atol = Coluna.DEF_OPTIMALITY_ATOL, rtol = Coluna.DEF_OPTIMALITY_RTOL
)
    closed = lp_gap(ov) <= 0
    closed = closed || isapprox(
        get_lp_primal_bound(ov).value, get_lp_dual_bound(ov).value, atol = atol, rtol = rtol
    )
    return closed
end

function set_lp_primal_bound!(ov::ObjValues{S}, b::PrimalBound{S}) where {S}
    ov.lp_primal_bound = b
    return
end

function set_lp_dual_bound!(ov::ObjValues{S}, b::DualBound{S}) where {S}
    ov.lp_dual_bound = b
    return
end

function set_ip_primal_bound!(ov::ObjValues{S}, b::PrimalBound{S}) where {S}
    ov.ip_primal_bound = b
    return
end

function set_ip_dual_bound!(ov::ObjValues{S}, b::DualBound{S}) where {S}
    ov.ip_dual_bound = b
    return
end

"""
Update the primal bound of the linear program if the new one is better than the
current one according to the objective sense.
"""
function update_lp_primal_bound!(ov::ObjValues{S}, b::PrimalBound{S}) where {S}
    if isbetter(b, get_lp_primal_bound(ov))
        ov.lp_primal_bound = b
        return true
    end
    return false
end

"""
Update the dual bound of the linear program if the new one is better than the 
current one according to the objective sense.
"""
function update_lp_dual_bound!(ov::ObjValues{S}, b::DualBound{S}) where {S}
    if isbetter(b, get_lp_dual_bound(ov))
        ov.lp_dual_bound = b
        return true
    end
    return false
end

"""
Update the primal bound of the mixed-integer program if the new one is better
than the current one according to the objective sense.
"""
function update_ip_primal_bound!(ov::ObjValues{S}, b::PrimalBound{S}) where {S}
    if isbetter(b, get_ip_primal_bound(ov))
        ov.ip_primal_bound = b
        return true
    end
    return false
end

"""
Update the dual bound of the mixed-integer program if the new one is better than
the current one according to the objective sense.
"""
function update_ip_dual_bound!(ov::ObjValues{S}, b::DualBound{S}) where {S}
    if isbetter(b, get_ip_dual_bound(ov))
        ov.ip_dual_bound = b
        return true
    end
    return false
end

function update!(dest::ObjValues{S}, src::ObjValues{S}) where {S}
    update_ip_dual_bound!(dest, get_ip_dual_bound(src))
    update_ip_primal_bound!(dest, get_ip_primal_bound(src))
    update_lp_dual_bound!(dest, get_lp_dual_bound(src))
    update_lp_primal_bound!(dest, get_lp_primal_bound(src))
    return
end

function Base.show(io::IO, ov::ObjValues{S}) where {S}
    println(io, "ObjValues{", S, "}:")
    println(io, "ip_primal_bound : ", ov.ip_primal_bound)
    println(io, "ip_dual_bound : ", ov.ip_dual_bound)
    println(io, "lp_primal_bound : ", ov.lp_primal_bound)
    println(io, "lp_dual_bound : ", ov.lp_dual_bound)
end