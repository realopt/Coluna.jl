"""
    VarConstrBuffer{T<:AbstractVarConstr}

A `VarConstrBuffer{T}` stores the ids of the entities to be added and removed from the formulation where it belongs.
"""
mutable struct VarConstrBuffer{T<:AbstractVarConstr}
    added::Set{Id{T}}
    removed::Set{Id{T}}
end

"""
    VarConstrBuffer{T}() where {T<:AbstractVarConstr}

Constructs an empty `VarConstrBuffer{T}` for entities of type `T`.
"""
VarConstrBuffer{T}() where {T<:AbstractVarConstr} = VarConstrBuffer{T}(Set{T}(), Set{T}())

function add!(buffer::VarConstrBuffer{VC}, id::Id{VC}) where {VC<:AbstractVarConstr}
    !(id in buffer.removed) && push!(buffer.added, id)
    delete!(buffer.removed, id)
    return
end

function remove!(buffer::VarConstrBuffer{VC}, id::Id{VC}) where {VC<:AbstractVarConstr}
    !(id in buffer.added) && push!(buffer.removed, id)
    delete!(buffer.added, id)
    return
end

"""
    FormulationBuffer()

A `FormulationBuffer` stores all changes done to a `Formulation` `f` since last call to `sync_solver!`.
When function `sync_solver!` is called, the optimizer is synched with all changes in FormulationBuffer

When `f` is modified, such modification should not be passed directly to its optimizer, but instead should be passed to `f.buffer`.

The concerned modificatios are:
1. Cost change in a variable
2. Bound change in a variable
3. Right-hand side change in a Constraint
4. Variable is removed
5. Variable is added
6. Constraint is removed
7. Constraint is added
8. Coefficient in the matrix is modified (reset)
"""
mutable struct FormulationBuffer
    changed_cost::Set{Id{Variable}}
    changed_bound::Set{Id{Variable}}
    changed_var_kind::Set{Id{Variable}}
    changed_constr_kind::Set{Id{Constraint}}
    changed_rhs::Set{Id{Constraint}}
    var_buffer::VarConstrBuffer{Variable}
    constr_buffer::VarConstrBuffer{Constraint}
    reset_coeffs::Dict{Pair{Id{Constraint},Id{Variable}},Float64}
    #reset_partial_sols::Dict{Pair{Id{Variable},Id{Variable}},Float64}
end
"""
    FormulationBuffer()

Constructs an empty `FormulationBuffer`.
"""
FormulationBuffer() = FormulationBuffer(
    Set{Id{Variable}}(), Set{Id{Variable}}(), Set{Id{Variable}}(),
    Set{Id{Constraint}}(), Set{Id{Constraint}}(), VarConstrBuffer{Variable}(),
    VarConstrBuffer{Constraint}(),
    Dict{Pair{Id{Constraint},Id{Variable}},Float64}()
    # , Dict{Pair{Id{Variable},Id{Variable}},Float64}()
)

add!(b::FormulationBuffer, varid::VarId) = add!(b.var_buffer, varid)
add!(b::FormulationBuffer, constrid::ConstrId) = add!(b.constr_buffer, constrid)

remove!(buffer::FormulationBuffer, varid::VarId) = remove!(
    buffer.var_buffer, varid
)

remove!(buffer::FormulationBuffer, constrid::ConstrId) = remove!(
    buffer.constr_buffer, constrid
)

function change_rhs!(buffer::FormulationBuffer, constrid::ConstrId)
    push!(buffer.changed_rhs, constrid)
    return
end

function change_kind!(buffer::FormulationBuffer, constrid::ConstrId)
    push!(buffer.changed_constr_kind, constrid)
    return
end

function change_cost!(buffer::FormulationBuffer, varid::VarId)
    push!(buffer.changed_cost, varid)
    return
end

function change_bound!(buffer::FormulationBuffer, varid::VarId)
    push!(buffer.changed_bound, varid)
    return
end

function change_kind!(buffer::FormulationBuffer, varid::VarId)
    push!(buffer.changed_var_kind, varid)
    return
end

function set_matrix_coeff!(buffer::FormulationBuffer, varid::VarId,
                           constrid::ConstrId, new_coeff::Float64)
    buffer.reset_coeffs[Pair(constrid, varid)] = new_coeff
end
