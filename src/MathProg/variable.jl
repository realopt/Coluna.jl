abstract type AbstractVarData <: AbstractVcData end

mutable struct VarData <: AbstractVcData
    cost::Float64
    lb::Float64
    ub::Float64
    kind::VarKind
    inc_val::Float64
    is_active::Bool
    is_explicit::Bool
end

function _set_bounds_acc_kind!(vdata::VarData, kind::VarKind)
    if kind == Binary
        if vdata.lb < 0
            vdata.lb = 0
        end
        if vdata.ub > 1
            vdata.ub = 0
        end
    elseif kind == Integer
        vdata.lb = ceil(vdata.lb)
        vdata.ub = floor(vdata.ub)
    end
    return
end

"""
    VarData

Information that defines a state of a variable.
"""
function VarData(
    ;cost::Float64 = 0.0, lb::Float64 = 0.0, ub::Float64 = Inf, kind::VarKind = Continuous,
    inc_val::Float64 = -1.0, is_active::Bool = true, is_explicit::Bool = true
)
    vc = VarData(cost, lb, ub, kind, inc_val, is_active, is_explicit)
    _set_bounds_acc_kind!(vc, kind)
    return vc
end

VarData(vd::VarData) = VarData(
    vd.cost, vd.lb, vd.ub, vd.kind, vd.inc_val, vd.is_active, vd.is_explicit
)

"""
    MoiVarRecord

Structure to hold the pointers to the MOI representation of a Coluna Variable.
"""
mutable struct MoiVarRecord
    index::MoiVarIndex
    bounds::MoiVarBound
    kind::MoiVarKind
end

function MoiVarRecord(;index::MoiVarIndex = MoiVarIndex())
    return MoiVarRecord(index, MoiVarBound(), MoiVarKind())
end

getindex(record::MoiVarRecord) = record.index
getbounds(record::MoiVarRecord) = record.bounds
getkind(record::MoiVarRecord) = record.kind

setindex!(record::MoiVarRecord, index::MoiVarIndex) = record.index = index
setbounds!(record::MoiVarRecord, bounds::MoiVarBound) = record.bounds = bounds
setkind!(record::MoiVarRecord, kind::MoiVarKind) = record.kind = kind

"""
    Variable

Representation of a variable in Coluna.
"""
mutable struct Variable <: AbstractVarConstr
    id::Id{Variable}
    name::String
    perendata::VarData
    curdata::VarData
    moirecord::MoiVarRecord
end

const VarId = Id{Variable}

getid(var::Variable) = var.id

function Variable(
    id::VarId, name::String; var_data = VarData(), moi_index::MoiVarIndex = MoiVarIndex()
)
    return Variable(
        id, name, var_data, VarData(var_data),
        MoiVarRecord(index = moi_index)
    )
end
