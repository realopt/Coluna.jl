# alias only used in this file
const Primal = AbstractPrimalSpace
const Dual = AbstractDualSpace
const MinSense = AbstractMinSense
const MaxSense = AbstractMaxSense

# Bounds
struct Bound{Space<:AbstractSpace,Sense<:AbstractSense} <: Real
    value::Float64
    Bound{Space,Sense}(x::Number) where {Space,Sense} = new(x === NaN ? _defaultboundvalue(Space, Sense) : x)
end

_defaultboundvalue(::Type{<:Primal}, ::Type{<:MinSense}) = Inf
_defaultboundvalue(::Type{<:Primal}, ::Type{<:MaxSense}) = -Inf
_defaultboundvalue(::Type{<:Dual}, ::Type{<:MinSense}) = -Inf
_defaultboundvalue(::Type{<:Dual}, ::Type{<:MaxSense}) = Inf

"""
    Bound{Space, Sense}

doc todo
"""
function Bound{Space,Sense}() where {Space<:AbstractSpace,Sense<:AbstractSense}
    val = _defaultboundvalue(Space, Sense)
    return Bound{Space,Sense}(val)
end

getvalue(b::Bound) = b.value
Base.float(b::Bound) = b.value

"""
    isbetter

doc todo
"""
isbetter(b1::Bound{Sp,Se}, b2::Bound{Sp,Se}) where {Sp<:Primal,Se<:MinSense} = b1.value < b2.value
isbetter(b1::Bound{Sp,Se}, b2::Bound{Sp,Se}) where {Sp<:Primal,Se<:MaxSense} = b1.value > b2.value
isbetter(b1::Bound{Sp,Se}, b2::Bound{Sp,Se}) where {Sp<:Dual,Se<:MinSense} = b1.value > b2.value
isbetter(b1::Bound{Sp,Se}, b2::Bound{Sp,Se}) where {Sp<:Dual,Se<:MaxSense} = b1.value < b2.value

"""
    diff 

distance to reach the dual bound from the primal bound;
non-positive if dual bound reached.
"""
function Base.diff(pb::Bound{<:Primal,<:MinSense}, db::Bound{<:Dual,<:MinSense})
    return pb.value - db.value
end

function Base.diff(db::Bound{<:Dual,<:MinSense}, pb::Bound{<:Primal,<:MinSense})
    return pb.value - db.value
end

function Base.diff(pb::Bound{<:Primal,<:MaxSense}, db::Bound{<:Dual,<:MaxSense})
    return db.value - pb.value
end

function Base.diff(db::Bound{<:Dual,<:MaxSense}, pb::Bound{<:Primal,<:MaxSense})
    return db.value - pb.value
end

"""
    gap

relative gap. Gap is non-positive if pb reached db
"""
function gap(pb::Bound{<:Primal,<:MinSense}, db::Bound{<:Dual,<:MinSense})
    return diff(pb, db) / abs(db.value)
end

function gap(db::Bound{<:Dual,<:MinSense}, pb::Bound{<:Primal,<:MinSense})
    return diff(pb, db) / abs(db.value)
end

function gap(pb::Bound{<:Primal,<:MaxSense}, db::Bound{<:Dual,<:MaxSense})
    return diff(pb, db) / abs(pb.value)
end

function gap(db::Bound{<:Dual,<:MaxSense}, pb::Bound{<:Primal,<:MaxSense})
    return diff(pb, db) / abs(pb.value)
end

"""
    printbounds

doc todo
"""
function printbounds(db::Bound{<:Dual,S}, pb::Bound{<:Primal,S}) where {S<:MinSense}
    Printf.@printf "[ %.4f , %.4f ]" getvalue(db) getvalue(pb)
end

function printbounds(db::Bound{<:Dual,S}, pb::Bound{<:Primal,S}) where {S<:MaxSense}
    Printf.@printf "[ %.4f , %.4f ]" getvalue(pb) getvalue(db)
end

function Base.show(io::IO, b::Bound)
    print(io, getvalue(b))
end

Base.promote_rule(B::Type{<:Bound}, ::Type{<:AbstractFloat}) = B
Base.promote_rule(B::Type{<:Bound}, ::Type{<:Integer}) = B
Base.promote_rule(B::Type{<:Bound}, ::Type{<:AbstractIrrational}) = B
Base.promote_rule(::Type{<:Bound{<:Primal,Se}}, ::Type{<:Bound{<:Dual,Se}}) where {Se<:AbstractSense} = Float64

Base.convert(::Type{<:AbstractFloat}, b::Bound) = b.value
Base.convert(::Type{<:Integer}, b::Bound) = b.value
Base.convert(::Type{<:AbstractIrrational}, b::Bound) = b.value
Base.convert(B::Type{<:Bound}, f::AbstractFloat) = B(f)
Base.convert(B::Type{<:Bound}, i::Integer) = B(i)
Base.convert(B::Type{<:Bound}, i::AbstractIrrational) = B(i)

Base.:-(b::B) where {B<:Bound} = B(-b.value)
Base.:+(b1::B, b2::B) where {B<:Bound} = B(b1.value + b2.value)
Base.:-(b1::B, b2::B) where {B<:Bound} = B(b1.value - b2.value)
Base.:*(b1::B, b2::B) where {B<:Bound} = B(b1.value * b2.value)
Base.:/(b1::B, b2::B) where {B<:Bound} = B(b1.value / b2.value)
Base.:(==)(b1::B, b2::B) where {B<:Bound} = b1.value == b2.value
Base.:<(b1::B, b2::B) where {B<:Bound} = b1.value < b2.value
Base.:(<=)(b1::B, b2::B) where {B<:Bound} = b1.value <= b2.value
Base.:(>=)(b1::B, b2::B) where {B<:Bound} = b1.value >= b2.value
Base.:>(b1::B, b2::B) where {B<:Bound} = b1.value > b2.value
Base.isapprox(b1::B, b2::B) where {B<:Bound} = isapprox(b1.value, b2.value) # TODO : rm ?
Base.isapprox(b::B, val::Number) where {B<:Bound} = isapprox(b.value, val) # TODO : rm ?
Base.isapprox(val::Number, b::B) where {B<:Bound} = isapprox(b.value, val) # TODO : rm ?

#extremum(bounds::Vector{Bound{Sp,Se}}) where {Sp<:Primal,Se<:MinSense} = minimum(bounds) # TODO : use worst or best instead ?
#extremum(bounds::Vector{Bound{Sp,Se}}) where {Sp<:Dual,Se<:MinSense} = maximum(bounds)
#extremum(bounds::Vector{Bound{Sp,Se}}) where {Sp<:Primal,Se<:MaxSense} = maximum(bounds)
#extremum(bounds::Vector{Bound{Sp,Se}}) where {Sp<:Dual,Se<:MaxSense} = minimum(bounds)

"""
    TerminationStatus

Theses statuses are the possible reasons why an algorithm stopped the optimization. 
When a subsolver is called through MOI, the
MOI [`TerminationStatusCode`](https://jump.dev/MathOptInterface.jl/stable/apireference/#MathOptInterface.TerminationStatusCode)
is translated into a Coluna `TerminationStatus`.

Description of the termination statuses: 
- `OPTIMAL` : the algorithm found a global optimal solution given the optimality tolerance.
- `INFEASIBLE` : the algorithm proved infeasibility
- `TIME_LIMIT` : the algorithm stopped because of the time limit
- `NODE_LIMIT` : the branch-and-bound based algorithm stopped due to the node limit
- `OTHER_LIMIT` : the algorithm stopped because of a limit that is neither the time limit 
nor the node limit

If the algorithm has not been called, the default value of the termination status should be:
- `UNKNOWN_TERMINATION_STATUS`

If the subsolver called through MOI returns a 
`TerminationStatusCode` that is not `MOI.OPTIMAL`, `MOI.INFEASIBLE`, `MOI.TIME_LIMIT`, `MOI.NODE_LIMIT`, or 
`MOI.OTHER_LIMIT`:
- `UNCOVERED_TERMINATION_STATUS` : should not be used by a Coluna algorithm
"""
@enum(
    TerminationStatus, OPTIMAL, INFEASIBLE, TIME_LIMIT, NODE_LIMIT, OTHER_LIMIT, 
    UNKNOWN_TERMINATION_STATUS, UNCOVERED_TERMINATION_STATUS
)

"""
    SolutionStatus

todo
"""
@enum(
    SolutionStatus, FEASIBLE_SOL, INFEASIBLE_SOL, UNKNOWN_FEASIBILITY, UNKNOWN_SOLUTION_STATUS,
    UNCOVERED_SOLUTION_STATUS
)

# Solution
struct Solution{Model<:AbstractModel,Decision,Value} <: AbstractDict{Decision,Value}
    model::Model
    bound::Float64
    status::SolutionStatus
    sol::DynamicSparseArrays.PackedMemoryArray{Decision,Value}
end

"""
    Solution

doc todo. Solution is immutable.
"""
function Solution{Mo,De,Va}(model::Mo, decisions::Vector{De}, vals::Vector{Va}, value::Float64, status::SolutionStatus) where {Mo<:AbstractModel,De,Va}
    sol = DynamicSparseArrays.dynamicsparsevec(decisions, vals)
    return Solution(model, value, status, sol)
end

getsol(s::Solution) = s.sol
getvalue(s::Solution) = float(s.bound)

Base.iterate(s::Solution) = iterate(s.sol)
Base.iterate(s::Solution, state) = iterate(s.sol, state)
Base.length(s::Solution) = length(s.sol)
Base.get(s::Solution{Mo,De,Va}, id::De, default) where {Mo,De,Va} = s.sol[id] # TODO
Base.setindex!(s::Solution{Mo,De,Va}, val::Va, id::De) where {Mo,De,Va} = s.sol[id] = val

# todo: move in DynamicSparseArrays or avoid using filter ?
function Base.filter(f::Function, pma::DynamicSparseArrays.PackedMemoryArray{K,T,P}) where {K,T,P}
    ids = Vector{K}()
    vals = Vector{T}()
    for e in pma
        if f(e)
            push!(ids, e[1])
            push!(vals, e[2])
        end
    end
    return DynamicSparseArrays.dynamicsparsevec(ids, vals)
end

function Base.filter(f::Function, s::S) where {S <: Solution}
    return S(s.model, s.bound, s.status, filter(f, s.sol))
end

function Base.show(io::IO, solution::Solution{Mo,De,Va}) where {Mo,De,Va}
    println(io, "Solution")
    for (decision, value) in solution
        println(io, "| ", decision, " = ", value)
    end
    Printf.@printf(io, "└ value = %.2f \n", getvalue(solution))
end

Base.copy(s::S) where {S<:Solution} = S(s.bound, copy(s.sol))