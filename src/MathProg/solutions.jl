# new structures for solutions

# Constructors for Primal & Dual Solutions
const PrimalSolution{M} = Solution{M, VarId, Float64}
const DualSolution{M} = Solution{M, ConstrId, Float64}

function PrimalSolution(
    form::M, decisions::Vector{De}, vals::Vector{Va}, val::Float64, status::SolutionStatus
) where {M<:AbstractFormulation,De,Va}
    return Solution{M,De,Va}(form, decisions, vals, val, status)
end

function DualSolution(
    form::M, decisions::Vector{De}, vals::Vector{Va}, val::Float64, status::SolutionStatus
) where {M<:AbstractFormulation,De,Va}
    return Solution{M,De,Va}(form, decisions, vals, val, status)
end

function Base.isinteger(sol::Solution)
    for (vc_id, val) in sol
        #if getperenkind(sol.model, vc_id) != Continuous
            !isinteger(val) && return false
        #end
    end
    return true
end

isfractional(sol::Solution) = !Base.isinteger(sol)

function contains(sol::PrimalSolution, f::Function)
    for (varid, val) in sol
        f(varid) && return true
    end
    return false
end

function contains(sol::DualSolution, f::Function)
    for (constrid, val) in sol
        f(constrid) && return true
    end
    return false
end

function _assert_same_model(sols::NTuple{N, Solution{M, I, Float64}}) where {N,M,I}
    for i in 2:length(sols)
        sols[i-1].model != sols[i].model && return false
    end
    return true
end

function Base.cat(sols::Solution{M, I, Float64}...) where {M,I}
    _assert_same_model(sols) || error("Cannot concatenate solutions not attached to the same model.")

    ids = I[]
    vals = Float64[]
    for sol in sols, (id, value) in sol
        push!(ids, id)
        push!(vals, value)
    end

    return Solution{M,I,Float64}(
        sols[1].model, ids, vals, sum(getvalue.(sols)), sols[1].status
    )
end

function Base.print(io::IO, form::AbstractFormulation, sol::Solution)
    println(io, "Solution")
    for (id, val) in sol
        println(io, getname(form, id), " = ", val)
    end
    return
end

function Base.show(io::IO, solution::DualSolution{M}) where {M}
    println(io, "Dual solution")
    for (constrid, value) in solution
        println(io, "| ", getname(solution.model, constrid), " = ", value)
    end
    Printf.@printf(io, "└ value = %.2f \n", getvalue(solution))
end

function Base.show(io::IO, solution::PrimalSolution{M}) where {M}
    println(io, "Primal solution")
    for (varid, value) in solution
        println(io, "| ", getname(solution.model, varid), " = ", value)
    end
    Printf.@printf(io, "└ value = %.2f \n", getvalue(solution))
end
