"""
    SelectionCriterion
"""
@enum SelectionCriterion FirstFoundCriterion MostFractionalCriterion 


"""
    BranchingRuleInput

    Input of a branching rule (branching separation algorithm)
    Contains current solution, max number of candidates and local candidate id.
"""
struct BranchingRuleInput <: AbstractInput 
    solution::PrimalSolution 
    isoriginalsol::Bool
    max_nb_candidates::Int64
    criterion::SelectionCriterion
    local_id::Int64
    int_tol::Float64
end

"""
    BranchingRuleOutput

    Input of a branching rule (branching separation algorithm)
    Contains current incumbents, infeasibility status, and the record of its storage.
"""
struct BranchingRuleOutput <: AbstractOutput 
    local_id::Int64
    groups::Vector{BranchingGroup}
end

function BranchingRuleOuput(input::BranchingRuleInput)
    return BranchingRuleOuput(input.local_id, Vector{BranchingGroup}())
end

getlocalid(output::BranchingRuleOutput) = output.local_id
getgroups(output::BranchingRuleOutput) = output.groups

"""
    AbstractBranchingRule

    Branching rules are algorithms which find branching candidates 
    (branching separation algorithms).
"""
abstract type AbstractBranchingRule <: AbstractAlgorithm end

# branching rules are always manager algorithms (they manage storing and restoring storages)
ismanager(algo::AbstractBranchingRule) = true

function run!(
    rule::AbstractBranchingRule, data::AbstractData, input::BranchingRuleInput
)::BranchingRuleOutput
    algotype = typeof(rule)
    error("Method run! in not defined for branching rule $(typeof(rule)), data $(typeof(data)), and input $(typeof(input)).")
end
