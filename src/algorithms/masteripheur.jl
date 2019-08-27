struct MasterIpHeuristic <: AbstractAlgorithm end

struct MasterIpHeuristicData
    incumbents::Incumbents
end
MasterIpHeuristicData(S::Type{<:AbstractObjSense}) = MasterIpHeuristicData(Incumbents(S))

struct MasterIpHeuristicRecord <: AbstractAlgorithmRecord
    incumbents::Incumbents
end

function prepare!(::Type{MasterIpHeuristic}, form, node, strategy_rec, params)
    @logmsg LogLevel(-1) "Prepare MasterIpHeuristic."
    return
end

function run!(::Type{MasterIpHeuristic}, form, node, strategy_rec, params)
    @logmsg LogLevel(1) "Applying Master IP heuristic"
    master = getmaster(form)
    algorithm_data = MasterIpHeuristicData(getobjsense(master))
    deactivate!(master, MasterArtVar)
    enforce_integrality!(master)
    opt_result = optimize!(master)
    relax_integrality!(master)
    activate!(master, MasterArtVar)
    set_ip_primal_sol!(algorithm_data.incumbents, getbestprimalsol(opt_result))
    @logmsg LogLevel(1) string("Found primal solution of ", get_ip_primal_bound(algorithm_data.incumbents))
    @logmsg LogLevel(-3) get_ip_primal_sol(algorithm_data.incumbents)
    # Record data 
    set_ip_primal_sol!(node.incumbents, get_ip_primal_sol(algorithm_data.incumbents))
    return MasterIpHeuristicRecord(algorithm_data.incumbents)
end