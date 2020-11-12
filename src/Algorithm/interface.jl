"""
    About algorithms
    ----------------

    An algorithm is a procedure with a known interface (input and output) applied to a data.
    An algorithm can use storages inside the data to keep its computed data between different
    runs of the algorithm or between runs of different algorithms.
    The algorithm itself contains only its parameters. 

    Parameters of an algorithm may contain its child algorithms which used by it. Therefore, 
    the algoirthm tree is formed, in which the root is the algorithm called to solver the model 
    (root algorithm should be an optimization algorithm, see below). 

    Algorithms are divided into two types : "manager algorithms" and "worker algorithms". 
    Worker algorithms just continue the calculation. They do not store and restore storages 
    as they suppose it is done by their master algorithms. Manager algorithms may divide 
    the calculation flow into parts. Therefore, they store and restore storages to make sure 
    that their child worker algorithms have storages prepared. 
    A worker algorithm cannot have child manager algorithms. 

    Examples of manager algorithms : TreeSearchAlgorithm (which covers both BCP algorithm and 
    diving algorithm), conquer algorithms, strong branching, branching rule algorithms 
    (which create child nodes). Examples of worker algorithms : column generation, SolveIpForm, 
    SolveLpForm, cut separation, pricing algorithms, etc.

"""

"""
    AbstractInput

    Input of an algorithm.     
"""
abstract type AbstractInput end 

struct EmptyInput <: AbstractInput end

"""
    AbstractOutput

Output of an algorithm.     
"""
abstract type AbstractOutput end 

"""
    AbstractAlgorithm

"""
abstract type AbstractAlgorithm end

ismanager(algo::AbstractAlgorithm) = false

"""
    get_child_algorithms(::AbstractAlgorithm, ::AbstractModel)::Vector{Tuple{AbstractAlgorithm, AbstractModel}}

    Every algorithm should communicate its child algorithms and the model to which 
    each child algorithm is applied. 
"""
get_child_algorithms(::AbstractAlgorithm, ::AbstractModel) = Tuple{AbstractAlgorithm, AbstractModel}[]

"""
    get_storages_usage(algo::AbstractAlgorithm, model::AbstractModel)::Vector{Tuple{AbstractModel, StorageTypePair, StorageAccessMode}}

    Every algorithm should communicate the storages it uses (so that these storages 
    are created in the beginning) and the usage mode (read only or read-and-write). Usage mode is needed for 
    in order to restore storages before running a worker algorithm.
"""
get_storages_usage(algo::AbstractAlgorithm, model::AbstractModel) = Tuple{AbstractModel, StorageTypePair, StorageAccessMode}[] 

"""
    run!(algo::AbstractAlgorithm, model::AbstractData, input::AbstractInput)::AbstractOutput

    Runs the algorithm. The storage of the algorithm can be obtained from the data
    Returns algorithm's output.    
"""
function run!(algo::AbstractAlgorithm, data::AbstractData, input::AbstractInput)::AbstractOutput
    error("run! not defined for algorithm $(typeof(algo)), data $(typeof(data)), and input $(typeof(input)).")
end

"""
    OptimizationInput

    Contains OptimizationState
"""
struct OptimizationInput{F,S} <: AbstractInput
    optstate::OptimizationState{F,S}
end

getoptstate(input::OptimizationInput) =  input.optstate

"""
    OptimizationOutput

    Contain OptimizationState, PrimalSolution (solution to relaxation), and 
    DualBound (dual bound value)
"""
struct OptimizationOutput{F,S} <: AbstractOutput
    optstate::OptimizationState{F,S}    
end

getoptstate(output::OptimizationOutput)::OptimizationState = output.optstate

"""
    AbstractOptimizationAlgorithm

    This type of algorithm is used to "bound" a model, i.e. to improve primal
    and dual bounds of the model. Solving to optimality is a special case of "bounding".
    The input of such algorithm should be of type Incumbents.    
    The output of such algorithm should be of type OptimizationState.    
"""
abstract type AbstractOptimizationAlgorithm <: AbstractAlgorithm end

exploits_primal_solutions(algo::AbstractOptimizationAlgorithm) = false

# this function collects storages to restore for an algorithm and all its child worker algorithms,
# child manager algorithms are skipped, as their restore storages themselves
function collect_storages_to_restore!(
    global_storages_usage::StoragesUsageDict, algo::AbstractAlgorithm, model::AbstractModel
)
    local_storages_usage = get_storages_usage(algo, model)
    for (stor_model, stor_pair, stor_usage) in local_storages_usage
        add_storage_pair_usage!(global_storages_usage, stor_model, stor_pair, stor_usage)
    end

    child_algos = get_child_algorithms(algo, model)
    for (childalgo, childmodel) in child_algos
        !ismanager(childalgo) && collect_storages_to_restore!(global_storages_usage, childalgo, childmodel)
    end
end

# this function collects storages to create for an algorithm and all its child algorithms
# this function is used only the function initialize_storages!() below
function collect_storages_to_create!(
    storages_to_create::Dict{AbstractModel,Set{StorageTypePair}}, algo::AbstractAlgorithm, model::AbstractModel
)
    storages_usage = get_storages_usage(algo, model)
    for (stor_model, stor_pair, stor_usage) in storages_usage
        if !haskey(storages_to_create, stor_model)
            storages_to_create[stor_model] = Set{StorageTypePair}()
        end
        push!(storages_to_create[stor_model], stor_pair)
    end

    child_algos = get_child_algorithms(algo, model)
    for (childalgo, childmodel) in child_algos
        collect_storages_to_create!(storages_to_create, childalgo, childmodel)
    end
end

# this function initializes all the storages
function initialize_storages!(data::AbstractData, algo::AbstractOptimizationAlgorithm)
    storages_to_create = Dict{AbstractModel,Set{StorageTypePair}}()
    collect_storages_to_create!(storages_to_create, algo, getmodel(data)) 

    for (model, type_pair_set) in storages_to_create        
        #println(IOContext(stdout, :compact => true), model, " ", type_pair_set)
        ModelType = typeof(model)
        storagedict = get_model_storage_dict(data, model)
        if storagedict === nothing
            error(string("Model of type $(typeof(model)) with id $(getuid(model)) ",
                         "is not contained in $(getnicename(data))")                        
            )
        end   
        for type_pair in type_pair_set
            (StorageType, StorageStateType) = type_pair
            storagedict[type_pair] = 
                StorageContainer{ModelType, StorageType, StorageStateType}(model)
        end
    end
end
