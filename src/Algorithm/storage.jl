@enum(StorageAccessMode, READ_AND_WRITE, READ_ONLY, NOT_USED)

"""
    About storages
    --------------

    Storages keep user data (a model) and computed data between different runs 
    of an algorithm or between runs of different algorithms. 
    Models are storages themselves. Each storage is associated with a
    model. Thus a storage adds computed data to a model.  

    Storage states are useful to store states of storages at some point 
    of the calculation flow so that we can later return to this point and 
    restore the storages. For example, the calculation flow may return to
    some saved node in the search tree.

    Some storages can have different parts which are stored in different 
    storage states. Thus, we operate with triples (model, storage, storage state).
    For every model there may be only one storage for each couple 
    (storage type, storage state type). 
    
    To store all storages of a data, we use functions 
    "store_states!(::AbstractData)::StorageStatesVector" or
    "copy_states(::StorageStatesVector)::StorageStatesVector"
  
    Every stored state should be removed or restored using functions 
    "restore_states!(::StorageStatesVector,::StoragesUsageDict)" 
    and "remove_states!(::StorageStatesVector)"
  
    After storing current states, if we write to some storage, we should restore 
    it for writing using "restore_states!(...)" 
    After storing current states, if we read from a storage, 
    no particular precautions should be taken.   
"""

"""
    AbstractStorage 

    For every storage a constructor should be defined which
    takes a model as a parameter. This constructor is 
    called when the formulation is completely known so the data
    can be safely computed.
"""

abstract type AbstractStorage end

"""
    AbstractStorageState
    
    For each storage state, a constructor should be defined which
    takes a model and a storage as parameters. This constructor
    is called during storing a storage. 
    
"""

abstract type AbstractStorageState end

"""
    function restorefromstate!(model, storage, storage state)
    
    Should be defined for every triple (model type, storage type, storage state type)
    used by an algorithm.     
"""

restorefromstate!(model::AbstractModel, storage::AbstractStorage, state::AbstractStorageState) =
    error(string(
        "Method restorefromstate!() is not defined for model type $(typeof(model)), ",
        "storage type $(typeof(storage)), and storage state type $(typeof(state))"
    ))    


"""
    EmptyStorageState

    If a storage is not changed after initialization, then 
    the empty storage state should be used with it.
"""

struct EmptyStorageState <: AbstractStorageState end

EmptyStorageState(model::AbstractModel, storage::AbstractStorage) = nothing

restorefromstate!(::AbstractModel, ::AbstractStorage, ::EmptyStorageState) = nothing

# StorageTypePair = Pair{Type{<:AbstractStorage}, Type{<:AbstractStorageState}}.
# see https://github.com/atoptima/Coluna.jl/pull/323#discussion_r418972805
const StorageTypePair = Pair{DataType, DataType}

# TO DO : replace with the set of StorageTypePair, should only contain storages which should 
#         be restored for writing (all other storages are restored anyway but just for reading)
const StoragesUsageDict = Dict{Tuple{AbstractModel, StorageTypePair}, StorageAccessMode}


"""
    function add_storage_pair_usage!(::StoragesUsageDict, ::AbstractModel, ::StorageTypePair, ::StorageAccessMode)

    An auxiliary function to be used when adding storage usage to a StorageUsageDict
"""
function add_storage_pair_usage!(
    dict::StoragesUsageDict, model::AbstractModel, pair::StorageTypePair, mode::StorageAccessMode
)
    current_mode = get(dict, (model, pair), NOT_USED) 
    if current_mode == NOT_USED && mode != NOT_USED
        dict[(model, pair)] = mode
    else
        if mode == READ_AND_WRITE && current_mode == READ_ONLY
            dict[(model, pair)] = READ_AND_WRITE
        end    
    end
end

"""
    StorageStateContainer

    This container keeps additional storage state information needed for 
    keeping the number of times the state has been stored. When 
    this number drops to zero, the state can be deleted. 
"""

const StateId = Int

mutable struct StorageStateContainer{SS<:AbstractStorageState}
    id::StateId
    participation::Int
    state::Union{Nothing, SS}
end

StorageStateContainer{SS}(stateid::StateId, participation::Int) where {SS<:AbstractStorageState} =
    StorageStateContainer{SS}(stateid, participation, nothing)

getstateid(ssc::StorageStateContainer) = ssc.id
stateisempty(ssc::StorageStateContainer) = ssc.state === nothing
getparticipation(ssc::StorageStateContainer) = ssc.participation
getstate(ssc::StorageStateContainer) = ssc.state
increaseparticipation!(ssc::StorageStateContainer) = ssc.participation += 1
decreaseparticipation!(ssc::StorageStateContainer) = ssc.participation -= 1

function setstate!(ssc::StorageStateContainer{SS}, state_to_set::SS) where {SS<:AbstractStorageState}
    ssc.state = state_to_set
end

function Base.show(io::IO, statecont::StorageStateContainer{SS}) where {SS<:AbstractStorageState}
    print(io, "state ", remove_until_last_point(string(SS)))
    print(io, " with id=", getstateid(statecont), " part=", getparticipation(statecont))
    if getstate(statecont) === nothing
        print(io, " empty")
    else
        print(io, " ", getstate(statecont))
    end
end

"""
    EmptyStorageStateContainer
"""

const EmptyStorageStateContainer = StorageStateContainer{EmptyStorageState}

EmptyStorageStateContainer(stateid::StateId, participation::Int) =
    EmptyStorageStateContainer(1, 0)

getstateid(essc::EmptyStorageStateContainer) = 1
stateisempty(essc::EmptyStorageStateContainer) = true 
getparticipation(essc::EmptyStorageStateContainer) = 0
increaseparticipation!(essc::EmptyStorageStateContainer) = nothing
decreaseparticipation!(essc::EmptyStorageStateContainer) = nothing

"""
    StorageContainer

    This container keeps storages and all states which have been 
    stored. It implements storing and restoring states in an 
    efficient way. 
"""

mutable struct StorageContainer{M<:AbstractModel, S<:AbstractStorage, SS<:AbstractStorageState}
    model::M
    curstatecont::StorageStateContainer{SS}
    maxstateid::StateId
    storage::S
    typepair::StorageTypePair
    statesdict::Dict{StateId, StorageStateContainer{SS}}
end 

const StorageStatesVector = Vector{Pair{StorageContainer, StateId}}

const StorageDict = Dict{StorageTypePair, StorageContainer}

function StorageContainer{M,S,SS}(model::M) where {M,S,SS}
    return StorageContainer{M,S,SS}(
        model, StorageStateContainer{SS}(1, 0), 1, S(model), 
        S => SS, Dict{StateId, StorageStateContainer{SS}}()
    )
end    

getmodel(sc::StorageContainer) = sc.model
getcurstatecont(sc::StorageContainer) = sc.curstatecont
getmaxstateid(sc::StorageContainer) = sc.maxstateid
getstatesdict(sc::StorageContainer) = sc.statesdict
getstorage(sc::StorageContainer) = sc.storage
gettypepair(sc::StorageContainer) = sc.typepair

function Base.show(io::IO, storagecont::StorageContainer)
    print(io, "storage (")
    print(IOContext(io, :compact => true), getmodel(storagecont))
    (StorageType, StorageStateType) = gettypepair(storagecont)    
    print(io, ", ", remove_until_last_point(string(StorageType)))    
    print(io, ", ", remove_until_last_point(string(StorageStateType)), ")")        
end

function setcurstate!(
    storagecont::StorageContainer{M,S,SS}, statecont::StorageStateContainer{SS}
) where {M,S,SS} 
    # we delete the current state container from the dictionary if necessary
    curstatecont = getcurstatecont(storagecont)
    if !stateisempty(curstatecont) && getparticipation(curstatecont) == 0
        delete!(getstatesdict(storagecont), getstateid(curstatecont))
        @logmsg LogLevel(-2) string("Removed state with id ", getstateid(curstatecont), " for ", storagecont)
    end
    storagecont.curstatecont = statecont
    if getmaxstateid(storagecont) < getstateid(statecont) 
        storagecont.maxstateid = getstateid(statecont)
    end
end

function increaseparticipation!(storagecont::StorageContainer, stateid::StateId)
    statecont = getcurstatecont(storagecont)
    if (getstateid(statecont) == stateid)
        increaseparticipation!(statecont)
    else
        statesdict = getstatesdict(storagecont)
        if !haskey(statesdict, stateid) 
            error(string("State with id $stateid does not exist for ", storagecont))
        end
        increaseparticipation!(statesdict[stateid])
    end
end

function retrieve_from_statesdict(storagecont::StorageContainer, stateid::StateId)
    statesdict = getstatesdict(storagecont)
    if !haskey(statesdict, stateid)
        error(string("State with id $stateid does not exist for ", storagecont))
    end
    statecont = statesdict[stateid]
    decreaseparticipation!(statecont)
    if getparticipation(statecont) < 0
        error(string("Participation is below zero for state with id $stateid of ", storagecont))
    end
    return statecont
end

function save_to_statesdict!(
    storagecont::StorageContainer{M,S,SS}, statecont::StorageStateContainer{SS}
) where {M,S,SS}
    if getparticipation(statecont) > 0 && stateisempty(statecont)
        state = SS(getmodel(storagecont), getstorage(storagecont))
        @logmsg LogLevel(-2) string("Created state with id ", getstateid(statecont), " for ", storagecont)
        setstate!(statecont, state)
        statesdict = getstatesdict(storagecont)
        statesdict[getstateid(statecont)] = statecont
    end
end

function storestate!(storagecont::StorageContainer)::StateId 
    statecont = getcurstatecont(storagecont)
    increaseparticipation!(statecont)
    return getstateid(statecont)
end

function restorestate!(
    storagecont::StorageContainer{M,S,SS}, stateid::StateId, mode::StorageAccessMode
) where {M,S,SS}
    statecont = getcurstatecont(storagecont)
    if getstateid(statecont) == stateid 
        decreaseparticipation!(statecont)
        if getparticipation(statecont) < 0
            error(string("Participation is below zero for state with id $stateid of ", getnicename(storagecont)))
        end
        if mode == READ_AND_WRITE 
            save_to_statesdict!(storagecont, statecont)
            statecont = StorageStateContainer{SS}(getmaxstateid(storagecont) + 1, 0)
            setcurstate!(storagecont, statecont)
        end
        return
    elseif mode != NOT_USED
        # we save current state to dictionary if necessary
        save_to_statesdict!(storagecont, statecont)
    end

    statecont = retrieve_from_statesdict(storagecont, stateid)

    if mode == NOT_USED
        if !stateisempty(statecont) && getparticipation(statecont) == 0
            delete!(getstatesdict(storagecont), getstateid(statecont))
            @logmsg LogLevel(-2) string("Removed state with id ", getstateid(statecont), " for ", storagecont)
        end
    else 
        restorefromstate!(getmodel(storagecont), getstorage(storagecont), getstate(statecont))
        @logmsg LogLevel(-2) string("Restored state with id ", getstateid(statecont), " for ", storagecont)
        if mode == READ_AND_WRITE 
            statecont = StorageStateContainer{SS}(getmaxstateid(storagecont) + 1, 0)
        end 
        setcurstate!(storagecont, statecont)
    end
end

"""
    Storage functions used by Coluna
"""

# this is a "lighter" alternative to restore_states!() function below
# not used for the moment as it has impact on the code readability
# we keep this function for a while for the case when function restore_states!()
# happens to be a bottleneck
# function reserve_for_writing!(storagecont::StorageContainer{M,S,SS}) where {M,S,SS}
#     statecont = getcurstatecont(storagecont)
#     save_to_statesdict!(storagecont, statecont)
#     statecont = StorageStateContainer{SS}(getmaxstateid(storagecont) + 1, 0)
#     setcurstate!(storagecont, statecont)
# end

function restore_states!(ssvector::StorageStatesVector, storages_to_restore::StoragesUsageDict)
    TO.@timeit Coluna._to "Restore/remove states" begin
        for (storagecont, stateid) in ssvector
            mode = get(
                storages_to_restore, 
                (getmodel(storagecont), gettypepair(storagecont)), 
                READ_ONLY
            )
            restorestate!(storagecont, stateid, mode)
        end
    end    
    empty!(ssvector) # vector of states should be emptied 
end

function remove_states!(ssvector::StorageStatesVector)
    TO.@timeit Coluna._to "Restore/remove states" begin
        for (storagecont, stateid) in ssvector
            restorestate!(storagecont, stateid, NOT_USED)
        end
    end
    empty!(ssvector) # vector of states should be emptied 
end

function copy_states(states::StorageStatesVector)::StorageStatesVector
    statescopy = StorageStatesVector()
    for (storagecont, stateid) in states
        push!(statescopy, storagecont => stateid)
        increaseparticipation!(storagecont, stateid)
    end
    return statescopy
end

function check_storage_states_participation(storagecont::StorageContainer)
    curstatecont = getcurstatecont(storagecont)
    if getparticipation(curstatecont) > 0
        @warn string("Positive participation of state ", curstatecont)
    end
    statesdict = getstatesdict(storagecont)
    for (stateid, statecont) in statesdict
        if getparticipation(statecont) > 0
            @warn string("Positive participation of state ", statecont)
        end
    end
end

"""
    IMPORTANT!

    Every stored or copied state should be either restored or removed so that it's 
    participation is correctly computed and memory correctly controlled
"""


