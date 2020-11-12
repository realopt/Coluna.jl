mutable struct ColGenStabilizationStorage <: AbstractStorage
    basealpha::Float64 # "global" alpha parameter
    curalpha::Float64 # alpha parameter during the current misprice sequence
    nb_misprices::Int64 # number of misprices during the current misprice sequence
    pseudo_dual_bound::DualBound # pseudo dual bound, may be non-valid, f.e. when the pricing problem solved heuristically
    valid_dual_bound::DualBound # valid dual bound
    stabcenter::Union{Nothing, DualSolution} # current stability center, correspond to cur_dual_bound
    newstabcenter::Union{Nothing, DualSolution} # to keep temporarily stab. center after update
    basestabcenter::Union{Nothing, DualSolution} # stability center, corresponding to valid_dual_bound
end

function ColGenStabilizationStorage(master::Formulation) 
    return ColGenStabilizationStorage(
        0.5, 0.0, 0, DualBound(master), DualBound(master), nothing, nothing, nothing
    )
end

smoothing_is_active(storage::ColGenStabilizationStorage) = !iszero(storage.curalpha)

subgradient_is_needed(storage::ColGenStabilizationStorage, smoothparam::Float64) =  
    smoothparam == 1.0 && storage.nb_misprices == 0

mutable struct ColGenStabStorageState <: AbstractStorageState
    alpha::Float64
    dualbound::DualBound
    stabcenter::Union{Nothing, DualSolution}
end

function ColGenStabStorageState(master::Formulation, storage::ColGenStabilizationStorage)
    alpha = storage.basealpha < 0.5 ? 0.5 : storage.basealpha
    return ColGenStabStorageState(alpha, storage.valid_dual_bound, storage.basestabcenter)
end

function restorefromstate!(
    master::Formulation, storage::ColGenStabilizationStorage, state::ColGenStabStorageState
)
    storage.basealpha = state.alpha
    storage.valid_dual_bound = state.dualbound
    storage.basestabcenter = state.stabcenter
    return
end

const ColGenStabilizationStoragePair = (ColGenStabilizationStorage => ColGenStabStorageState)

function init_stab_before_colgen_loop!(storage::ColGenStabilizationStorage)
    storage.stabcenter = storage.basestabcenter
    storage.pseudo_dual_bound = storage.valid_dual_bound
    return
end

function componentwisefunction(in_dual_sol::DualSolution, out_dual_sol::DualSolution, f::Function)
    form = out_dual_sol.model
    constrids = Vector{ConstrId}()
    constrvals = Vector{Float64}()
    value::Float64 = 0.0
    out_next = iterate(out_dual_sol)
    in_next = iterate(in_dual_sol)
    while out_next !== nothing || in_next !== nothing
        if out_next !== nothing 
            ((out_constrid, out_val), out_state) = out_next
            if in_next !== nothing
                ((in_constrid, in_val), in_state) = in_next
                if in_constrid < out_constrid
                    value = f(in_val, 0.0)
                    push!(constrids, in_constrid)
                    push!(constrvals, value)
                    in_next = iterate(in_dual_sol, in_state)
                elseif out_constrid < in_constrid    
                    value = f(0.0, out_val)
                    push!(constrids, out_constrid)
                    push!(constrvals, value)
                    out_next = iterate(out_dual_sol, out_state)    
                else
                    push!(constrids, out_constrid)
                    value = f(in_val, out_val)
                    push!(constrvals, value)
                    in_next = iterate(in_dual_sol, in_state)
                    out_next = iterate(out_dual_sol, out_state)    
                end
            else
                value = f(0.0, out_val)
                push!(constrids, out_constrid)
                push!(constrvals, value)
                out_next = iterate(out_dual_sol, out_state)    
            end
        else    
            ((in_constrid, in_val), in_state) = in_next
            value = f(in_val, 0.0)
            push!(constrids, in_constrid)
            push!(constrvals, value)
            in_next = iterate(in_dual_sol, in_state)
        end
    end
    return (constrids, constrvals)
end

function linear_combination(in_dual_sol::DualSolution, out_dual_sol::DualSolution, coeff::Float64)
    constrids, constrvals = componentwisefunction(
        in_dual_sol, out_dual_sol, 
        (x, y) -> coeff * x + (1.0 - coeff) * y
    )

    form = in_dual_sol.model
    bound = 0.0
    for (i, constrid) in enumerate(constrids)
        bound += constrvals[i] * getcurrhs(form, constrid) 
    end
    return DualSolution(form, constrids, constrvals, bound, UNKNOWN_FEASIBILITY)
end

function update_stab_after_rm_solve!(
    storage::ColGenStabilizationStorage, smoothparam::Float64, lp_dual_sol::DualSolution
)
    iszero(smoothparam) && return lp_dual_sol

    storage.curalpha = 0.0
    storage.nb_misprices = 0

    if storage.stabcenter === nothing
        # cannot do smoothing, just return the current dual solution
        storage.stabcenter = lp_dual_sol
        return lp_dual_sol
    end

    storage.curalpha = smoothparam == 1.0 ? storage.basealpha : smoothparam

    return linear_combination(storage.stabcenter, lp_dual_sol, storage.curalpha)
end

function norm(dualsol::DualSolution)
    product_sum = 0.0
    for (constrid, val) in dualsol
        product_sum += val * val
    end
    return sqrt(product_sum)
end

function update_alpha_automatically!(
    storage::ColGenStabilizationStorage, nb_new_col::Int64, lp_dual_sol::DualSolution{M},  
    smooth_dual_sol::DualSolution{M}, subgradient_contribution::DualSolution{M}
) where {M}    

    master = lp_dual_sol.model 

    # first we calculate the in-sep direction
    constrids, constrvals = componentwisefunction(smooth_dual_sol, storage.stabcenter, -)
    in_sep_direction = DualSolution(master, constrids, constrvals, 0.0, UNKNOWN_FEASIBILITY)
    in_sep_dir_norm = norm(in_sep_direction)

    # we initialize the subgradient with the right-hand-side of all master constraints
    # except the convexity constraints
    constrids = Vector{ConstrId}()
    constrrhs = Vector{Float64}() 
    for (constrid, constr) in getconstrs(master)
        if !(getduty(constrid) <= MasterConvexityConstr) && 
           iscuractive(master, constr) && isexplicit(master, constr)
            push!(constrids, constrid)
            push!(constrrhs, getcurrhs(master, constrid))
        end 
    end
    subgradient = DualSolution(master, constrids, constrrhs, 0.0, UNKNOWN_FEASIBILITY)
    
    # we calculate the subgradient at the sep point
    for (constrid, value) in subgradient_contribution
        subgradient[constrid] = subgradient[constrid] - value
    end
    subgradient_norm = norm(subgradient)

    # we now calculate the angle between the in-sep direction and the subgradient 
    constrids, constrvals = componentwisefunction(in_sep_direction, subgradient, *)
    angle = sum(constrvals) / (in_sep_dir_norm * subgradient_norm)
    if getobjsense(master) == MaxSense
        angle = -angle
    end

    # we modify the alpha parameter based on the calculated angle
    if nb_new_col == 0 || angle > 1e-12 
        storage.basealpha -= 0.1
    elseif angle < -1e-12 && storage.basealpha < 0.999
        storage.basealpha += (1.0 - storage.basealpha) * 0.1
    end   
    return 
end

function update_stab_after_gencols!(
    storage::ColGenStabilizationStorage, smoothparam::Float64, nb_new_col::Int64, 
    lp_dual_sol::DualSolution{M}, smooth_dual_sol::DualSolution{M}, 
    subgradient_contribution::DualSolution{M}
) where {M}

    iszero(smoothparam) && return nothing

    if smoothparam == 1.0 && storage.nb_misprices == 0
        update_alpha_automatically!(
            storage, nb_new_col, lp_dual_sol, smooth_dual_sol, subgradient_contribution
        )
    end

    if nb_new_col > 0 || !smoothing_is_active(storage)
        return nothing
    end

    if smoothparam < 0
        storage.curalpha = 1.0 - (storage.nb_misprices + 1) * (1 - smoothparam)
    else
        storage.curalpha = 1.0 - (1.0 - storage.curalpha) * 2
    end

    storage.nb_misprices += 1

    if storage.nb_misprices > 10 || storage.curalpha <= 0.0
        storage.curalpha = 0.0
        return lp_dual_sol
    end

    return linear_combination(storage.stabcenter, lp_dual_sol, storage.curalpha)
end

function update_stability_center!(
    storage::ColGenStabilizationStorage, dual_sol::DualSolution, 
    valid_lagr_bound::DualBound, pseudo_lagr_bound::DualBound 
)
    if isbetter(valid_lagr_bound, storage.valid_dual_bound)
        storage.basestabcenter = dual_sol
        storage.valid_dual_bound = valid_lagr_bound
    end
    if isbetter(pseudo_lagr_bound, storage.pseudo_dual_bound)
        storage.newstabcenter = dual_sol
        storage.pseudo_dual_bound = pseudo_lagr_bound
    end
    return
end

function update_stab_after_colgen_iteration!(storage::ColGenStabilizationStorage)
    if storage.newstabcenter !== nothing
        storage.stabcenter = storage.newstabcenter
    end
    storage.curalpha = 0.0 
    storage.newstabcenter = nothing
    return
end
