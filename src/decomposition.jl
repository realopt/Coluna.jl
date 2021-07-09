function set_glob_art_var(form::Formulation, is_pos::Bool, env::Env)
    name = string("global_", (is_pos ? "pos" : "neg"), "_art_var")
    cost = env.params.global_art_var_cost
    cost *= getobjsense(form) == MinSense ? 1.0 : -1.0
    return setvar!(
        form, name, MasterArtVar;
        cost = cost, lb = 0.0, ub = Inf, kind = Continuous
    )
end

function create_global_art_vars!(masterform::Formulation, env::Env)
    global_pos = set_glob_art_var(masterform, true, env)
    global_neg = set_glob_art_var(masterform, false, env)
    matrix = getcoefmatrix(masterform)
    for (constrid, constr) in getconstrs(masterform)
        iscuractive(masterform, constrid) || continue
        getduty(constrid) <= AbstractMasterOriginConstr || continue
        if getcursense(masterform, constr) == Greater
            matrix[constrid, getid(global_pos)] = 1.0
        elseif getcursense(masterform, constr) == Less
            matrix[constrid, getid(global_neg)] = -1.0
        else # Equal
            matrix[constrid, getid(global_pos)] = 1.0
            matrix[constrid, getid(global_neg)] = -1.0
        end
    end
end

function instantiatemaster!(
    env::Env, prob::Problem, reform::Reformulation, ::Type{BD.Master},
    ::Type{BD.DantzigWolfe}
)
    form = create_formulation!(
        env,
        MathProg.DwMaster();
        parent_formulation = reform,
        obj_sense = getobjsense(get_original_formulation(prob))
    )
    setobjconst!(form, getobjconst(get_original_formulation(prob)))
    setmaster!(reform, form)
    return form
end

function instantiatemaster!(
    env::Env, prob::Problem, reform::Reformulation, ::Type{BD.Master}, ::Type{BD.Benders}
)
    masterform = create_formulation!(
        env,
        MathProg.BendersMaster();
        parent_formulation = reform,
        obj_sense = getobjsense(get_original_formulation(prob))
    )
    setmaster!(reform, masterform)
    return masterform
end

function instantiatesp!(
    env::Env, reform::Reformulation, masterform::Formulation{DwMaster},
    ::Type{BD.DwPricingSp}, ::Type{BD.DantzigWolfe}
)
    spform = create_formulation!(
        env,
        MathProg.DwSp(nothing, 1, 1, Integ);
        parent_formulation = masterform,
        obj_sense = getobjsense(masterform)
    )
    add_dw_pricing_sp!(reform, spform)
    return spform
end

function instantiatesp!(
    env::Env, reform::Reformulation, masterform::Formulation{BendersMaster},
    ::Type{BD.BendersSepSp}, ::Type{BD.Benders}
)
    spform = create_formulation!(
        env,
        MathProg.BendersSp();
        parent_formulation = masterform,
        obj_sense = getobjsense(masterform)
    )
    add_benders_sep_sp!(reform, spform)
    return spform
end

# Master of Dantzig-Wolfe decomposition
function instantiate_orig_vars!(
    masterform::Formulation{DwMaster},
    origform::Formulation,
    annotations::Annotations,
    mast_ann
)
    vars_per_ann = annotations.vars_per_ann
    for (ann, vars) in vars_per_ann
        formtype = BD.getformulation(ann)
        if formtype <: BD.Master
            for (_, var) in vars
                clonevar!(origform, masterform, masterform, var, MasterPureVar, is_explicit = true)
            end
        end
    end
    return
end

function instantiate_orig_constrs!(
    masterform::Formulation{DwMaster},
    origform::Formulation{Original},
    env::Env,
    annotations::Annotations,
    mast_ann
)
    !haskey(annotations.constrs_per_ann, mast_ann) && return
    constrs = annotations.constrs_per_ann[mast_ann]
    for (_, constr) in constrs
        cloneconstr!(
            origform, masterform, masterform, constr, MasterMixedConstr, 
            loc_art_var_abs_cost = env.params.local_art_var_cost
        ) # TODO distinguish Pure versus Mixed
    end
    # Cut generation callbacks
    for constrgen in get_robust_constr_generators(origform)
        set_robust_constr_generator!(masterform, constrgen.kind, constrgen.separation_alg)
    end
    return
end

function create_side_vars_constrs!(
    masterform::Formulation{DwMaster},
    origform::Formulation{Original},
    env::Env,
    annotations::Annotations
)
    coefmatrix = getcoefmatrix(masterform)
    for (spuid, spform) in get_dw_pricing_sps(masterform.parent_formulation)
        ann = get(annotations, spform)
        setupvars = filter(v -> getduty(v.first) == DwSpSetupVar, getvars(spform))
        @assert length(setupvars) == 1
        setupvar = collect(values(setupvars))[1]
        setuprepvar = clonevar!(origform, masterform, spform, setupvar, MasterRepPricingSetupVar, is_explicit = false)

        # create convexity constraint & storing information about the convexity constraint
        # in the duty data of the formulation
        lb_mult = Float64(BD.getlowermultiplicity(ann))
        name = string("sp_lb_", spuid)
        lb_conv_constr = setconstr!(
            masterform, name, MasterConvexityConstr;
            rhs = lb_mult, kind = Essential, sense = Greater, inc_val = 100.0, 
            loc_art_var_abs_cost = env.params.local_art_var_cost
        )
        masterform.parent_formulation.dw_pricing_sp_lb[spuid] = getid(lb_conv_constr)
        coefmatrix[getid(lb_conv_constr), getid(setuprepvar)] = 1.0

        ub_mult =  Float64(BD.getuppermultiplicity(ann))
        name = string("sp_ub_", spuid)
        ub_conv_constr = setconstr!(
            masterform, name, MasterConvexityConstr; rhs = ub_mult,
            kind = Essential, sense = Less, inc_val = 100.0, 
            loc_art_var_abs_cost = env.params.local_art_var_cost
        )
        masterform.parent_formulation.dw_pricing_sp_ub[spuid] = getid(ub_conv_constr)
        coefmatrix[getid(ub_conv_constr), getid(setuprepvar)] = 1.0

        spform.duty_data.lower_multiplicity = lb_mult
        spform.duty_data.upper_multiplicity = ub_mult
        spform.duty_data.setup_var = getid(setupvar)

        # If pricing subproblem variables are continuous, the master columns generated by
        # the subproblem must have a continuous perenkind.
        # This piece of information is stored in the duty data of the formulation.
        continuous_columns = true
        for (varid, var) in getvars(spform)
            if getduty(varid) <= DwSpPricingVar && getperenkind(spform, var) !== Continuous
                continuous_columns = false
                break
            end
        end
        spform.duty_data.column_var_kind = continuous_columns ? Continuous : Integ
    end
    return
end

function create_artificial_vars!(masterform::Formulation{DwMaster}, env::Env)
    create_global_art_vars!(masterform, env)
    return
end

# Pricing subproblem of Danztig-Wolfe decomposition
function instantiate_orig_vars!(
    spform::Formulation{DwSp},
    origform::Formulation{Original},
    annotations::Annotations,
    sp_ann
)
    !haskey(annotations.vars_per_ann, sp_ann) && return
    vars = annotations.vars_per_ann[sp_ann]
    masterform = spform.parent_formulation
    for (id, var) in vars
        # An original variable annotated in a subproblem is a DwSpPricingVar
        clonevar!(origform, spform, spform, var, DwSpPricingVar, is_explicit = true)
        clonevar!(origform, masterform, spform, var, MasterRepPricingVar,
                  is_explicit = false)#, members = getcoefmatrix(origform)[:,id])
    end
    return
end

function instantiate_orig_constrs!(
    spform::Formulation{DwSp},
    origform::Formulation{Original},
    ::Env,
    annotations::Annotations,
    sp_ann
)
    !haskey(annotations.constrs_per_ann, sp_ann) && return
    constrs = annotations.constrs_per_ann[sp_ann]
    for (_, constr) in constrs
        cloneconstr!(origform, spform, spform, constr, DwSpPureConstr)
    end
    return
end

function create_side_vars_constrs!(
    spform::Formulation{DwSp},
    ::Formulation{Original},
    ::Env,
    ::Annotations
)
    name = "PricingSetupVar_sp_$(getuid(spform))"
    setvar!(
        spform, name, DwSpSetupVar; cost = 0.0, lb = 1.0, ub = 1.0, kind = Integ,
        is_explicit = true
    )
    return
end

function _dutyexpofbendmastvar(
    var::Variable, annotations::Annotations, origform::Formulation{Original}
)
    orig_coef = getcoefmatrix(origform)
    for (constrid, _) in @view orig_coef[:, getid(var)]
        constr_ann = annotations.ann_per_constr[constrid]
        #if coef != 0 && BD.getformulation(constr_ann) == BD.Benders  # TODO use haskey instead testing != 0
        if BD.getformulation(constr_ann) == BD.BendersSepSp
            return MasterBendFirstStageVar, true
        end
    end
    return MasterPureVar, true
end

# Master of Benders decomposition

function instantiate_orig_vars!(
    masterform::Formulation{BendersMaster},
    origform::Formulation{Original},
    annotations::Annotations,
    mast_ann
)
    !haskey(annotations.vars_per_ann, mast_ann) && return
    vars = annotations.vars_per_ann[mast_ann]
    for (_, var) in vars
        clonevar!(origform, masterform,  masterform, var, MasterPureVar, is_explicit = true)
    end
    return
end

function instantiate_orig_constrs!(
    masterform::Formulation{BendersMaster},
    origform::Formulation{Original},
    ::Env,
    annotations::Annotations,
    mast_ann
)
    !haskey(annotations.constrs_per_ann, mast_ann) && return
    constrs = annotations.constrs_per_ann[mast_ann]
    for (_, constr) in constrs
        cloneconstr!(
            origform, masterform, masterform, constr, MasterPureConstr, is_explicit = true
        )
    end
    return
end

function create_side_vars_constrs!(
    masterform::Formulation{BendersMaster},
    ::Formulation{Original},
    ::Env,
    ::Annotations
)
    for (_, spform) in get_benders_sep_sps(masterform.parent_formulation)
        nu_var = collect(values(filter(
            v -> getduty(v.first) == BendSpSlackSecondStageCostVar,
            getvars(spform)
        )))[1]

        name = "η[$(split(getname(spform, nu_var), "[")[end])"
        setvar!(
            masterform, name, MasterBendSecondStageCostVar;
            cost = 1.0,
            lb = getperenlb(spform, nu_var),
            ub = getperenub(spform, nu_var),
            kind = Continuous,
            is_explicit = true,
            id = Id{Variable}(MasterBendSecondStageCostVar, getid(nu_var), getuid(masterform))
        )
    end
    return
end

create_artificial_vars!(::Formulation{BendersMaster}, ::Env) = return

function instantiate_orig_vars!(
    spform::Formulation{BendersSp},
    origform::Formulation{Original},
    annotations::Annotations,
    sp_ann
)
    masterform = getmaster(spform)
    if haskey(annotations.vars_per_ann, sp_ann)
        vars = annotations.vars_per_ann[sp_ann]
        for (_, var) in vars
            clonevar!(origform, spform, spform, var, BendSpSepVar, cost = 0.0)
        end
    end
    mast_ann = get(annotations, masterform)
    if haskey(annotations.vars_per_ann, mast_ann)
        vars = annotations.vars_per_ann[mast_ann]
        for (id, var) in vars
            duty, _ = _dutyexpofbendmastvar(var, annotations, origform)
            if duty == MasterBendFirstStageVar
                name = "μ[$(split(getname(origform, var), "[")[end])"
                setvar!(
                    spform, name, BendSpSlackFirstStageVar;
                    cost = getcurcost(origform, var),
                    lb = getcurlb(origform, var),
                    ub = getcurub(origform, var),
                    kind = Continuous,
                    is_explicit = true,
                    id = Id{Variable}(BendSpSlackFirstStageVar, id, getuid(masterform))
                )
            end
        end
    end
    return
end

function _dutyexpofbendspconstr(constr, annotations::Annotations, origform)
    orig_coef = getcoefmatrix(origform)
    for (varid, _) in orig_coef[getid(constr), :]
        var_ann = annotations.ann_per_var[varid]
        if BD.getformulation(var_ann) == BD.Master
            return BendSpTechnologicalConstr, true
        end
    end
    return BendSpPureConstr, true
end

function instantiate_orig_constrs!(
    spform::Formulation{BendersSp},
    origform::Formulation{Original},
    ::Env,
    annotations::Annotations,
    sp_ann
)
    !haskey(annotations.constrs_per_ann, sp_ann) && return
    constrs = annotations.constrs_per_ann[sp_ann]
    for (_, constr) in constrs
        duty, explicit  = _dutyexpofbendspconstr(constr, annotations, origform)
        cloneconstr!(origform, spform, spform, constr, duty, is_explicit = explicit)
    end
    return
end

function create_side_vars_constrs!(
    spform::Formulation{BendersSp},
    origform::Formulation{Original},
    ::Env,
    ::Annotations
)
    sp_has_second_stage_cost = false
    global_costprofit_ub = 0.0
    global_costprofit_lb = 0.0
    for (varid, _) in getvars(spform)
        getduty(varid) == BendSpSepVar || continue
        orig_var = getvar(origform, varid)
        cost =  getperencost(origform, orig_var)
        if cost > 0.00001
            global_costprofit_ub += cost * getcurub(origform, orig_var)
            global_costprofit_lb += cost * getcurlb(origform, orig_var)
        elseif cost < - 0.00001
            global_costprofit_ub += cost * getcurlb(origform, orig_var)
            global_costprofit_lb += cost * getcurub(origform, orig_var)
        end
    end

    if global_costprofit_ub > 0.00001  || global_costprofit_lb < - 0.00001
        sp_has_second_stage_cost = true
    end

    if sp_has_second_stage_cost
        sp_coef = getcoefmatrix(spform)
        sp_id = getuid(spform)
        # Cost constraint
        nu = setvar!(
            spform, "ν[$sp_id]", BendSpSlackSecondStageCostVar;
            cost = 1.0,
            lb = - global_costprofit_lb,
            ub = global_costprofit_ub,
            kind = Continuous,
            is_explicit = true
        )
        setcurlb!(spform, nu, 0.0)
        setcurub!(spform, nu, Inf)

        cost = setconstr!(
            spform, "cost[$sp_id]", BendSpSecondStageCostConstr;
            rhs = 0.0,
            kind = Essential,
            sense = Greater,
            is_explicit = true
        )
        sp_coef[getid(cost), getid(nu)] = 1.0

        for (varid, _) in getvars(spform)
            getduty(varid) == BendSpSepVar || continue
            sp_coef[getid(cost), varid] = - getperencost(origform, varid)
        end
    end
    return
end

function assign_orig_vars_constrs!(
    destform::Formulation,
    origform::Formulation{Original},
    env::Env,
    annotations::Annotations,
    ann
)
    instantiate_orig_vars!(destform, origform, annotations, ann)
    instantiate_orig_constrs!(destform, origform, env, annotations, ann)
    clonecoeffs!(origform, destform)
end

_optimizerbuilder(opt::Function) = () -> UserOptimizer(opt)
_optimizerbuilder(opt::MOI.AbstractOptimizer) = () -> MoiOptimizer(opt)
_optimizerbuilder(opt::BD.AbstractCustomOptimizer) = () -> CustomOptimizer(opt)

function getoptimizerbuilders(prob::Problem, ann::BD.Annotation)
    optimizers = BD.getoptimizerbuilders(ann)
    if length(optimizers) > 0
        return map(o -> _optimizerbuilder(o), optimizers)
    end
    return [prob.default_optimizer_builder]
end

function buildformulations!(
    prob::Problem, reform::Reformulation, env::Env, annotations::Annotations, parent,
    node::BD.Root
)
    ann = BD.annotation(node)
    form_type = BD.getformulation(ann)
    dec_type = BD.getdecomposition(ann)
    masterform = instantiatemaster!(env, prob, reform, form_type, dec_type)
    store!(annotations, masterform, ann)
    origform = get_original_formulation(prob)
     for (_, child) in BD.subproblems(node)
        buildformulations!(prob, reform, env, annotations, node, child)
    end
    assign_orig_vars_constrs!(masterform, origform, env, annotations, ann)
    create_side_vars_constrs!(masterform, origform, env, annotations)
    create_artificial_vars!(masterform, env)
    closefillmode!(getcoefmatrix(masterform))
    push_optimizer!.(Ref(masterform), getoptimizerbuilders(prob, ann))
    push_optimizer!.(Ref(origform), getoptimizerbuilders(prob, ann))
    return
end

function buildformulations!(
    prob::Problem, reform::Reformulation, env::Env, annotations::Annotations,
    parent, node::BD.Leaf
)
    ann = BD.annotation(node)
    form_type = BD.getformulation(ann)
    dec_type = BD.getdecomposition(ann)
    masterform = getmaster(reform)
    spform = instantiatesp!(env, reform, masterform, form_type, dec_type)
    store!(annotations, spform, ann)
    origform = get_original_formulation(prob)
    assign_orig_vars_constrs!(spform, origform, env, annotations, ann)
    create_side_vars_constrs!(spform, origform, env, annotations)
    closefillmode!(getcoefmatrix(spform))
    push_optimizer!.(Ref(spform), getoptimizerbuilders(prob, ann))
    return
end

function reformulate!(prob::Problem, annotations::Annotations, env::Env)
    if getcoefmatrix(prob.original_formulation).fillmode
        closefillmode!(getcoefmatrix(prob.original_formulation))
    end
    decomposition_tree = annotations.tree
    if decomposition_tree !== nothing
        root = BD.getroot(decomposition_tree)
        reform = Reformulation()
        set_reformulation!(prob, reform)
        buildformulations!(prob, reform, env, annotations, reform, root)
        relax_integrality!(getmaster(reform))
    else
        push_optimizer!(
            prob.original_formulation,
            prob.default_optimizer_builder
        )
    end
    return
end
