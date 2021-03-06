import Base: == 
==(f1::Funnel, f2::Funnel) = f1.name == f2.name && f1.end_region[:yr] == f2.end_region[:yr] && f1.end_region[:xr] == f2.end_region[:xr]
function create_funnel_graph(domain_name::String, problem_name::String; 
    max_levels=10, has_placement_constraint=true)
    graph = Graph()
    domain = load_domain(domain_name)
    problem = load_problem(problem_name)
    init_propositions = get_init_propositions(domain, problem)
    ex_constraints = get_external_constraints(domain, problem) 
    graph.props[1] = init_propositions
    graph.initprops = init_propositions
    graph.goalprops = get_goal_propositions(domain, problem)
    graph.has_placement_constraint = has_placement_constraint
    generate_safe_poses!(graph)
    graph.num_levels = 1

    for _=1:max_levels
        expand!(graph, domain, problem, ex_constraints)
        if goal_reached(graph, domain, problem) break end 
    end
    return graph 
end

function expand!(graph::Graph, domain::HPD.Parser.GenericDomain, 
    problem::HPD.Parser.GenericProblem, constraints)
    k = graph.num_levels  
    if k ≥ 2 
        graph.μprops[k] = get_proposition_mutexes(graph, k)
    end
    action_list = []
    graph.props[k+1] = Dict()
    graph.props[k+1][:discrete] = []
    graph.props[k+1][:continuous] = []
    graph.acts[k] = [] 
    actions = get_applicable_actions(graph, k, constraints, domain, problem)
    repped = false
    for act in actions
        propositions = get_action_propositions(act)
        push!(graph.acts[k], act)
        if !isempty(propositions[:discrete])
            push!(graph.props[k+1][:discrete], propositions[:discrete]...)
        end
        if !isempty(propositions[:continuous])
            if k>1 && !action_in_prev(graph, act, k) && (act.name == :pick || act.name == :place || act.name == :throw) 
                graph.props[k+1][:continuous]=propositions[:continuous]
                repped = true
            elseif !repped && !(propositions[:continuous][1] in graph.props[k][:continuous]) 
                push!(graph.props[k+1][:continuous], propositions[:continuous]...)
            end
        end
        push!(action_list, act)
    end
    graph.acts[k] = collect(Set(graph.acts[k])) 
    for prop in graph.props[k][:discrete]
        noop = get_noop_action(prop)
        push!(graph.acts[k], noop)
        push!(graph.props[k+1][:discrete], prop)
    end 
    graph.props[k+1][:discrete] = collect(Set(graph.props[k+1][:discrete]))
    graph.props[k+1][:continuous] = collect(Set(graph.props[k+1][:continuous]))
    graph.μacts[k] = get_action_mutexes(graph, k)
    graph.num_levels +=1 
    return graph
end

function action_in_prev(graph::Graph, action::Funnel, level::Int)
    for act in graph.acts[level-1]
        if act.name == action.name 
            return true 
        end
    end
    return false
end

function terminal_props(domain, problem; init=true)
    termprops=[]
    if init props = collect(initstate(domain, problem).facts) else props = collect(goalstate(domain, problem).facts) end
    for prop in props
        objs = prop.args
        push!(termprops, fill_proposition(prop, objs))
    end
    return termprops

end

function get_init_propositions(domain::HPD.Parser.GenericDomain, 
                              problem::HPD.Parser.GenericProblem)
    init_props = Dict()
    init = initstate(domain, problem)
    init_props[:discrete] = terminal_props(domain, problem) 
    init_props[:continuous] = [init.values]
    return init_props
end

function get_goal_propositions(domain::HPD.Parser.GenericDomain, 
                              problem::HPD.Parser.GenericProblem)
    goal_props = Dict()
    goal = goalstate(domain, problem)
    goal_props[:discrete] = terminal_props(domain, problem;init=false) 
    goal_props[:continuous] = goal.values
    return goal_props
end

function goal_reached(graph::Graph, domain::HPD.Parser.GenericDomain, 
                                    problem::HPD.Parser.GenericProblem)
    goals = get_goal_propositions(domain, problem)
    index = graph.num_levels-1
    props = graph.props[index]
    μprops = graph.μprops[index]
    goal_found = false
    if issubset(goals[:discrete], props[:discrete])
        goal_found = true  
    elseif index > 1 && graph.props[index-1] == props 
        graph.leveled = true
    end
    if goal_found
        goal_found = false
        for cont in props[:continuous]
            if in_region(goals[:continuous], cont)
                goal_found = true
                break
            end
        end
    end
    return goal_found
end

function in_region(subregion , region )  
    for k in keys(subregion)
        if isa(subregion[k], Float64) && isa(region[k], Float64)
            if subregion[k] != region[k] return false end 
        elseif isa(subregion[k], Float64) && isa(region[k], Array)  
            if !(region[k][1] < subregion[k] < region[k][2]) return false; end 
        elseif isa(subregion[k], Array) && isa(region[k], Array)
            if subregion[k][2] < region[k][1] || region[k][2] < subregion[k][1]
                return false
            end
        elseif isa(subregion[k], Array) && isa(region[k], Float64)
            if !(subregion[k][1] < region[k] < subregion[k][2]) return false end
        end
    end 
    return true            
end

function get_external_constraints(domain::HPD.Parser.GenericDomain, 
    problem::HPD.Parser.GenericProblem) 
    ext_consts = [Meta.parse(c.name) for c in problem.external_constraints.args]
    return ext_consts
end

function fill_proposition(proposition, objs) 
    if isempty(objs) objs=Term[] end
    prop = Compound(Symbol(proposition.name), objs)
    return prop
end

function get_preconditions(domain, action, vars) 
    vardict = Dict()
    [vardict[action.args[i]] = vars[i] for i=1:length(action.args)]
    pos, neg = [], []
    for pre in action.precond.args 
        if pre.name == :not 
            vs = []
            [push!(vs, vardict[arg]) for arg in pre.args[1].args]
            prop = fill_proposition(pre.args[1], vs)
            push!(neg, prop)
        else 
            vs = []
            [push!(vs, vardict[arg]) for arg in pre.args]
            prop = fill_proposition(pre, vs)
            push!(pos, prop)
        end
    end
    return pos, neg
end

function get_effects(domain, action, vars) 
    vardict = Dict()
    [vardict[action.args[i]] = vars[i] for i=1:length(action.args)]
    pos, neg = [], []
    for eff in action.effect.args 
        if eff.name == :not 
            vs = []
            [push!(vs, vardict[arg]) for arg in eff.args[1].args]
            prop = fill_proposition(eff.args[1], vs)
            push!(neg, prop)
        else 
            vs = []
            [push!(vs, vardict[arg]) for arg in eff.args]
            prop = fill_proposition(eff, vs)
            push!(pos, prop)
        end
    end
    return pos, neg
end

function get_action_causal_fluents(act::HPD.Parser.GenericAction, 
    graph::Graph, domain::HPD.Parser.GenericDomain, 
    problem::HPD.Parser.GenericProblem)
    obs = collect(HPD.Parser.get_objects(problem))
    ob = obs[graph.indexes[1]]
    pos_prec, neg_prec = get_preconditions(domain, act, [ob])
    pos_eff, neg_eff = get_effects(domain, act, [ob])
    return pos_prec, neg_prec, pos_eff, neg_eff, [ob] 
end
 
function satisfies_precondition(action, vars, cont_prop)
    if action.name == :move || action.name == :move_holding  return true end
    var = vars[1]
    string_exp = action.cont_precond.args[1].name 
    precond_exp = Meta.parse(string_exp) 
    if action.name == :pick 
        xr, yr = cont_prop[:xr], cont_prop[:yr]
        box_vars = ["x"*string(var.name), "y"*string(var.name)] 
        xo, yo = cont_prop[Symbol(box_vars[1])], cont_prop[Symbol(box_vars[2])] 
        if isa(xr, Array) && isa(xo, Float64) 
            return xr[1]<=xo<=xr[2] && yr[1]<=yo<=yr[2]
        elseif isa(xr, Float64) && isa(xo, Array) 
            return xo[1]<=xr<=xo[2] && yo[1]<=yr<=yo[2] 
        elseif isa(xr, Array) && isa(xo, Array) 
            return xr[1]<=xo[1] && xo[2]<=xr[2] && yr[1]<=yo[1] && yo[2]<=yr[2]             
        else 
            fun = @eval (xr, yr, xo, yo) -> $precond_exp 
        end
        if Base.invokelatest(fun, xr, yr, xo, yo) return true end 
    elseif action.name == :place || action.name == :throw
        xr, yr = cont_prop[:xr], cont_prop[:yr]
        xg, yg = cont_prop[:xg], cont_prop[:yg]
        if isa(xr, Array) func(xr,yr,xg, yg) = xr[1]<=xg<=xr[2] && yr[1]<=yg<=yr[2]
        else func = @eval (xr, yr, xg, yg) -> $precond_exp end
        if Base.invokelatest(func, xr, yr, xg, yg) return true end
    end
    return false
end

function get_satisfying_index(action, var, cont_props)
    if action.name == :move || action.name == :move_holding return 1 end
    string_exp = action.cont_precond.args[1].name
    precond_exp = Meta.parse(string_exp) 
    for (i,cont_prop) in enumerate(cont_props)
        xr, yr = cont_prop[:xr], cont_prop[:yr] 
        box_vars = ["x"*string(var.name), "y"*string(var.name)]
        xo, yo = cont_prop[Symbol(box_vars[1])], cont_prop[Symbol(box_vars[2])]
        if isa(xr, Array) && isa(xo, Float64) 
            if xr[1]<=xo<=xr[2] && yr[1]<=yo<=yr[2] return i end 
        elseif isa(xr, Array) && isa(xo, Array)
            if xr[1]<=xo[1] && xo[2]<=xr[2] && yr[1]<=yo[1] && yo[2]<=yr[2] return i end 
        else 
            fun = @eval (xr, yr, xg, yg) -> $precond_exp  
            if Base.invokelatest(fun, xr, yr, xo, yo) return i end 
        end
    end
    return -1
end

function instantiate_action!(action::HPD.Parser.GenericAction, graph::Graph,
            domain::HPD.Parser.GenericDomain, problem::HPD.Parser.GenericProblem, ind::Int, contprop)
    funnel = Funnel(action.name)
    d = 1.0
    funnel.pos_prec, funnel.neg_prec, funnel.pos_eff, funnel.neg_eff, var = 
                    get_action_causal_fluents(action, graph, domain, problem)
    if funnel.name == :pick || funnel.name == :place || funnel.name == :throw funnel.is_continuous = false end 
    if funnel.name == :move funnel.params = [] else funnel.params = var end
    if length(action.cont_precond.args)>0 
        push!(funnel.continuous_prec, Meta.parse(action.cont_precond.args[1].name)) 
    end 
    funnel.cont_ind = ind 
    if action.dynamics != true 
        xr = contprop[:xr]; yr = contprop[:yr]
        vxmin = graph.props[1][:continuous][1][:vxmin]
        vxmax = graph.props[1][:continuous][1][:vxmax]
        vymin = graph.props[1][:continuous][1][:vymin]
        vymax = graph.props[1][:continuous][1][:vymax]
        if isa(xr, Array)
            xrmin = xr[1] + d*vxmin
            xrmax = xr[2] + d*vxmax 
            yrmin = yr[1] + d*vymin 
            yrmax = yr[2] + d*vymax 
        else
            xrmin = xr + d*vxmin 
            xrmax = xr + d*vxmax 
            yrmin = yr + d*vymin 
            yrmax = yr + d*vymax
        end
        funnel.end_region[:xr] = [xrmin, xrmax]
        funnel.end_region[:yr] = [yrmin, yrmax]
    else
        funnel.end_region[:xr] = contprop[:xr]
        funnel.end_region[:yr] = contprop[:yr]
    end
    #propagate other continuous values
    for key in keys(contprop) 
        if !(key in keys(funnel.end_region))
            funnel.end_region[key] = contprop[key]
        end
    end
    if action.cont_effect != true 
        if action.name == :pick 
            funnel.end_region[:xr] = contprop[Symbol("x"*string(var[1]))]   
            funnel.end_region[:yr] = contprop[Symbol("y"*string(var[1]))] 
        elseif action.name == :place || action.name == :throw
            funnel.end_region[:xr] = contprop[:xg]  
            funnel.end_region[:yr] = contprop[:yg]    
            funnel.end_region[Symbol("x"*string(var[1]))] = contprop[:xg]   
            funnel.end_region[Symbol("y"*string(var[1]))] = contprop[:yg]  
        end  
        if action.name == :move_holding 
            funnel.end_region[Symbol("x"*string(var[1]))] = contprop[:xr]
            funnel.end_region[Symbol("y"*string(var[1]))] = contprop[:yr] 
        end  
    end
    return funnel
end

function working_instantiate_action!(action::HPD.Parser.GenericAction, graph::Graph,
    domain::HPD.Parser.GenericDomain, problem::HPD.Parser.GenericProblem, ind::Int, contprop)
    funnel = Funnel(action.name)
    d = 1.0
    funnel.pos_prec, funnel.neg_prec, funnel.pos_eff, funnel.neg_eff, var = 
                get_action_causal_fluents(action, graph, domain, problem)
    if funnel.name == :pick || funnel.name == :place  || funnel.name == :throw funnel.is_continuous = false end 
    if funnel.name == :move funnel.params = [] else funnel.params = var end
    if length(action.cont_precond.args)>0 
        push!(funnel.continuous_prec, Meta.parse(action.cont_precond.args[1].name)) 
    end
    cont_props = graph.props[graph.num_levels][:continuous] 
    if action.dynamics != true 
        xr = contprop[:xr]; yr = contprop[:yr]
        vxmin = graph.props[1][:continuous][1][:vxmin]
        vxmax = graph.props[1][:continuous][1][:vxmax]
        vymin = graph.props[1][:continuous][1][:vymin]
        vymax = graph.props[1][:continuous][1][:vymax]
        if isa(xr, Array)
            xrmin = xr[1] + d*vxmin
            xrmax = xr[2] + d*vxmax 
            yrmin = yr[1] + d*vymin 
            yrmax = yr[2] + d*vymax 
        else
            xrmin = xr + d*vxmin 
            xrmax = xr + d*vxmax 
            yrmin = yr + d*vymin 
            yrmax = yr + d*vymax
        end
        funnel.end_region[:xr] = [xrmin, xrmax]
        funnel.end_region[:yr] = [yrmin, yrmax]
    else
        funnel.end_region[:xr] = contprop[:xr]
        funnel.end_region[:yr] = contprop[:yr]
    end
    #propagate other continuous values
    for key in keys(contprop) 
        if !(key in keys(funnel.end_region))
            funnel.end_region[key] = contprop[key]
        end
    end
    if action.cont_effect != true 
        if action.name == :move_holding 
            funnel.end_region[Symbol("x"*string(var[1]))] = contprop[:xr]
            funnel.end_region[Symbol("y"*string(var[1]))] = contprop[:yr] 
        end  
    end
    return funnel
end

function get_applicable_actions(graph::Graph, level::Int, constraints::Vector{Expr}, domain::HPD.Parser.GenericDomain,  problem::HPD.Parser.GenericProblem)
    applicable_actions = []
    actions = collect(values(HPD.Parser.get_actions(domain))) 
    props = graph.props[level]
    num_obs = length(problem.objects)
    update_indexes!(graph, level, num_obs)
    for act in actions 
        fluents = get_action_causal_fluents(act, graph, domain, problem)
        pos_prec, neg_prec, pos_eff, neg_eff, vs = fluents  
        if issubset(pos_prec, props[:discrete]) 
            for (ind,contprop) in enumerate(props[:continuous])
                pose = get_placement_pose(graph, constraints, vs) 
                contprop[:xg] = pose[1]; contprop[:yg] = pose[2]
                if satisfies_precondition(act, vs, contprop)  
                    instantiated_action = instantiate_action!(act, graph, domain, problem,ind,contprop)
                    if !(instantiated_action in applicable_actions) 
                            push!(applicable_actions, instantiated_action) 
                    end 
                end
            end
        end
    end  
    applicable_actions = prune_funnels(applicable_actions) 
    return applicable_actions
end

function int_range(funnel::Union{Funnel, Nothing})
    if isnothing(funnel) return -99 end
    array = funnel.end_region[:yr]
    return array[2]-array[1]
end

function prune_funnels(applicable_actions::Vector{Any})
    pruned_actions = []
    fun_dicts = Dict(:move=>[], :move_holding=>[], :pick=>[], :place=>[], :throw=>[])
    for action in applicable_actions push!(fun_dicts[action.name], action) end
    for key in keys(fun_dicts)
        if !isempty(fun_dicts[key])
            max_fun = nothing
            dim1 = []
            for act in fun_dicts[key]
                if isa(act.end_region[:yr], Array) && (act.name == :pick || act.name == :place || act.name == :throw)
                    continue 
                end
                if isa(act.end_region[:yr], Array) 
                    if int_range(act) > int_range(max_fun) max_fun = act  end
                else
                    push!(dim1, act)
                end
            end
            fun_dicts[key] = [max_fun, dim1...]
        end
    end
    pruned_actions = reduce(vcat, collect(values(fun_dicts)))
    filter!(x->x !== nothing, pruned_actions)
    return pruned_actions
end

function update_indexes!(graph::Graph, level::Int, num_obs::Int)
    props = graph.props[level]
    if graph.has_placement_constraint donesymb = :laid else donesymb = :stocked end
    laid_fluent = Compound(donesymb, [Const(Symbol("b"*string(graph.indexes[1])))])
    if laid_fluent in props[:discrete] && graph.indexes[1] < num_obs
        graph.indexes[1]+=1
    end
end

function get_proposition_mutexes(graph::Graph, level::Int)
    disc_props = graph.props[level][:discrete]
    μprops=[] 
    action_list = graph.acts[level-1]
    for a in action_list
        for (p,q) in collect(permutations(disc_props, 2)) 
            if p in a.pos_eff && q in a.neg_eff 
                push!(μprops, [p,q])
            end 
        end
        
    end  
    return μprops
end

function action_pairs_are_mutex(a, b, μprops) 
    if !isempty(intersect(a.neg_eff, union(b.pos_prec, b.pos_eff))) ||
        !isempty(intersect(b.neg_eff, union(a.pos_prec, a.pos_eff)))
         return true
    end
    if !isempty(μprops)
        for μ in μprops 
            p = μ[1]
            q = μ[2] 
            if (p in a.pos_prec && q in b.pos_prec)
                return true
            end 
        end
    end
    return false
end

function get_action_mutexes(graph::Graph, level::Int)
    actions = graph.acts[level]
    μacts = []
    μprops = graph.μprops[level]
    for (a, b) in collect(permutations(actions, 2)) 
        if action_pairs_are_mutex(a, b, μprops)
            push!(μacts, [a, b])
        end 
    end
    return μacts
end

function get_action_propositions(action::Funnel)
    propositions = Dict()
    propositions[:discrete] = action.pos_eff
    propositions[:continuous] = [action.end_region]
    return propositions
end

function get_noop_action(proposition)
    noop = Funnel(:noop)
    noop.pos_prec = [proposition]
    noop.pos_eff = [proposition]
    noop.is_continuous = false
    return noop
end

function get_external_constraint(graph::Graph, level::Int, constraints)
    return constraints #[graph.indexes[1]]
end
 
function generate_safe_poses!(graph::Graph)
    xmin = graph.initprops[:continuous][1][:xbmin]
    ymin = graph.initprops[:continuous][1][:ybmin]
    xmax = graph.initprops[:continuous][1][:xbmax]
    ymax = graph.initprops[:continuous][1][:ybmax]
    Δ = graph.initprops[:continuous][1][:bdelta]
    xs = xmin:Δ:xmax
    ys = ymin:Δ:ymax
    graph.safe_poses = [xs, ys]
end

function get_placement_pose(graph::Graph, constraints, vs)
    if graph.has_placement_constraint 
        ls = get_poses_from_constraints(constraints, graph)
        return ls[graph.indexes[1]] 
    else
        ob = vs[1].name 
        kx = Symbol("x"*string(ob))
        ky = Symbol("y"*string(ob))
        goal_pose = [graph.goalprops[:continuous][kx], 
                graph.goalprops[:continuous][ky]]
        return goal_pose
    end
end 

function get_pose_from_fxn!(exp::Expr, ls::Vector{Any}, graph::Graph) 
    xs = graph.safe_poses[1]
    ys = graph.safe_poses[2]
    fun = @eval (xg, yg) -> $exp
    for x in xs 
        for y in ys  
            if Base.invokelatest(fun, x, y) 
                push!(ls, (x, y))
            end
        end
    end  
end

function get_poses_from_constraints(exps::Vector{Expr}, graph::Graph)
    ls = []
    for exp in exps 
        get_pose_from_fxn!(exp, ls, graph) 
    end
    return collect(Set(ls)) 
end 