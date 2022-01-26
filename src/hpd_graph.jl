function create_funnel_graph(domain_name, problem_name; max_levels=10)
    graph = Graph()
    domain = load_domain(domain_name)
    problem = load_problem(problem_name)
    init_propositions = get_init_propositions(domain, problem)
    ex_constraints = get_external_constraints(domain, problem) 
    graph.props[1] = init_propositions
    graph.initprops = init_propositions
    graph.goalprops = get_goal_propositions(domain, problem)
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
        graph.μprops = get_proposition_mutexes(graph, k)
    end
    action_list = []
    graph.props[k+1] = Dict()
    graph.grops[k+1][:discrete] = []
    graph.props[k+1][:continuous] = []
    graph.acts[k] = []
    ext_constraint = get_external_constraint(graph, k, constraints)
    actions = get_applicable_actions(graph, k, ext_constraint, domain, problem)
    for act in actions
        propositions = get_action_propositions(act)
        push!(graph.acts[k], act)
        if !isempty(propositions[:discrete])
            push!(graph.props[k+1][:discrete], propositions[:discrete]...)
        end
        if !isempty(propositions[:continuous])
            push!(graph.props[k+1][:continuous], propositions[:continuous]...)
        end
        push!(action_list, act)
    end
    for prop in graph.props[k][:discrete]
        noop = get_noop_action(prop)
        push!(graph.acts[k], noop)
        push!(graph.props[k+1][:discrete], prop)
    end
    for cont in graph.props[k][:continuous]
        if !(cont in graph.props[k+1][:continuous])
            push!(graph.props[k+1][:continuous], cont)
        end
    end
    graph.props[k+1][:discrete] = collect(Set(graph.props[k+1][:discrete]))
    graph.μacts[k] = get_action_mutexes(graph, k)
    graph.num_levels +=1
    return graph
end

function get_init_propositions(domain::HPD.Parser.GenericDomain, 
                              problem::HPD.Parser.GenericProblem)
    init_props = Dict()
    init = initstate(domain, problem)
    init_props[:discrete] = collect(init.facts)
    init_props[:continuous] = [init.values]
    return init_props
end

function get_goal_propositions(domain::HPD.Parser.GenericDomain, 
                              problem::HPD.Parser.GenericDomain)
    goal_props = Dict()
    goal = goalstate(domain, problem)
    goal_props[:discrete] = collect(goal.facts)
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
        for goal_pair in collect(permutations(goals[:discrete], 2))
            if goal_pair in μprops 
                goal_found = false 
                break 
            end
        end
    elseif index > 1 && graph.props[index-1] == props 
        graph.leveled = true
    end
    if goal_found
        for cont in props[:continuous]
            if !in_region(goals[:continuous], cont)
                goal_found = false
            end
        end
    end
    return goal_found
end

function in_region(subregion::Dict{Any}, region::Dict{Any}) 
    for k in keys(region)
        if isa(subregion[k], Float64) && isa(region[k], Float64)
            if subregion[k] != region[k] return false end 
        elseif isa(subregion[k], Float64) && isa(region[k], Array)
            if !(region[k][1] < subregion[k] < region[k][2]) return false end 
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

function get_action_causal_fluents(act::HPD.Parser.GenericAction, 
    graph::Graph, domain::HPD.Parser.GenericDomain, 
    problem::HPD.Parser.GenericDomain)
    obs = collect(HPD.Parser.get_objects(problem))
    ob = obs[graph.indexes[2]]
    pos_prec, neg_prec = get_preconditions(domain, act, [obs])
    pos_eff, neg_eff = get_effects(domain, act, [obs])
    return pos_prec, neg_prec, pos_eff, neg_eff, [ob] 
end

function get_expression_from_string(string_exp, var) 
    string_exp = replace(string_exp, "xo"=>string(var.name))
    exp = Meta.parse(string_exp)
    return exp
end

# make more general
# only works with one object
function satisfies_precondition(action, vars, cont_props)
    if action.name == "move" || action.name == "move_holding"  return true end
    var = vars[1]
    string_exp = action.cont_precond.args[1].name
    string_exp = replace(string_exp, "xo"=>string(var.name))
    precond_exp = Meta.parse(string_exp)
    # τ = precond_exp.args[3]
    for cont_prop in cont_props 
        xr, yr = cont_prop[:xr], cont_prop[:yr]
        box_vars = ["x"*string(var.name), "y"*string(var.name)]
        for xb = box_vars
            @eval $xb = cont_prob[Symbol(xb)]
        end
        if eval(precond_exp) return true end
    end
    return false
end

function get_satisfying_index(action, var, cont_props)
    for (i,cont_prop) in enumerate(cont_props)
        xr, yr = cont_prop[:xr], cont_prop[:yr]
        box_vars = ["x"*string(var.name), "y"*string(var.name)]
        for xb = box_vars
            @eval $xb = cont_prob[Symbol(xb)]
        end
        if eval(precond_exp) return i end
    end
    return -1
end

function instantiate_action(action::HPD.Parser.GenericAction, graph::Graph,
            domain::HPD.Parser.GenericDomain, problem::HPD.Parser.GenericDomain)
    act = Funnel(action.name)
    d = 1.0
    act.pos_prec, act.neg_prec, act.pos_eff, act.neg_eff, act.params = 
                    get_action_causal_fluents(action, graph, domain, problem)
    if act.name == :pick || act.name == :place act.is_continuous = false end 
    push!(act.continuous_prec, 
      get_expression_from_string(action.cont_precond.args[1].name, act.params))
    cont_props = graph.props[graph.num_levels][:continuous]
    if action.dynamics != true
        ind = get_satisfying_index(action, act.params[1], cont_props)
        xr = cont_props[ind][:xr]; yr = cont_props[ind][:yr]
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
    end
    if action.cont_effect != true
        # for ce in action.cont_effect.args
        # cexp = Meta.parse(ce.name)
        act.end_region[Symbol("x"*string(act.params[1]))] = [xrmin, xrmax]
        act.end_region[Symbol("y"*string(act.params[1]))] = [yrmin, yrmax]
    end
    return act
    
end

function get_applicable_actions(graph::Graph, level::Int, constraint, 
        domain::HPD.Parser.GenericDomain,  problem::HPD.Parser.GenericDomain)
    applicable_actions = []
    actions = collect(values(HPD.Parser.get_actions(domain)))
    props = graph.props[level]
    for act in actions 
        fluents = get_action_causal_fluents(act, graph, domain, problem)
        pos_prec, neg_prec, pos_eff, neg_eff, vs = fluents
        if issubset(pos_prec, props[:discrete])
            if satisfies_precondition(act, vs, props[:continuous])
                instantiated_action = instantiate_action(act, graph, domain, problem)
                push!(applicable_actions, instantiated_action)
            end
        end
    end
    return applicable_actions
end

function get_proposition_mutexes(graph::Graph, level::Int)
end

function get_action_mutexes(graph::Graph, level::Int)
end

function get_action_propositions(action::Funnel)
end

function get_noop_action(proposition)
end

function get_external_constraint(graph::Graph, level::Int, constraints)
    
end