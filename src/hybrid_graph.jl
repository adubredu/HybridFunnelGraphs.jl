#TODO: Adds continuous constraints
function add_continuous_constraints!(hybrid_action)
    name = hybrid_action.name 
    if name == :pick 

    elseif name == :place 

    elseif name == :move 

    end
end


function get_all_actions(domain, problem)
    actions = []
    obs = collect(PDDL.get_objects(problem))
    for act in values(PDDL.get_actions(domain))
        vars = act.args   
        for vs in collect(permutations(obs, length(vars)))
            a = HybridAction(act.name)
            a.pos_prec, a.neg_prec = get_preconditions(domain, act, vs)
            a.pos_eff, a.neg_eff = get_effects(domain, act, vs) 
            a.params = vs
            add_continuous_constraints!(a)
            push!(actions, a)
        end
    end
    return actions
end


function get_preconditions(domain, act, vars)
    action = PDDL.get_action(domain, act.name)
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


function get_effects(domain, act, vars)
    action = PDDL.get_action(domain, act.name)
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


function fill_proposition(proposition, objs) 
    if isempty(objs) objs=Term[] end
    prop = Compound(Symbol(proposition.name), objs)
    return prop
end

#TODO: Augment with continuous constraint checking
function goal_reached!(domain, problem, graph)
    goal_set = goalstate(domain, problem).facts
    goal_set = [goal_set...]
    index = graph.num_levels - 1
    props = graph.props[index]
    μprops = graph.μprops[index]
    goal_found = false 
    if issubset(goal_set, props)
        goal_found = true
        for goal_pair in collect(permutations(goal_set, 2))
            if goal_pair in μprops 
                goal_found = false
                break 
            end
        end
    elseif index > 0 && graph.props[index-1] == props 
        graph.leveled = true
    end 
    return goal_found
end

#TODO: Add continuous mutex checking
function is_mutex_acts(act_pair, μprops)
    a = act_pair[1]
    b = act_pair[2]

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
 
#TODO: Add continuous mutex checking
function is_mutex_props(prop_pair, action_list, μacts)
    p = prop_pair[1]
    q = prop_pair[2]

    for a in action_list
        if p in a.pos_eff && q in a.pos_eff 
            return false
        end
    end

    actions_with_p = Set()
    for a in action_list 
        if p in a.pos_eff 
            push!(actions_with_p, a)
        end
    end

    actions_with_q = Set()
    for a in action_list 
        if q in a.pos_eff 
            push!(actions_with_q, a)
        end
    end

    μall = true 
    for p_action in actions_with_p
        for q_action in actions_with_q 
            if p_action == q_action return false end 
            if !([p_action, q_action] in μacts)
                μall = false
                break 
            end
        end
        if !μall break end 
    end
    return μall 
end

#TODO: Add continuous prop expansion
function expand!(domain, problem, graph)
    level = graph.num_levels 
    #As 
    action_list = []
    for action in get_all_actions(domain, problem)
        if action_is_applicable(action, graph.props[level-1], graph.μprops[level-1]) 
            push!(action_list, action)
        end
    end
    for prop in graph.props[level-1]
        push!(action_list, HybridNoOp(prop))
    end
    graph.acts[level] = action_list

    #Ps 
    proposition_list = Set()
    for action in action_list 
        for eff in action.pos_eff 
            push!(proposition_list, eff)
        end
    end
    graph.props[level] = collect(proposition_list)
    proposition_list = collect(proposition_list)

    #μA
    action_μ_list = []
    for act_pair in collect(permutations(action_list, 2))
        if is_mutex_acts(act_pair, graph.μprops[level-1])
            push!(action_μ_list, act_pair)
        end
    end
    graph.μacts[level] = action_μ_list

    #μP 
    proposition_μ_list = []
    for prop_pair in collect(permutations(proposition_list, 2))
        if is_mutex_props(prop_pair, action_list,action_μ_list)
            if !(prop_pair in proposition_μ_list)
                swapped = [prop_pair[2], prop_pair[1]]
                if !(swapped in proposition_μ_list)
                    push!(proposition_μ_list, prop_pair)
                end
            end
        end
    end
    graph.μprops[level] = proposition_μ_list
    
    graph.num_levels = level + 1 
    if graph.props[level-1] == proposition_list
        graph.leveled = true
    end

    return graph 
end

#TODO: add continuous init state
function get_init_propositions(domain, problem)
    initprops=[]
    inits = collect(initstate(domain, problem).facts)
    for init in inits
        objs = init.args
        push!(initprops, fill_proposition(init, objs))
    end
    return initprops

end


#TODO: add continuous goal state
function get_goal_propositions(domain, problem)
    goalprops=[]
    goals = collect(goalstate(domain, problem).facts)
    for goal in goals
        objs = goal.args
        push!(goalprops, fill_proposition(goal, objs))
    end
    return goalprops

end


function create_graph(domain, problem; max_levels=10)
    graph = Graph()
    graph.num_levels = 1 
    graph.props[0] =  get_init_propositions(domain, problem)

    for _ in 1:max_levels
        expand!(domain, problem, graph)
        if goal_reached!(domain, problem, graph) break end 
        if graph.leveled break end
    end
    graph.initprops = get_init_propositions(domain, problem)
    graph.goalprops = get_goal_propositions(domain, problem) 
    return graph  
end