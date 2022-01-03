function create_funnel_graph(domain_name, problem_name; max_levels=10)
    graph = Graph()
    domain = load_domain(domain_name)
    problem = load_problem(problem_name)
    init_propositions = get_init_propositions(domain, problem)
    ex_constraints = get_external_constraints(domain, problem)
    graph.props[1] = init_propositions
    graph.num_levels = 1

    for _=1:max_levels
        expand!(graph, domain, problem, ex_constraints)
        if goal_reached(graph, domain, problem) break end 
    end
    return graph 
end


function expand!(graph, domain, problem, constraints)
    k = graph.num_levels
    if k ≥ 2
        graph.μprops[k] = get_proposition_mutexes(graph, k)
    end
    actions = get_all_actions(domain, problem)
    for act in actions 
        if is_applicable(act, graph, constraints, k)
            instantiated_action = instantiate_action(act, graph, constraints, k)
            propositions = get_action_propositions(instantiated_action)
            push!(graph.acts[k], instantiated_action)
            if !isempty(propositions[:discrete])
                push!(graph.props[k+1][:discrete], propositions[:discrete]...) 
            end 
            if !isempty(propositions[:continuous]) 
                push!(graph.props[k+1][:continuous], propositions[:continuous]...) 
            end
        end
    end
    for prop in graph.props[k]
        if !is_continuous_proposition(prop)
            noop = get_noop_action(prop)
            push!(graph.acts[k], noop)
            push!(graph.props[k+1][:discrete], prop)
        end
    end
    graph.μacts[k] = get_action_mutexes(graph, k)

    return graph 
end 


function get_init_propositions(domain, problem)
    init_props = Dict()
    init_props[:discrete] =  collect(initstate(domain, problem).facts)

    robot_region = Region(:robot, [Ineq(1,0,0,0), Ineq(-1,0,0,0.5), Ineq(0,1,0,0), Ineq(0,-1,0,0.5)], 0.0)
    b1_region = Region(:b1, [Ineq(1,0,0,10), Ineq(-1,0,0,11), Ineq(0,1,0,10), Ineq(0,-1,0,11)], 0.0)
    init_props[:continuous] = [robot_region, b1_region]

    return init_props
end

function get_goal_propositions(domain, problem)
    goal_props = Dict()
    goal_props[:discrete] = goalstate(domain, problem).facts

    robot_region = Region(:robot, [Ineq(1,0,0,0), Ineq(-1,0,0,0.5), Ineq(0,1,0,0), Ineq(0,-1,0,0.5)], 0.0)
    b1_region = Region(:b1, [Ineq(1,0,0,20), Ineq(-1,0,0,21), Ineq(0,1,0,20), Ineq(0,-1,0,21)], 0.0)
    goal_props[:continuous] = [robot_region, b1_region]

    return goal_props
end
 

function get_external_constraints(domain, problem)
    constraints = []
    return constraints
end


function goal_reached(graph, domain, problem)
    goals = get_goal_propositions(domain, problem)
    goal_set = goals[:discrete]
    index = graph.num_levels - 1
    props = graph.props[index]
    μprops = graph.μprops[index]
    goal_found = false 
    if issubset(goal_set, props[:discrete])
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
    if goal_found
        if !regions_intersect(goals[:continuous], props[:continuous])
            goal_found = false
        end
    end
    return goal_found
end

function get_min_max(region)
    x=[0,0]; y=[0,0]
    x[1] = region.r[1].þ
    x[2] = region.r[2].þ
    y[1] = region.r[3].þ
    y[2] = region.r[4].þ
    return x,y

end

function overlaps(x1r, y1r, x2r, y2r)
    if x1r[2] < x2r[1] || x2r[2] < x1r[1] || y1r[2] < y2r[1] || y2r[2] < y1r[1]
        return false
    end
    return true 
end


function regions_intersect(regions1, regions2) 
    for r2 in regions2 
        x2r, y2r = get_min_max(r2)
        for r1 in regions1
            if r2.name == r1.name 
                x1r, y1r = get_min_max(r1)
                if !overlaps(x1r, y1r, x2r, y2r) return false end
            end
        end
    end
    return true 
end 


function get_proposition_mutexes(graph, level)
end


function get_action_mutexes(graph, level)
end


function get_all_actions(domain, problem)
end


function is_applicable(action, graph, constraints, level)
end


function instantiate_action(action, graph, constraints, level)
end


function get_action_propositions(action)
end


function is_continuous_proposition(proposition)
end


function get_noop_action(proposition)
end
