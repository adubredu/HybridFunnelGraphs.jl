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
    graph.props[k+1] = Dict()
    graph.grops[k+1][:discrete] = []
    graph.props[k+1][:continuous] = Dict()
    graph.acts[k] = []
    ext_constraint = get_external_constraint(graph, k, constraints)
    actions = get_applicable_actions(graph, k, ext_constraint)
    for act in actions 
        # get act propositions 


end

function get_init_propositions(domain::HPD.Parser.GenericDomain, 
                              problem::HPD.Parser.GenericProblem)
    init_props = Dict()
    init = initstate(domain, problem)
    init_props[:discrete] = collect(init.facts)
    init_props[:continuous] = init.values
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
        if !in_region(goals[:continuous], props[:continuous])
            goal_found = false
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

function get_external_constraint(graph, level, constraints)
end

function get_applicable_actions(graph, level, constraint)
end

function get_proposition_mutexes(graph, level)
end