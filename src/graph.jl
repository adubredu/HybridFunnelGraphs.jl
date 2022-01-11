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


function expand!(graph, domain, problem, constraints)
    k = graph.num_levels
    if k ≥ 2
        graph.μprops[k] = get_proposition_mutexes(graph, k)
    end 
    actions = get_all_actions(domain, problem)
    action_list=[]
    graph.props[k+1] = Dict()
    graph.props[k+1][:discrete]=[]
    graph.props[k+1][:continuous]=[]
    graph.acts[k]=[]
    for act in actions 
        if is_applicable(act, graph, constraints, k)
            instantiated_action = compute_funnel(act, graph, constraints, k)
            propositions = get_action_propositions(instantiated_action)
            push!(graph.acts[k], instantiated_action)
            if !isempty(propositions[:discrete])
                push!(graph.props[k+1][:discrete], propositions[:discrete]...) 
            end 
            if !isempty(propositions[:continuous]) 
                push!(graph.props[k+1][:continuous], propositions[:continuous]...) 
            end
            push!(action_list, instantiated_action)
        end
    end
    for prop in graph.props[k][:discrete] 
            noop = get_noop_action(prop)
            push!(graph.acts[k], noop)
            push!(graph.props[k+1][:discrete], prop) 
    end
    for region in graph.props[k][:continuous] 
        if !(region in graph.props[k+1][:continuous]) && region.name == :b1
            mn = get_maintain_funnel(region) 
            push!(graph.acts[k], mn) 
            push!(graph.props[k+1][:continuous], region) 
        end
    end
    graph.props[k+1][:discrete] = collect(Set(graph.props[k+1][:discrete]))
    graph.props[k+1][:continuous] = get_unique_regions(graph.props[k+1][:continuous])
    graph.μacts[k] = get_action_mutexes(graph, k)
    graph.num_levels+=1
    # println(" ")
    return graph 
end 

function get_unique_regions(regions)
    unique_regions = []
    for reg in regions 
        if !any([reg.r[1].ϕ₁ == x.r[1].ϕ₁ && reg.r[1].ϕ₂ == x.r[1].ϕ₂ && reg.r[1].ϕ₃ == x.r[1].ϕ₃ && reg.r[1].þ == x.r[1].þ && 
            reg.r[2].ϕ₁ == x.r[2].ϕ₁ && reg.r[2].ϕ₂ == x.r[2].ϕ₂ && reg.r[2].ϕ₃ == x.r[2].ϕ₃ && reg.r[2].þ == x.r[2].þ && 
            reg.r[3].ϕ₁ == x.r[3].ϕ₁ && reg.r[3].ϕ₂ == x.r[3].ϕ₂ && reg.r[3].ϕ₃ == x.r[3].ϕ₃ && reg.r[3].þ == x.r[3].þ &&  
            reg.r[4].ϕ₁ == x.r[4].ϕ₁ && reg.r[4].ϕ₂ == x.r[4].ϕ₂ && reg.r[4].ϕ₃ == x.r[4].ϕ₃ && reg.r[4].þ == x.r[4].þ  
            for x in unique_regions])
            push!(unique_regions, reg)
        end
    end
    return unique_regions
end


function fill_proposition(proposition, objs) 
    if isempty(objs) objs=Term[] end
    prop = Compound(Symbol(proposition.name), objs)
    return prop
end



function term_props(domain, problem; init=true)
    termprops=[]
    if init props = collect(initstate(domain, problem).facts) else props = collect(goalstate(domain, problem).facts) end
    for prop in props
        objs = prop.args
        push!(termprops, fill_proposition(prop, objs))
    end
    return termprops

end


function get_init_propositions(domain, problem)
    init_props = Dict()
    init_props[:discrete] =  term_props(domain, problem)

    robot_region = Region(:robot, [Ineq(1,0,0,0), Ineq(-1,0,0,0.5), Ineq(0,1,0,0), Ineq(0,-1,0,0.5)], 0.0)
    b1_region = Region(:b1, [Ineq(1,0,0,10), Ineq(-1,0,0,11), Ineq(0,1,0,10), Ineq(0,-1,0,11)], 0.0)
    init_props[:continuous] = [robot_region, b1_region]
    # init_props[:robot_position] = (0.25,0.25,0.0)
    # init_props[:b1_position] = [10.5, 10.5, 0.0]

    return init_props
end


function get_goal_propositions(domain, problem)
    goal_props = Dict()
    goal_props[:discrete] = term_props(domain, problem;init=false)

    # robot_region = Region(:robot, [Ineq(1,0,0,0), Ineq(-1,0,0,0.5), Ineq(0,1,0,0), Ineq(0,-1,0,0.5)], 0.0)
    b1_region = Region(:b1, [Ineq(1,0,0,20), Ineq(-1,0,0,21), Ineq(0,1,0,20), Ineq(0,-1,0,21)], 0.0)
    robot_region = Region(:robot, [Ineq(1,0,0,40), Ineq(-1,0,0,41), Ineq(0,1,0,40), Ineq(0,-1,0,41)], 0.0)
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


function get_min_max(region)
    x=[0.,0.]; y=[0.,0.]
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

function intersects(reg1, reg2)
    x1r, y1r = get_min_max(reg1)
    x2r, y2r = get_min_max(reg2)
    return overlaps(x1r, y1r, x2r, y2r)
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

# function in_region(regions1, regions2) 
#     for r2 in regions2 
#         x2r, y2r = get_min_max(r2)
#         for r1 in regions1
#             if r2.name == r1.name 
#                 x1r, y1r = get_min_max(r1)
#                 if overlaps(x1r, y1r, x2r, y2r) return true end
#             end
#         end
#     end
#     return false 
# end 

function in_region(regions1, regions2)
    allbool = []
    for r1 in regions1 
        x1r, y1r = get_min_max(r1)
        isin = false 
        for r2 in regions2 
            if r1.name == r2.name  
                x2r, y2r = get_min_max(r2)
                if overlaps(x1r, y1r, x2r, y2r) isin = true; break end 
            end
        end
        push!(allbool, isin)
    end
    return all(allbool)
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


function get_all_actions(domain, problem)
    actions = []
    obs = collect(PDDL.get_objects(problem))
    for act in values(PDDL.get_actions(domain))
        vars = act.args
        # if act.name != :move
        for vs in collect(permutations(obs, length(vars)))
            a = Funnel(act.name)
            a.pos_prec, a.neg_prec = get_preconditions(domain, act, vs)
            a.pos_eff, a.neg_eff = get_effects(domain, act, vs)
            a.params = vs 
            if act.name == :pick #make more general. preferrably  from prob.pddl
                push!(a.continuous_prec, Region(:robot, [Ineq(1,0,0,20), Ineq(-1,0,0,21), Ineq(0,1,0,20), Ineq(0,-1,0,21)], 0.0))
                a.is_continuous = false

            elseif act.name == :place #make more general later
                push!(a.continuous_prec, Region(:b1, [Ineq(1,0,0,0), Ineq(-1,0,0,0.5), Ineq(0,1,0,0), Ineq(0,-1,0,0.5)], 0.0))
                push!(a.continuous_prec, Region(:robot, [Ineq(1,0,0,0), Ineq(-1,0,0,0.5), Ineq(0,1,0,0), Ineq(0,-1,0,0.5)], 0.0))
                a.is_continuous = false

            elseif act.name == :move_holding #compute funnel at instantitation
                push!(a.continuous_prec, Region(:b1, [Ineq(1,0,0,-Inf), Ineq(-1,0,0,Inf), Ineq(0,1,0,-Inf), Ineq(0,-1,0,Inf)], 0.0))
                a.dynamics = Dynamics(I(2), I(2), [-5.,5], [-5, 5], 1)
                a.is_continuous = true
            elseif act.name == :move 
                push!(a.continuous_prec, Region(:b1, [Ineq(1,0,0,-Inf), Ineq(-1,0,0,Inf), Ineq(0,1,0,-Inf), Ineq(0,-1,0,Inf)], 0.0))
                a.dynamics = Dynamics(I(2), I(2), [-5.,5], [-5, 5], 1)
                a.is_continuous = true
            end
            push!(actions, a)
        end 
    end
    return actions
end


function apply_constraints(xmin, xmax, ymin, ymax, constraints)
    isempty(constraints) && return xmin, xmax, ymin, ymax
end


#TODO: get init position ranges from problem.pddl 
function compute_funnel(action, graph, constraints, level)  
    if !action.is_continuous return action end 
    regions = graph.props[level][:continuous] 
    xyr_init = get_min_max(regions[1])
     
    a = action
    d = a.dynamics.d
    if action.name == :move 
        xmin = xyr_init[1][1] + d*a.dynamics.vx_range[1]
        xmax = xyr_init[1][2] + d*a.dynamics.vx_range[2]
        ymin = xyr_init[2][1] + d*a.dynamics.vy_range[1]
        ymax = xyr_init[2][2] + d*a.dynamics.vy_range[2]
        xmin, xmax, ymin, ymax = apply_constraints(xmin, xmax, ymin, ymax, constraints)
        push!(action.end_region, Region(:robot, [Ineq(1,0,0,xmin), Ineq(-1,0,0,xmax), Ineq(0,1,0,ymin), Ineq(0,-1,0,ymax)], 0.0))

    elseif action.name == :move_holding 
        xyo_init = get_min_max(regions[2])
        xmin = xyr_init[1][1] + d*a.dynamics.vx_range[1]
        xmax = xyr_init[1][2] + d*a.dynamics.vx_range[2]
        ymin = xyr_init[2][1] + d*a.dynamics.vy_range[1]
        ymax = xyr_init[2][2] + d*a.dynamics.vy_range[2]

        xomin = xyo_init[1][1] + d*a.dynamics.vx_range[1]
        xomax = xyo_init[1][2] + d*a.dynamics.vx_range[2]
        yomin = xyo_init[2][1] + d*a.dynamics.vy_range[1]
        yomax = xyo_init[2][2] + d*a.dynamics.vy_range[2]
        xomin, xomax, yomin, yomax = apply_constraints(xomin, xomax, yomin, yomax, constraints)

        push!(action.end_region, Region(:robot, [Ineq(1,0,0,xmin), Ineq(-1,0,0,xmax), Ineq(0,1,0,ymin), Ineq(0,-1,0,ymax)], 0.0))
        push!(action.end_region, Region(:b1, [Ineq(1,0,0,xomin), Ineq(-1,0,0,xomax), Ineq(0,1,0,yomin), Ineq(0,-1,0,yomax)], 0.0))
    end 
    return action

end


function get_noop_action(proposition)
    noop = Funnel(:noop)
    noop.pos_prec = [proposition]
    noop.pos_eff = [proposition]
    noop.is_continuous = false
    return noop
end

function get_maintain_funnel(region)
    maintain = Funnel(:maintain)
    maintain.continuous_prec=[region]
    maintain.end_region=[region]
    return maintain     
end


function is_applicable(action, graph, constraints, level) 
    props = graph.props[level][:discrete]
    μprops = graph.μprops[level]
    if issubset(action.pos_prec, props) #&& isdisjoint(action.neg_prec, props)
        app = true 
        # println("Level ",level," action ",action.name," is subset to ",props)
        if !isempty(μprops)
            for precondition in collect(permutations(action.pos_prec, 2))
                if precondition in μprops
                    app = false
                    break 
                end
            end
        end
        
    else
        # println("Level ",level," action ",action.name, " ",action.pos_prec," is NOT subset to ",props)
        app = false
    end 
    if app
        regions2 = graph.props[level][:continuous]
        regions1 = action.continuous_prec
        
        if in_region(regions1, regions2) 
            app = true 
            # println("Level ",level," action ",action.name," intersects prop region")
        else 
            # println("Level ",level," action ",action.name," DOESN'T intersect prop region")
            app = false 
        end
    end
    return app  
end

#for continuous, resolved region might have something to say about it
function get_action_propositions(action)
    propositions = Dict()
    propositions[:discrete]=[]
    propositions[:continuous]=[]
    for prop in action.pos_eff 
        push!(propositions[:discrete], prop)
        # println("pushed ",action.name,"'s ",prop)
    end
    for region in action.end_region 
        push!(propositions[:continuous], region)
    end
    return propositions 
end


function get_proposition_mutexes(graph, level)
    props = vcat(graph.props[level][:discrete], graph.props[level][:continuous])
    μprops=[]
    actions_with_p = Set()
    actions_with_q = Set()
    action_list = graph.acts[level-1]
    for a in action_list
        for (p,q) in collect(permutations(props, 2)) 
            if p in a.pos_eff && q in a.neg_eff 
                push!(μprops, [p,q])
            end
            if p in a.pos_eff push!(actions_with_p, a) end 
            if q in a.pos_eff push!(actions_with_q, a) end
        end
        
    end
    for (p,q) in collect(permutations(graph.props[level][:continuous], 2)) 
        if !regions_intersect([p],[q]) push!(μprops, [p,q]) end 
    # end
        μacts = graph.μacts[level-1]
        for p_action in actions_with_p 
            for q_action in actions_with_q
                if p_action != q_action
                    if [p_action, q_action] in μacts 
                        push!(μprops, [p,q])
                    end
                end
            end
        end
    end
    return μprops 
end


function action_pairs_are_mutex(a, b, μprops) 
    #1
    if !isempty(intersect(a.neg_eff, union(b.pos_prec, b.pos_eff))) ||
        !isempty(intersect(b.neg_eff, union(a.pos_prec, a.pos_eff)))
         return true
    end
    if !isempty(μprops)
        for μ in μprops 
            p = μ[1]
            q = μ[2]
            if typeof(p) != Region && typeof(q) != Region 
                if (p in a.pos_prec && q in b.pos_prec)
                    return true
                end
            end 
        end
    end
    #2
    μall = true
    for r1 in a.continuous_prec
        for r2 in b.continuous_prec
            if r1 != r2 
                if !([r1, r2] in μprops) μall = false end
            end
        end
    end
    if !μall return false end 
    #3
    for r1 in a.continuous_prec
        μall = false 
        for d2 in b.pos_prec 
            if [r1, d2] in μprops  μall = true end
        end
        if !μall return false end
    end
    #have to add resolved_intersection mutex conditions at some point
    return false
end


function get_action_mutexes(graph, level)
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

 

