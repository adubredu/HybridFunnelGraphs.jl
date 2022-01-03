using Revise 
using PDDL 

domain = load_domain("test/pddl/domain.pddl")
problem = load_problem("test/pddl/problem.pddl")

inits = collect(initstate(domain, problem).facts)