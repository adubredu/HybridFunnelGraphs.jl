using Revise 
using HybridFunnelGraphs

graph = create_funnel_graph("test/pddl/domain.pddl", "test/pddl/problem.pddl"; max_levels=60)
graph.num_levels