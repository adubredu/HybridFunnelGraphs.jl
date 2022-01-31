using Revise 
using HybridFunnelGraphs

graph = create_funnel_graph("test/hpd/domain.hpd", "test/hpd/problem.hpd"; max_levels=25)
graph.num_levels