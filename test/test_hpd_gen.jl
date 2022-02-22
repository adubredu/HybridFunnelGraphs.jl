using Revise 
using HybridFunnelGraphs

# graph = create_funnel_graph("test/hpd/domain.hpd", "test/hpd/problem.hpd"; max_levels=170)
graph = create_funnel_graph("test/hpd/dom2.hpd", "test/hpd/prob2.hpd"; max_levels=300, has_placement_constraint=false)
graph.num_levels