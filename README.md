# HybridFunnelGraphs.jl
Julia Package for generating Hybrid Funnel Graphs with discrete and continuous variables and constraints.

## Installation
1. Open your Julia REPL by typing  `julia` in your terminal.
2. Press `]` on your keyboard to enter the package manager
3. Enter command `add https://github.com/adubredu/HybridFunnelGraphs.jl` and press 
`Enter` on your keyboard to install this package.
4. Press the `Backspace` key on your keyboard to return to the REPL

## Usage
A Hybrid Funnel Graph is generated from [.hpd](https://github.com/adubredu/HPD.jl) files that describe an optimal constrained task planning problem. Two `.hpd` files are required to generate a Hybrid Funnel Graph. The first `.hpd` file is one which contains a description of the planning domain and contains specifications of action primitives. The second `.hpd` file required is one that describes the initial and final states of continuous and symbolic propositions as well as any external constraints and an objective function.

Here's an example usage of this package:
```julia
using HybridFunnelGraphs

# gets path to .hpd files
path = joinpath(dirname(pathof(HybridFunnelGraphs)), "..", "test/hpd")
domain_path = joinpath(path, "domain.hpd")
problem_path = joinpath(path, "problem.hpd")

#creates a Hybrid Funnel Graph struct from .hpd files
graph = create_funnel_graph(domain_path, problem_path)
```

A Hybrid Funnel Graph (HFG) has the following key attributes:

1. `num_levels`: Integer-valued number of state levels in the HFG
2. `acts`: A dictionary whose keys are level indices and values are lists of funnels in the corresponding level.
3. `props`: A dictionary whose keys are level indices and values are lists of dictionaries of discrete and continuous propositions in the corresponding level.
4. `initprops`: A dictionary of discrete and continuous initial propositions.
5. `goalprops`: A dictionary of discrete and continuous goal propositions.