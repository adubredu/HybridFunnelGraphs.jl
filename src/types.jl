# x(t+1) = A*x(t) + B*u(t)
struct Dynamics 
    A 
    B 
end

# ϕ₁*x + ϕ₂*y + ϕ₃*z ≥ -þ
# ϕ₁*x + ϕ₂*y + ϕ₃*z ≥ -þ == -ϕ₁*x - ϕ₂*y - ϕ₃*z ≤ þ
struct Ineq 
    ϕ₁::Float64
    ϕ₂::Float64
    ϕ₃::Float64
    þ::Float64
end


struct Region 
    name::String
    r::Vector{Ineq}
    θ::Float64
end


mutable struct DiscreteAction 
    name
    args
    pos_prec 
    neg_prec 
    pos_eff 
    neg_eff      
    function Discrete_Action(name)
        new(name, [], [], [], [], [])
    end
end

mutable struct Funnel 
    name 
    params 
    pos_prec 
    neg_prec 
    continuous_prec 
    dynamics 
    pos_eff 
    neg_eff 
    function Funnel(name)
        new(name, [],[],[],[],[],[],[])
    end

end

mutable struct NoOp  
    name 
    args
    pos_prec 
    neg_prec 
    pos_eff 
    neg_eff 
    function NoOp(prop)
        new(:NoOp, [], [prop], [], [prop], [])
    end
end


mutable struct Graph 
    num_levels::Int64 
    acts 
    μacts
    props 
    μprops 
    leveled 
    initprops 
    goalprops 
    function Graph()
        new(0, Dict(1=>[]), Dict(1=>[]),  Dict(1=>[]), Dict(1=>Dict()), 
        false, [], [])
    end
end