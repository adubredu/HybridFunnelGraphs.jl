# x(t+1) = A*x(t) + B*u(t)
struct Dynamics 
    A 
    B 
end

# || ϕ₁ - ϕ₂ || ≤ þ
struct Continuous_condition 
    ϕ₁
    ϕ₂
    þ
end


struct Region 
    xᵣ::Vector{Float64}
    yᵣ::Vector{Float64}
    θ::Float64
end


mutable struct Discrete_Action 
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

