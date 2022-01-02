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

mutable struct HybridAction 
    name 
    params 
    pos_prec 
    neg_prec 
    continuous_prec 
    dynamics 
    pos_eff 
    neg_eff 
    function HybridAction(name)
        new(name, [],[],[],[],[],[],[])
    end

end

mutable struct HybridNoOp 
    name 
    params 
    pos_prec 
    neg_prec 
    continuous_prec::Vector{Continuous_condition} 
    dynamics 
    pos_eff 
    neg_eff 
    function HybridNoOp(prop)
        new(:NoOp, [],[prop],[],[],[],[prop],[])
    end

end