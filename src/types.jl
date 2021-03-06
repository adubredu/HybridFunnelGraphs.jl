# x(t+1) = A*x(t) + B*u(t)
struct Dynamics 
    A 
    B
    vx_range 
    vy_range 
    d
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
    name 
    r::Vector{Ineq}
    θ::Float64
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
    end_region 
    is_continuous
    cont_ind 
    function Funnel(name)
        new(name, [],[],[],[],[],[],[],Dict(),true, 1)
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
    indexes #[object_index, place_pose_index]
    safe_poses
    has_placement_constraint
    function Graph()
        new(0, Dict(1=>[]), Dict(1=>[]), Dict(1=>Dict()),  Dict(1=>[]), false, Dict(), Dict(), [1,1],[],true)
    end
end