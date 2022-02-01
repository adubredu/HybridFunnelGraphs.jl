function orig_eval(exp::Expr) 
    xs = 0.0:0.1:5
    ys = 0.0:0.1:5 
    ls = []
    fun = @eval (xg, yg) -> $exp 
    # fun = eval(exp) 
    for x in xs 
        for y in ys  
            if Base.invokelatest(fun, x, y) 
                push!(ls, (x, y))
            end
        end
    end 
    return ls 
end

function get_pose_from_fxn!(exp::Expr, ls::Vector{Any}) 
    xs = 0.0:0.1:5
    ys = 0.0:0.1:5 
    fun = @eval (xg, yg) -> $exp
    for x in xs 
        for y in ys  
            if Base.invokelatest(fun, x, y) 
                push!(ls, (x, y))
            end
        end
    end  
end

function get_poses_from_fxn(exps::Vector{Expr})
    ls = []
    for exp in exps 
        get_pose_from_fxn!(exp, ls) 
    end
    return collect(Set(ls)) 
end

function and_get_pose_from_fxn(exps::Vector{Expr}) 
    xs = 0.0:0.1:5
    ys = 0.0:0.1:5
    ls = [] 
    for xg in xs 
        for yg in ys 
            sat = true
            for exp in exps
                @eval f(xg, yg)=$exp   
                sat = sat && f(xg, yg)  
            end
            if sat push!(ls, (xg, yg)) end
        end
    end
    return ls
end

function or_get_pose_from_fxn(exps::Vector{Expr}) 
    xs = 0.0:0.1:5
    ys = 0.0:0.1:5
    ls = [] 
    for exp in exps
        for xg in xs 
            for yg in ys 
                @eval f(xg, yg)=$exp  
                if f(xg, yg) 
                    push!(ls, (xg, yg))
                end
            end
        end
        println("length ls ",length(ls))
    end
    return ls
end
# xg=0.; yg=0.
get_poses_from_fxn([:(0.0<=xg<=0.05&&0.0<=yg<=1.5) , :(2xg+3yg>=0.7&&2xg+3yg<=1.7)] )

# orig_eval(:(0.0<=xg<=0.05&&0.0<=yg<=1.5))