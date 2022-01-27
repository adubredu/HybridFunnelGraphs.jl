
function get_pose_from_fxn(exp) 
    xs = 0.1:0.1:5
    ys = 0.1:0.1:5
    ls = [] 
    for xg in xs 
        for yg in ys 
            @eval f(xg, yg)=$exp  
            if f(xg, yg) 
                push!(ls, (xg, yg))
            end
        end
    end
    return ls
end
# xg=0.; yg=0.
get_pose_from_fxn(:(2xg+3yg>=0.7&&2xg+3yg<=1.7))