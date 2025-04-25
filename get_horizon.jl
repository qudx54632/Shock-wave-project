using Roots

function find_all_roots(grid::InitializeGrid, t_index::Int, sig::Int)
    Num = length(grid.x)
    roots = Float64[]
    
    # Define the function at discrete points
    f(x_idx) = grid.x[x_idx] * grid.sol_a[t_index][x_idx] + 1*sig

    for x_index in 1:Num-1
        f1 = f(x_index)
        f2 = f(x_index + 1)

        # Check for a sign change between adjacent points
        if f1 * f2 < 0
            # Linear interpolation to estimate the root location
            x1 = grid.x[x_index]
            x2 = grid.x[x_index + 1]
            root = x1 - f1 * (x2 - x1) / (f2 - f1)
            push!(roots, root)
        end
    end

    return roots
end

function find_all_roots_negative(Xvalue::Vector{Float64}, xvalues::Vector{Float64}, tvalues::Vector{Float64})
    Num = length(Xvalue)
    roots = Float64[]
    rootstime = Float64[]
    
    # Define the function at discrete points
    f(x_idx) = xvalues[x_idx] * Xvalue[x_idx] - 1

    for x_index in 1:Num-1
        f1 = f(x_index)
        f2 = f(x_index + 1)

        # Check for a sign change between adjacent points
        if f1 * f2 < 0
            # Linear interpolation to estimate the root location
            x1 = xvalues[x_index]
            x2 = xvalues[x_index + 1]
            root = x1 - f1 * (x2 - x1) / (f2 - f1)
            push!(roots, root)
            push!(rootstime, (tvalues[x_index]+tvalues[x_index+1])/2)
        end
    end

    return roots, rootstime
end
