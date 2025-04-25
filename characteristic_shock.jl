function find_nearest_index(value::Float64, array::Vector{Float64})::Int
    return argmin(abs.(array .- value))  # Find index of closest value
end

function dr_dt(r::Float64, t_idx::Int, grid::InitializeGrid)
    x_idx = find_nearest_index(r, grid.x)  # Get nearest x index
    soln_value = @inbounds grid.sol_a[t_idx][x_idx]  # Fast lookup from grid.sol_a
    sign_value = @inbounds grid.sign_list[t_idx][x_idx]
    drdt_value = -r^(-l) * soln_value
    return drdt_value, sign_value  # Modify equation as needed
end

@inline function euler_step(f::Function, r::Float64, t_idx::Int, dt::Float64, grid::InitializeGrid)
    drdt_value, sign_value = f(r, t_idx, grid)
    return r + dt * drdt_value, sign_value  # Single-line Euler step
end

function solve_r_ode(grid::InitializeGrid, r0_idx::Int)
    num_t = length(grid.t)
    r_values = zeros(num_t)  # Preallocate results
    sign_values = ones(num_t)
    r_values[1] = grid.x[r0_idx]  # Initial condition

    @inbounds for i in 1:num_t-1
        dt = grid.t[i+1] - grid.t[i]  # Compute time step
        r_values[i+1], sign_values[i+1] = euler_step(dr_dt, r_values[i], i, dt, grid)  # Update using Euler method
    end

    return r_values, sign_values
end

using Polynomials  # Required for polynomial fitting
function shock_surface(grid::InitializeGrid, tstart::Int, tend::Int)
    N = tend - tstart + 1  # Number of time steps
    xoutlist = Vector{Float64}(undef, N)  # Preallocate
    toutlist = Vector{Float64}(undef, N)  # Preallocate

    count = 0  # Track valid entries

    @inbounds for i in tstart:tend
        # Find shock position indices (nonzero elements)
        shock_indexall = findall(!iszero, grid.shock_sig[i])

        # Ensure we have a valid shock position
        shock_index = isempty(shock_indexall) ? nothing : last(shock_indexall)

        # Store results if a shock was found
        if shock_index !== nothing
            count += 1
            xoutlist[count] = grid.x[shock_index]
            toutlist[count] = grid.t[i]
        end
    end

    # Resize vectors to valid count (to remove unused preallocated elements)
    xoutlist = xoutlist[1:count]
    toutlist = toutlist[1:count]

    if isempty(toutlist) || isempty(xoutlist)
        error("No shock positions were detected in the given grids.")
    end

    # Fit a 3rd-degree polynomial (more stable than degree 10)
    x_t_shock = fit(Polynomial, toutlist, xoutlist, 10)
    errors = maximum([abs(x_t_shock(toutlist[i]) - xoutlist[i])  for i in 1:count])

    return x_t_shock, minimum(toutlist), maximum(toutlist), errors
end