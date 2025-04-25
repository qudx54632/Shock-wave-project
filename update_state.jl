function initialize_state(grid::InitializeGrid, C::Float64, r_left::Float64, mass_func::Function)
    piecewise_initial(grid, mass_func)
    
    a_tn = copy(grid.sol_a_t)  # Store previous valid state
    sign_list = ones(Int, length(grid.sol_a_t))
    ca_values = c_value(r_left, mass_func)
    
    dt = C * (grid.dx) / 10  # Initial stable dt
    max_retries = 50  
    flip_count = zeros(Int, length(grid.sol_a_t))

    return a_tn, sign_list, ca_values, dt, max_retries, flip_count
end

function compute_time_step(grid::InitializeGrid, C::Float64, Flux_func::Function, 
                           source_terms_func::Function, velocity_func::Function, 
                           threshold::Float64, sign_list::Vector{Int})
    
    max_wavespeed, k1, _ = flux_update(grid, Flux_func, source_terms_func, velocity_func, threshold, sign_list)
    dt = C * (grid.dx) / max_wavespeed
    
    return dt, k1
end

function correct_sign_and_nan(grid::InitializeGrid, a_tn::Vector{Float64}, sign_list::Vector{Int}, 
    flip_count::Vector{Int}, zeta::Float64, retry_count::Int, dt::Float64)

    nan_detected = false
    sign_changed = false

    @inbounds for i in grid.ilo:grid.ihi
        if isnan(grid.sol_a_t[i])
            if abs(a_tn[i] - 1/(2*zeta)) < 1e-7 || abs(a_tn[i] + 1/(2*zeta)) < 1e-7
                sign_list[i] = -sign_list[i]
                flip_count[i] += 1
                println("Flipped sign at position $i")
                sign_changed = true
            else
                nan_detected = true
                retry_count += 1
                dt = (grid.dx)^2 / retry_count
                break  # Exit early for retry
            end
        end
    end

    return nan_detected, sign_changed, retry_count, dt, sign_list
end

function check_shock_sign_changes(grid::InitializeGrid, sign_list::Vector{Int}, shock_sig::Vector{Bool}, t::Float64)
    shock_pos = findall(!iszero, shock_sig)
    sign_changed = false

    if !isempty(shock_pos)
        if is_shock_moving_right(grid)
            for i in shock_pos
                if sign_list[i] != sign_list[3]
                    println("Sign change at shock $i")
                    sign_list[i] = -sign_list[i]
                    sign_changed = true
                end
            end
        elseif is_shock_moving_left(grid) && t > 0.74
            for i in shock_pos
                if sign_list[i] == sign_list[3]
                    println("Sign change at shock $i")
                    sign_list[i] = -sign_list[i]
                    sign_changed = true
                end
            end
        end
    end

    return sign_changed, sign_list
end


function update_state(grid::InitializeGrid, C::Float64, tstart::Float64, tend::Float64, 
    Flux_func::Function, source_terms_func::Function, velocity_func::Function, mass_func::Function, threshold::Float64)

    # Initialize state variables
    a_tn, sign_list, ca_values, dt, max_retries, flip_count = initialize_state(grid, C, grid.x[grid.ilo], mass_func)

    t = tstart
    step_counter = 0
    slots = 10

    while t < tend
        if isempty(grid.sol_a)
            a_tn = grid.a0
        else
            a_tn = grid.sol_a[end]
        end

        valid_step = false
        retry_count = 0

        while !valid_step
            if retry_count >= max_retries
                error("Too many dt retries ($max_retries) at time $t, terminating.")
            end

            grid.sol_a_t = copy(a_tn)  # Reset state

            # Compute time step and first update
            dt, k1 = compute_time_step(grid, C, Flux_func, source_terms_func, velocity_func, threshold, sign_list)

            if dt + t > tend
                dt = tend - t
            end

            grid.sol_a_t += 0.5 * dt * k1
            fill_BCs(grid, t + 0.5 * dt, ca_values)

            # Compute k2 and update solution
            _, k2, shock_sig = flux_update(grid, Flux_func, source_terms_func, velocity_func, threshold, sign_list)
            grid.sol_a_t .= a_tn .+ dt .* k2

            # Check for NaNs and sign flips
            nan_detected, sign_changed, retry_count, dt = correct_sign_and_nan(grid, a_tn, sign_list, flip_count, zeta, retry_count, dt)

            # Check for shock sign changes
            shock_changed, sign_list = check_shock_sign_changes(grid, sign_list, shock_sig, t)

            if nan_detected
                continue  # Retry with smaller dt
            elseif sign_changed || shock_changed
                valid_step = false  # Repeat with updated sign changes
            else
            # Finalize update and save step
                t += dt
                fill_BCs(grid, t, ca_values)

                if step_counter % 2000 == 0
                    analytic_value = [X_inner(t, grid.x[i], ca_values) for i in 1:slots]
                    # error = maximum(abs.(analytic_value[1:slots] - grid.sol_a_t[1:slots]))
                    println("t = $t")
                end    

                push!(grid.sol_a, copy(grid.sol_a_t))
                push!(grid.t, t)
                push!(grid.shock_sig, copy(shock_sig))
                push!(grid.sign_list, copy(sign_list)) 
                
                step_counter += 1
                valid_step = true
            end
        end
    end

    return grid
end

