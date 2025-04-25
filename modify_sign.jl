using Statistics

function is_shock_moving_right(grid)
    shock_pos = findall(!iszero, grid.shock_sig[end])  # Current shock positions

    # Return false if no shocks exist
    if isempty(shock_pos) || length(grid.shock_sig) < 2
        return false
    end

    prev_shock_pos = findall(!iszero, grid.shock_sig[end-1])  # Previous shock positions

    # Return false if no previous shocks exist
    if isempty(prev_shock_pos)
        return false
    end

    # Compute mean positions of shocks
    mean_curr = mean(shock_pos)
    mean_prev = mean(prev_shock_pos)

    return mean_curr > mean_prev  # True if shock moved right
end

function is_shock_moving_left(grid)
    shock_pos = findall(!iszero, grid.shock_sig[end])  # Current shock positions

    # Return false if no shocks exist
    if isempty(shock_pos) || length(grid.shock_sig) < 2
        return false
    end

    prev_shock_pos = findall(!iszero, grid.shock_sig[end-1])  # Previous shock positions

    # Return false if no previous shocks exist
    if isempty(prev_shock_pos)
        return false
    end

    # Compute mean positions of shocks
    mean_curr = mean(shock_pos)
    mean_prev = mean(prev_shock_pos)

    return mean_curr < mean_prev  # True if shock moved left
end