# Optimized Minmod Function
@inline function minmod(a0::Float64, a1::Float64, a2::Float64)
    diff_left = a1 - a0
    diff_right = a2 - a1
    return sign(diff_left) == sign(diff_right) ? min(abs(diff_left), abs(diff_right)) * sign(diff_left) : 0.0
end

# Optimized Compute Interface States Function
function compute_interface_states(grid::InitializeGrid)
    a = grid.sol_a_t
    al = similar(a)  # Preallocate array without copying data
    ar = similar(a)

    @inbounds for i in grid.ilo-1:grid.ihi
        dal = minmod(a[i-1], a[i], a[i+1])
        dar = minmod(a[i], a[i+1], a[i+2])
        al[i] = a[i] + 0.5 * dal
        ar[i] = a[i+1] - 0.5 * dar
    end
    
    return al, ar
end