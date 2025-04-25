function flux_update(grid::InitializeGrid, Flux_func::Function, source_terms_func::Function, 
    velocity_func::Function, threshold::Float64, sign_list::Vector{Int})

    N = length(grid.x)

    wave_speed = zeros(N)
    source_term_value = zeros(N)
    flux = zeros(N)
    rhs = zeros(N)

    soln_t = grid.sol_a_t
    shock_sig = grid.shock_sig_t  

    al, ar = compute_interface_states(grid)

    @inbounds for i in grid.ilo:grid.ihi
        source_term_value[i] = 0.5 * (source_terms_func(ar[i-1], grid.xl[i], sign_list[i]) + source_terms_func(al[i], grid.xr[i], sign_list[i]))
        wave_speed[i] = velocity_func(soln_t[i], grid.x[i])
    end

    @inbounds for i in grid.ilo-1:grid.ihi
        flux[i], shock_sig[i] = combine_KT_Riemann_solver(grid.xr[i], al[i], ar[i], grid.dx, Flux_func, velocity_func, threshold)
    end

    @inbounds for i in grid.ilo:grid.ihi  
        rhs[i] = (-(flux[i] - flux[i-1]) / grid.dx + source_term_value[i])
    end

    maxwave = maximum(abs.(wave_speed[3:end-2]))

    return maxwave, rhs, shock_sig
end