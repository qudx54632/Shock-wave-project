using Roots

# Optimized ODEs_Original Function
function ODEs_Original(Xsol::Float64, r::Float64, sign::Int)
    dr_dt = -r^(-l) * Xsol
    term1 = r^(2 + (2 * l)) - 4 * zeta^2 * Xsol^2

    if term1 < 0
        return dr_dt, NaN  # Return NaN for invalid term1
    end

    dXsol_dt = (2 - l) * Xsol^2 / r^(l + 1) + (3 * (-r^(l + 1) + sign * sqrt(max(term1, 0)))) / (4 * zeta^2)
    return -dr_dt, -dXsol_dt
end

# Optimized RK4 step function
@inline function rk4_step_Original(f::Function, Bsol::Float64, x::Float64, dt::Float64, sign::Int)
    dx1, dBsol1 = f(Bsol, x, sign)
    dx2, dBsol2 = f(Bsol + 0.5 * dt * dBsol1, x + 0.5 * dt * dx1, sign)
    dx3, dBsol3 = f(Bsol + 0.5 * dt * dBsol2, x + 0.5 * dt * dx2, sign)
    dx4, dBsol4 = f(Bsol + dt * dBsol3, x + dt * dx3, sign)

    x_next = x + (dt / 6) * (dx1 + 2*dx2 + 2*dx3 + dx4)
    Bsol_next = Bsol + (dt / 6) * (dBsol1 + 2*dBsol2 + 2*dBsol3 + dBsol4)

    return x_next, Bsol_next
end

# Optimized Characteristic Solver
function Characteristic_Original_solver(grid::InitializeGrid, xini_index::Int, t_range::Vector{Float64})
    Bsol, x = grid.a0[xini_index], grid.x[xini_index]  # Initial conditions

    N = length(t_range)
    tvalue = Vector{Float64}(undef, N) 
    x_vals = Vector{Float64}(undef, N)
    Bsol_vals = Vector{Float64}(undef, N)

    sign = 1
    dt = t_range[2] - t_range[1]  # Time step

    @inbounds for i in 1:N
        valid_step = false
        while !valid_step
            # Perform an RK4 step
            x_next, Bsol_next = rk4_step_Original(ODEs_Original, Bsol, x, dt, sign)

            # Check for invalid values (NaN)
            if isnan(Bsol_next)
                sign = -sign  # Flip sign
            else
                # Store values and update for next step
                tvalue[i] = t_range[i]
                Bsol_vals[i] = Bsol
                x_vals[i] = x
                Bsol, x = Bsol_next, x_next
                valid_step = true
            end
        end
    end 

    return tvalue, x_vals, Bsol_vals
end
