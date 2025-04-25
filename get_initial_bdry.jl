using SpecialFunctions, Roots, Polylogarithms

# Effective Mass Function without normalization (raw function)
function Meff_func_raw(y)
    a, b, c, r = 0, rho0, 0, 100 

    # Compute Meff(y)
    exp_term = exp(r * (-r0 + c + y))  # Precompute the exponential term
    log_term = log1p(exp_term)  # Numerically stable log(1 + exp_term)

    return (4 * π / (3 * r^3)) * (
        r^2 * y^2 * (b * r * y + 3 * (a - b) * log_term) +
        6 * (a - b) * r * y * real(polylog(2, -exp_term)) -
        6 * (a - b) * real(polylog(3, -exp_term))
    )
end

# Effective Mass Function with automatic normalization
const Meff_0 = Meff_func_raw(0)  # Normalization factor at y=0
const Meff_5 = Meff_func_raw(5)  # Normalization factor at y=5

function Meff_func(y)
    result = Meff_func_raw(y)
    return M * (result - Meff_0) / Meff_5
end

# Mass Function for OS model
function Meff_sharp(y)
    return y <= r0 ? (4/3) * π * rho0 * y^3 : M
end

# Set the initial condition for the grid
function piecewise_initial(grid::InitializeGrid, mass_function::Function)
    @inbounds for i in eachindex(grid.x)
        Meff = mass_function(grid.x[i]) 

        term1 = 2.0 * Meff * grid.x[i]^(2*l-1) - 4.0 * Meff^2 * zeta^2 * grid.x[i]^(2*l-4)
        grid.sol_a_t[i] = sqrt(max(term1, 0))  # Avoid NaN from negative sqrt
    end
    grid.a0 .= grid.sol_a_t  # More efficient than `copy()`
end

# Define the inner solution function
function X_inner(t::Float64, r::Float64, c1::Float64)
    numerator = r^(1 + l) * (-6.0 * t - 8.0 * r^(1 + l) * zeta^2 * c1)
    denominator = 4.0 * zeta^2 + (3.0 * t + 4.0 * r^(1 + l) * zeta^2 * c1)^2
    return numerator / denominator
end

function c_value(r_left::Float64, mass_function::Function)
    mass = mass_function(r_left)
    X_left = sqrt(max(2.0 * mass * r_left^(2*l-1) - 
                      4.0 * mass^2 * zeta^2 * r_left^(2*l-4), 0)) 
    return find_zero(c -> X_left - X_inner(0.0, r_left, c), -1.2, atol=1e-15, rtol=1e-15)
end

# Define the vacuum solution function
function X_outer(r::Float64, c2::Float64)
    return sqrt(max(-c2 * (c2 - 2 * r^3), 0)) / r^3  # Avoid NaN from sqrt of negatives
end

function c_value_outer(r_right::Float64, mass_function::Function)
    mass = mass_function(r_right)
    X_right = sqrt(max(2.0 * mass * r_right^(2*l-1) - 
                       4.0 * mass^2 * zeta^2 * r_right^(2*l-4), 0))
    return find_zero(c -> X_right - X_outer(r_right, c), 0.8, atol=1e-15, rtol=1e-15)
end

# Fill Boundary Conditions
function fill_BCs(grid::InitializeGrid, t::Float64, c1::Float64)
    @inbounds for i in 1:grid.ng
        grid.sol_a_t[i] = X_inner(t, grid.x[i], c1)
    end
    @inbounds for i in grid.ihi+1:grid.ihi+grid.ng
        grid.sol_a_t[i] = grid.a0[i]
    end
end