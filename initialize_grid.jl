mutable struct InitializeGrid
    xmin::Float64                 # Minimum x value
    xmax::Float64                 # Maximum x value
    nx::Int                       # Number of intervals
    ng::Int                       # Number of ghost cells
    dx::Float64                   # Interval length
    ilo::Int                      # Index of the first physical cell
    ihi::Int                      # Index of the last physical cell
    x::Vector{Float64}            # Cell-center coordinates
    xl::Vector{Float64}           # Left-edge coordinates
    xr::Vector{Float64}           # Right-edge coordinates
    a0::Vector{Float64}           # Initial values
    sol_a_t::Vector{Float64}      # Solution at the current time step
    sol_a::Vector{Vector{Float64}} # Solution array over time
    shock_sig_t::Vector{Bool}      # Shock signal at current time
    shock_sig::Vector{Vector{Bool}} # Shock signal array over time
    sign_list::Vector{Vector{Int}} # Sign change list over time
    t::Vector{Float64}            # Time coordinates

    function InitializeGrid(nx::Int, ng::Int, xmin::Float64, xmax::Float64)
        dx = (xmax - xmin) / (nx + 2 * ng)
        ilo, ihi = ng + 1, ng + nx

        # Efficient array initialization
        x = range(xmin, step=dx, length=nx + 2 * ng) |> collect
        xl = x .- 0.5 * dx
        xr = x .+ 0.5 * dx

        # Preallocating memory
        size_total = nx + 2 * ng
        a0 = zeros(size_total)
        sol_a_t = zeros(size_total)
        shock_sig_t = falses(size_total)  # Efficient way to initialize a Bool array

        # Use empty vectors to store solutions over time
        sol_a = Vector{Vector{Float64}}()
        shock_sig = Vector{Vector{Bool}}()
        sign_list = Vector{Vector{Int}}()
        t = Float64[]

        # Construct the object
        return new(xmin, xmax, nx, ng, dx, ilo, ihi, x, xl, xr, a0, sol_a_t, sol_a, shock_sig_t, shock_sig, sign_list, t)
    end
end