# Optimized Kurganov-Tadmor (KT) Solver
@inline function KT_solver(xi::Float64, al::Float64, ar::Float64, Flux_func::Function, velocity_func::Function, dx::Float64)
    ul, ur = velocity_func(al, xi), velocity_func(ar, xi)
    flux_l, flux_r = Flux_func(al, xi), Flux_func(ar, xi)
    
    speed = max(abs(ul), abs(ur))
    flux = 0.5 * (flux_l + flux_r - speed * (ar - al))
    
    shock_sig = ul > ur + 10 * dx  # Detect shocks
    return flux, shock_sig
end

# Optimized Riemann Solver
@inline function Riemann_solver(xi::Float64, al::Float64, ar::Float64, Flux_func::Function, velocity_func::Function, dx::Float64)
    ul, ur = velocity_func(al, xi), velocity_func(ar, xi)
    shock_sig = ul > ur + 10 * dx  # Detect shocks

    if shock_sig
        flux_l, flux_r = Flux_func(al, xi), Flux_func(ar, xi)
        shockspeed = (flux_r - flux_l) / (ar - al)
        as = shockspeed > 0 ? al : shockspeed < 0 ? ar : 0.0  
    else
        as = ur < 0 ? ar : ul > 0 ? al : 0.0
    end

    return Flux_func(as, xi), shock_sig
end

# Optimized Combined KT-Riemann Solver
@inline function combine_KT_Riemann_solver(xi::Float64, al::Float64, ar::Float64, dx::Float64, 
                                           Flux_func::Function, velocity_func::Function, threshold::Float64)
    gradient = abs((al - ar) / dx)

    return gradient > threshold ? 
           Riemann_solver(xi, al, ar, Flux_func, velocity_func, dx) :
           KT_solver(xi, al, ar, Flux_func, velocity_func, dx)
end