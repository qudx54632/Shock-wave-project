function compute_flux(u::Float64, r::Float64)
    return -u^2 / (2.0 * r^l)  
end

function compute_source(u::Float64, r::Float64, sign::Int)
    term1 = r^(2 + (2 * l)) - 4 * zeta^2 * u^2

    if term1 <= 0.0    
        return NaN 
    else 
        return -u^2 * (l - 4) / (2 * r^(1 + l)) + 3.0 * (-r^(l + 1) + sign * sqrt(max(term1, 0))) / (4.0 * zeta^2)
    end
end

function compute_velocity(u::Float64, r::Float64)
    return -u * r^(-l)  
end