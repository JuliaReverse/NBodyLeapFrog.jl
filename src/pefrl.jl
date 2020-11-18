# Forest-Ruth and Position Extended Forest-Ruth Like algorithm
# Leapfrog is 2nd order, FR and PEFRL are 4th order.

function fr_step!(r::AbstractVecOrMat{V3{T}}, v::AbstractVecOrMat{V3{T}},
         m::AbstractVector; G, dt) where T
    nplanets = length(m)
    a = zeros(V3{T}, nplanets)
    θ = 1/(2-2^(1/3))
    @inbounds begin
        for j=1:nplanets
            r[j] += θ/2*dt*v[j]
        end
        update_acceleration!(a, r, m, G)
        for j=1:nplanets
            v[j] += θ*dt*a[j]
            r[j] += (1-θ)/2*dt*v[j]
        end
        update_acceleration!(a, r, m, G)
        for j=1:nplanets
            v[j] += (1-2θ)*dt*a[j]
            r[j] += (1-θ)/2*dt*v[j]
        end
        update_acceleration!(a, r, m, G)
        for j=1:nplanets
            v[j] += θ*dt*a[j]
            r[j] += θ/2*dt*v[j]
        end
    end
    return r, v
end

function pefrl_step!(r::AbstractVecOrMat{V3{T}}, v::AbstractVecOrMat{V3{T}},
         m::AbstractVector; G, dt) where T
    nplanets = length(m)
    a = zeros(V3{T}, nplanets)
    ξ = +0.1786178958448091E+00
    λ = -0.2123418310626054E+00
    χ = -0.6626458266981849E-01
    @inbounds begin
        for j=1:nplanets
            r[j] = r[j]
            v[j] = v[j]
            r[j] += ξ*dt*v[j]
        end
        update_acceleration!(a, r, m, G)
        for j=1:nplanets
            v[j] += (1-2λ)/2*dt*a[j]
            r[j] += χ*dt*v[j]
        end
        update_acceleration!(a, r, m, G)
        for j=1:nplanets
            v[j] += λ*dt*a[j]
            r[j] += (1-2(ξ+χ))*dt*v[j]
        end
        update_acceleration!(a, r, m, G)
        for j=1:nplanets
            v[j] += λ*dt*a[j]
            r[j] += χ*dt*v[j]
        end
        update_acceleration!(a, r, m, G)
        for j=1:nplanets
            v[j] += (1-2λ)/2*dt*a[j]
            r[j] += ξ*dt*v[j]
        end
    end
    return r, v
end

function fr!(r::AbstractVecOrMat{V3{T}}, v::AbstractVecOrMat{V3{T}},
         m::AbstractVector; G=G_year_AU, n, dt) where T
    leapfrog_looper!(FR,r, v, m; G=G, n=n, dt=dt)
end

function pefrl!(r::AbstractVecOrMat{V3{T}}, v::AbstractVecOrMat{V3{T}},
         m::AbstractVector; G=G_year_AU, n, dt) where T
    leapfrog_looper!(PEFRL,r, v, m; G=G, n=n, dt=dt)
end

function fr(planets::AbstractVector{Body{T}}; G=G_year_AU, n, dt, keep_history=false) where {T}
    leapfrog(FR, planets; G=G, n=n, dt=dt, keep_history=keep_history)
end

function pefrl(planets::AbstractVector{Body{T}}; G=G_year_AU, n, dt, keep_history=false) where {T}
    leapfrog(PEFRL, planets; G=G, n=n, dt=dt, keep_history=keep_history)
end

steper(::Type{FR}) = fr_step!
steper(::Type{PEFRL}) = pefrl_step!