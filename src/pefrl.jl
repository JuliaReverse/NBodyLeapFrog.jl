# Forest-Ruth and Position Extended Forest-Ruth Like algorithm
# Leapfrog is 2nd order, FR and PEFRL are 4th order.

function fr!(r::AbstractVecOrMat{V3{T}}, v::AbstractVecOrMat{V3{T}},
         m::AbstractVector; G, n, dt) where T
    nplanets = length(m)
    a = zeros(V3{T}, nplanets)
    θ = 1/(2-2^(1/3))
    @inbounds for i=1:n
        idx = r isa AbstractVector ? 1 : i
        idx_plus1 = r isa AbstractVector ? 1 : i+1
        for j=1:nplanets
            r[j,idx_plus1] = r[j,idx]
            v[j,idx_plus1] = v[j,idx]
            r[j,idx_plus1] += θ/2*dt*v[j,idx_plus1]
        end
        update_acceleration!(a, _view(r,:,idx_plus1), m, G)
        for j=1:nplanets
            v[j,idx_plus1] += θ*dt*a[j]
            r[j,idx_plus1] += (1-θ)/2*dt*v[j,idx_plus1]
        end
        update_acceleration!(a, _view(r,:,idx_plus1), m, G)
        for j=1:nplanets
            v[j,idx_plus1] += (1-2θ)*dt*a[j]
            r[j,idx_plus1] += (1-θ)/2*dt*v[j,idx_plus1]
        end
        update_acceleration!(a, _view(r,:,idx_plus1), m, G)
        for j=1:nplanets
            v[j,idx_plus1] += θ*dt*a[j]
            r[j,idx_plus1] += θ/2*dt*v[j,idx_plus1]
        end
    end
    return r, v
end

function pefrl!(r::AbstractVecOrMat{V3{T}}, v::AbstractVecOrMat{V3{T}},
         m::AbstractVector; G, n, dt) where T
    nplanets = length(m)
    a = zeros(V3{T}, nplanets)
    ξ = +0.1786178958448091E+00
    λ = -0.2123418310626054E+00
    χ = -0.6626458266981849E-01
    @inbounds for i=1:n
        idx = r isa AbstractVector ? 1 : i
        idx_plus1 = r isa AbstractVector ? 1 : i+1
        for j=1:nplanets
            r[j,idx_plus1] = r[j,idx]
            v[j,idx_plus1] = v[j,idx]
            r[j,idx_plus1] += ξ*dt*v[j,idx_plus1]
        end
        update_acceleration!(a, _view(r,:,idx_plus1), m, G)
        for j=1:nplanets
            v[j,idx_plus1] += (1-2λ)/2*dt*a[j]
            r[j,idx_plus1] += χ*dt*v[j,idx_plus1]
        end
        update_acceleration!(a, _view(r,:,idx_plus1), m, G)
        for j=1:nplanets
            v[j,idx_plus1] += λ*dt*a[j]
            r[j,idx_plus1] += (1-2(ξ+χ))*dt*v[j,idx_plus1]
        end
        update_acceleration!(a, _view(r,:,idx_plus1), m, G)
        for j=1:nplanets
            v[j,idx_plus1] += λ*dt*a[j]
            r[j,idx_plus1] += χ*dt*v[j,idx_plus1]
        end
        update_acceleration!(a, _view(r,:,idx_plus1), m, G)
        for j=1:nplanets
            v[j,idx_plus1] += (1-2λ)/2*dt*a[j]
            r[j,idx_plus1] += ξ*dt*v[j,idx_plus1]
        end
    end
    return r, v
end