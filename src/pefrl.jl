# Forest-Ruth and Position Extended Forest-Ruth Like algorithm
# Leapfrog is 2nd order, FR and PEFRL are 4th order.

function fr!(r::AbstractVecOrMat, v::AbstractVecOrMat, a::AbstractVecOrMat,
         planets::AbstractVector{Body{T}}; G, n, dt) where T
    nplanets = length(planets)
    @inbounds for i=1:n
        idx = r isa AbstractVector ? 1 : i
        idx_plus1 = r isa AbstractVector ? 1 : i+1
        # compute acceleration
        for j=1:nplanets
            a[j,idx] = zero(V3{T})
            for k=1:nplanets
                if j != k
                    a[j,idx] += acceleration(r[j,idx], r[k,idx], planets[k].m, G)
                end
            end
        end

        # update velocity
        if i==1
            for j=1:nplanets
                v[j,idx_plus1] = dt/2*a[j,idx] + v[j,idx]
                r[j,idx_plus1] = dt*v[j,idx_plus1] + r[j,idx]
            end
        else # Use leap frog method to update velocity
            for j=1:nplanets
                r[j,idx_plus1] = r[j,idx]
                v[j,idx_plus1] = v[j,idx]
                r[j,idx_plus1] += θ*dt/2*v[j,idx_plus1]
                v[j,idx_plus1] += θ*dt*a[j]
            end
        end
    end
    return r, v, a
end