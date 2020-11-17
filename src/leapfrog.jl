using .Bodies: G_year_AU, Body

@inline function acceleration(ra, rb, mb, G) # Get acceleration of celestial body
    d = distance(ra, rb)
    (G*mb/d^3) * (rb - ra)
end

function leapfrog(planets::AbstractVector{Body{T}}; G=G_year_AU, n, dt, keep_history=false) where T
    nplanets = length(planets)
    if keep_history
        a = zeros(V3{T}, (nplanets,n))
        v = zeros(V3{T}, (nplanets,n+1))
        r = zeros(V3{T}, (nplanets,n+1))
        v[:,1] .= getfield.(planets, :v)
        r[:,1] .= getfield.(planets, :r)
    else
        a = zeros(V3{T}, nplanets)
        v = zeros(V3{T}, nplanets)
        r = zeros(V3{T}, nplanets)
        v[:] .= getfield.(planets, :v)
        r[:] .= getfield.(planets, :r)
    end
    leapfrog!(r, v, a, planets; G=G, n=n, dt=dt)
end

function leapfrog!(r::AbstractVecOrMat, v::AbstractVecOrMat, a::AbstractVecOrMat,
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
            end
        else # Use leap frog method to update velocity
            for j=1:nplanets
                v[j,idx_plus1] = dt*a[j,idx] + v[j,idx]
            end
        end

        # update position
        for j=1:nplanets
            r[j,idx_plus1] = dt*v[j,idx_plus1] + r[j,idx]
        end
    end
    return r, v, a
end