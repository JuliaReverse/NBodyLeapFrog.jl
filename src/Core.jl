using .Bodies: G_year_AU, Body

@inline function acceleration(ra, rb, mb, G) # Get acceleration of celestial body
    d = distance(ra, rb)
    (G*mb/d^3) * (rb - ra)
end

function leapfrog(planets::AbstractVector{Body{T}}; G=G_year_AU, n, dt) where T
    nplanets = length(planets)
    a = zeros(V3{T}, (nplanets,n))
    v = zeros(V3{T}, (nplanets,n+1))
    r = zeros(V3{T}, (nplanets,n+1))
    v[:,1] .= getfield.(planets, :v)
    r[:,1] .= getfield.(planets, :r)
    @inbounds for i=1:n
        # compute acceleration
        for j=1:nplanets
            for k=1:nplanets
                if j != k
                    a[j,i] += acceleration(r[j,i], r[k,i], planets[k].m, G)
                end
            end
        end

        # update velocity
        if i==1
            for j=1:nplanets
                v[j,i+1] = dt/2*a[j,i] + v[j,i]
            end
        else # Use leap frog method to update velocity
            for j=1:nplanets
                v[j,i+1] = dt*a[j,i] + v[j,i]
            end
        end

        # update position
        for j=1:nplanets
            r[j,i+1] = dt*v[j,i+1] + r[j,i]
        end
    end
    return r, v, a
end

function fast_leapfrog(planets::AbstractVector{Body{T}}; G=G_year_AU, n, dt) where T
    nplanets = length(planets)
    a = zeros(V3{T}, nplanets)
    v = getfield.(planets, :v)
    r = getfield.(planets, :r)
    @inbounds for i=1:n
        # compute acceleration
        for j=1:nplanets
            a[j] = zero(V3{T})
            for k=1:nplanets
                if j != k
                    a[j] += acceleration(r[j], r[k], planets[k].m, G)
                end
            end
        end

        # update velocity
        if i==1
            for j=1:nplanets
                v[j] += dt/2*a[j]
            end
        else # Use leap frog method to update velocity
            for j=1:nplanets
                v[j] += dt*a[j]
            end
        end

        # update position
        for j=1:nplanets
            r[j] += dt*v[j]
        end
    end
    return r, v, a
end
