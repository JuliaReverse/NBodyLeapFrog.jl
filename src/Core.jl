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
                v[j,i+1] = dt*a[j,i] + v[j,i]
            end
        else # Use leap frog method to update velocity
            for j=1:nplanets
                v[j,i+1] = 2dt*a[j,i] + v[j,i-1]
            end
        end

        # update position
        if i ==1
            # Get position with Euler method for i==0
            for j=1:nplanets
                r[j,i+1] = dt*v[j,i] + r[j,i]
            end
        else
            # Update position with leap frog method
            for j=1:nplanets
                r[j,i+1] = 2dt*v[j,i] + r[j,i-1]
            end
        end
    end
    return r, v, a
end

function fast_leapfrog(planets::AbstractVector{Body{T}}; G=G_year_AU, n, dt) where T
    nplanets = length(planets)
    a = zeros(V3{T}, nplanets)
    v = zeros(V3{T}, nplanets,2)
    r = zeros(V3{T}, nplanets,2)
    v[:,1] .= getfield.(planets, :v)
    r[:,1] .= getfield.(planets, :r)
    @inbounds for i=1:n
        # compute acceleration
        for j=1:nplanets
            a[j] = zero(V3{T})
            for k=1:nplanets
                if j != k
                    a[j] += acceleration(r[j,mod1(i,2)], r[k,mod1(i,2)], planets[k].m, G)
                end
            end
        end

        # update velocity
        if i==1
            for j=1:nplanets
                v[j,2] = dt*a[j] + v[j,1]
            end
        else # Use leap frog method to update velocity
            for j=1:nplanets
                v[j,mod1(i+1,2)] += 2dt*a[j]
            end
        end

        # update position
        if i ==1
            # Get position with Euler method for i==0
            for j=1:nplanets
                r[j,2] += dt*v[j,1] + r[j,1]
            end
        else
            # Update position with leap frog method
            for j=1:nplanets
                r[j,mod1(i+1,2)] += 2dt*v[j,mod1(i,2)]
            end
        end
    end
    return r, v, a
end
