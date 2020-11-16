export i_leapfrog
using .Bodies: G_year_AU, Body

@i @inline function :(+=)(acceleration)(y!::V3{T}, ra::V3{T}, rb::V3{T}, mb::Real, G) where T
    @routine @invcheckoff begin
        @zeros T d anc1 anc2 anc3 anc4
        d += sqdistance(ra, rb)
        anc1 += sqrt(d)
        anc2 += anc1 ^ 3
        anc3 += G * mb
        anc4 += anc3 / anc2
        rb -= ra
    end
    y! += anc4 * rb
    ~@routine
end

@i function i_leapfrog(r::AbstractMatrix{V3{T}}, v::AbstractMatrix{V3{T}}, planets::AbstractVector{Body{T}}; G=G_year_AU, n, dt) where T
    @routine @invcheckoff begin
        nplanets ← length(planets)
        twodt ← zero(dt)
        twodt += 2 * dt
        m ← zeros(T, nplanets)
        a ← zeros(V3{T}, nplanets)
        for i=1:nplanets
            v[i,1] += planets[i].v
            r[i,1] += planets[i].r
            m[i] += planets[i].m
        end
    end
    @inbounds @invcheckoff for i=1:n
        # compute acceleration
        @routine for j=1:nplanets
            for k=1:nplanets
                if j != k
                    a[j] += acceleration(r[j,mod1(i,2)], r[k,mod1(i,2)], m[k], G)
                end
            end
        end

        # update velocity
        if i==1
            for j=1:nplanets
                v[j,2] += v[j,1]
                v[j,2] += dt*a[j]
            end
        else # Use leap frog method to update velocity
            for j=1:nplanets
                v[j,mod1(i+1,2)] += twodt*a[j]
            end
        end
        ~@routine

        # update position
        if i ==1
            # Get position with Euler method for i==0
            for j=1:nplanets
                r[j,2] += r[j,1]
                r[j,2] += dt*v[j,1]
            end
        else
            # Update position with leap frog method
            for j=1:nplanets
                r[j,mod1(i+1,2)] += twodt*v[j,mod1(i,2)]
            end
        end
    end
    ~@routine
end