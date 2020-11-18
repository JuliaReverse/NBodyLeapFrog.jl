export i_leapfrog!, i_leapfrog_reuse!

@i @inline function acceleration_reuse(y!::V3{T}, ra::V3{T}, rb::V3{T}, mb::Real, G) where T
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

@i @inline function :(+=)(acceleration)(y!::V3{T}, ra::V3{T}, rb::V3{T}, mb::Real, G) where T
    @routine @invcheckoff begin
        @zeros T d anc1 anc2 anc3 anc4
        rc ← zero(V3{T})
        d += sqdistance(ra, rb)
        anc1 += sqrt(d)
        anc2 += anc1 ^ 3
        anc3 += G * mb
        anc4 += anc3 / anc2
        rc += rb - ra
    end
    y! += anc4 * rc
    ~@routine
end

@i function i_update_acceleration!(a, r, m, G)
    @safe @assert length(a) == length(m) == length(r)
    @inbounds for j=1:length(a)
        for k=1:length(a)
            if j != k
                a[j] += acceleration(r[j], r[k], m[k], G)
            end
        end
    end
end

@i function i_leapfrog!(r::AbstractVector{V3{T}}, v::AbstractVector{V3{T}}, planets::AbstractVector{Body{T}}; G=G_year_AU, n, dt) where T
    @routine @invcheckoff begin
        nplanets ← length(planets)
        m ← zeros(T, nplanets)
        for i=1:nplanets
            m[i] += planets[i].m
        end
    end
    @invcheckoff @inbounds for i=1:nplanets
        v[i] += planets[i].v
        r[i] += planets[i].r
    end
    @invcheckoff for i=1:n
        i_leapfrog_step!(r, v, m; G=G, dt=dt)
    end
    ~@routine
end

@i function i_leapfrog_step!(r::AbstractVector{V3{T}}, v::AbstractVector{V3{T}}, m::AbstractVector; G=G_year_AU, dt) where T
    @routine @invcheckoff begin
        nplanets ← length(m)
        halfdt ← zero(dt)
        halfdt += dt/2
    end
    @inbounds @invcheckoff begin
        for j=1:nplanets
            r[j] += halfdt*v[j]
        end
        # compute acceleration and update velocity
        @routine begin
            a ← zeros(V3{T}, nplanets)
            i_update_acceleration!(a, r, m, G)
        end
        for j=1:nplanets
            v[j] += dt*a[j]
        end
        ~@routine

        for j=1:nplanets
            r[j] += halfdt*v[j]
        end
    end
    ~@routine
end

@i function i_leapfrog_reuse!(r::AbstractVector{V3{T}}, v::AbstractVector{V3{T}}, planets::AbstractVector{Body{T}}; G=G_year_AU, n, dt) where T
    @routine @invcheckoff begin
        nplanets ← length(planets)
        halfdt ← zero(dt)
        halfdt += dt/2
        m ← zeros(T, nplanets)
        a ← zeros(V3{T}, nplanets)
    end
    for i=1:nplanets
        m[i] += planets[i].m
    end
    for i=1:nplanets
        v[i] += planets[i].v
        r[i] += planets[i].r
    end
    @inbounds @invcheckoff for i=1:n
        # update position
        for j=1:nplanets
            r[j] += halfdt*v[j]
        end
        # compute acceleration
        @routine for j=1:nplanets
            for k=1:nplanets
                if j != k
                    acceleration_reuse(a[j], r[j], r[k], m[k], G)
                end
            end
        end

        # update velocity
        for j=1:nplanets
            v[j] += dt*a[j]
        end
        ~@routine

        # update position
        for j=1:nplanets
            r[j] += halfdt*v[j]
        end
    end
    ~@routine
end