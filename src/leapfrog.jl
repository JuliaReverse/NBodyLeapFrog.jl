using .Bodies: G_year_AU, Body

@inline function acceleration(ra, rb, mb, G) # Get acceleration of celestial body
    d = distance(ra, rb)
    (G*mb/d^3) * (rb - ra)
end

function leapfrog(planets::AbstractVector{Body{T}}; G=G_year_AU, n, dt, keep_history=false) where T
    nplanets = length(planets)
    if keep_history
        v = zeros(V3{T}, (nplanets,n+1))
        r = zeros(V3{T}, (nplanets,n+1))
        v[:,1] .= getfield.(planets, :v)
        r[:,1] .= getfield.(planets, :r)
    else
        v = zeros(V3{T}, nplanets)
        r = zeros(V3{T}, nplanets)
        v[:] .= getfield.(planets, :v)
        r[:] .= getfield.(planets, :r)
    end
    leapfrog!(r, v, planets; G=G, n=n, dt=dt)
end

function leapfrog!(r::AbstractVecOrMat, v::AbstractVecOrMat,
         planets::AbstractVector{Body{T}}; G, n, dt) where T
    nplanets = length(planets)
    a = zeros(V3{T}, nplanets)
    @inbounds for i=1:n
        idx = r isa AbstractVector ? 1 : i
        idx_plus1 = r isa AbstractVector ? 1 : i+1
        for j=1:nplanets
            # compute acceleration
            a[j] = zero(V3{T})
            for k=1:nplanets
                if j != k
                    a[j] += acceleration(r[j,idx], r[k,idx], planets[k].m, G)
                end
            end
        end
        for j=1:nplanets
            v[j,idx_plus1] = (i==1 ? dt/2 : dt)*a[j] + v[j,idx]
            r[j,idx_plus1] = dt*v[j,idx_plus1] + r[j,idx]
        end
    end
    return r, v
end